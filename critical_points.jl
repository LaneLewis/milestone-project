#using Pkg
#Pkg.add("NLsolve")
module Stability
using NLsolve
using LinearAlgebra
using Plots
function rectified_power_rate(voltage_vec,voltage_threshold_vec,k,power)
    rectified_voltage = map(x -> max(x,0),voltage_vec.-voltage_threshold_vec)
    return k*rectified_voltage.^power
end
function power_rate_derivative(v,k,power,v₀)
    if (v-v₀)<0
        return 0.0
    else
        return k*power*(v-v₀)^(power-1)
    end
end
function deterministic_network_func(;
    #Weights
    w_r1er1e::Float64=1.25,w_r1er1i::Float64=0.65,w_r1ir1i::Float64=0.5,w_r1ir1e::Float64=1.2,
    w_r2er2e::Float64=1.25,w_r2er2i::Float64=0.65,w_r2ir2i::Float64=0.5,w_r2ir2e::Float64=1.2,
    w_r1er2e::Float64=0.4,w_r1er2i::Float64=0.0,w_r1ir2i::Float64=0.0,w_r1ir2e::Float64=0.5,
    w_r2er1e::Float64=1.1,w_r2er1i::Float64=0.0,w_r2ir1i::Float64=0.0,w_r2ir1e::Float64=1.1,
    #resting voltages
    vᵣ_r1e::Float64=-70.0,vᵣ_r1i::Float64=-70.0,vᵣ_r2e::Float64=-70.0,vᵣ_r2i::Float64=-70.0,
    #voltage thresholds
    vₜ_r1e::Float64=-70.0,vₜ_r1i::Float64=-70.0,vₜ_r2e::Float64=-70.0,vₜ_r2i::Float64=-70.0,
    #voltage time constants
    τ_r1e::Float64 = 20.0,τ_r1i::Float64 = 10.0,τ_r2e::Float64 = 20.0,τ_r2i::Float64 = 10.0,
    #rate variables
    k::Float64=0.3,power::Float64=2.0,
    #input amplitudes as a function of time. dim(h1_matrix) === dim(h2_matrix)
    h=1.0)
    #builds the weights
    weight_matrix = Array([w_r1er1e -1.0*w_r1er1i w_r1er2e -1.0*w_r1er2i;
                        w_r1ir1e -1.0*w_r1ir1i w_r1ir2e -1.0*w_r1ir2i;
                        w_r2er1e -1.0*w_r2er1i w_r2er2e -1.0*w_r2er2i;
                        w_r2ir1e -1.0*w_r2ir1i w_r2ir2e -1.0*w_r2ir2i])
    #builds the resting voltages
    resting_voltages = Array([vᵣ_r1e;
                            vᵣ_r1i;
                            vᵣ_r2e;
                            vᵣ_r2i])
    #builds the voltage thresholds
    voltage_thresholds = Array([vₜ_r1e;
                                vₜ_r1i;
                                vₜ_r2e;
                                vₜ_r2i])
    input_vec = [h;
                 h;
                 0;
                 0]
    voltage_time_constants = Array([τ_r1e;
                 τ_r1i;
                 τ_r2e;
                 τ_r2i])
    function voltage_eq!(voltage_vec)
        #col_voltage = reshape(voltage_vec,:,1)
        #display(-1*col_voltage + resting_voltages+input_vec)
        voltage_vec_2 = -1*voltage_vec + resting_voltages + weight_matrix*rectified_power_rate(voltage_vec,voltage_thresholds,k,power) + input_vec
        for i in 1:4
            voltage_vec[i] = voltage_vec_2[i]
        end
        return voltage_vec
    end
    parameters = Dict("weights" => weight_matrix,"resting_voltages"=>resting_voltages,
                      "voltage_thresholds"=>voltage_thresholds,"input_vec"=>input_vec,
                      "voltage_time_constants"=> voltage_time_constants,
                      "k"=>k,"power"=>power,"α" => w_r1er2e,"β" => w_r1ir2e)
    return voltage_eq!,parameters
end

function find_fixed_points(;α=0.5,β=0.5,h=1.0,iterations=10000,search_square_size=100)
    voltage_sols = []
    voltage_eq!,model_parameters = deterministic_network_func(w_r1er2e=α,w_r1ir2e=β,h=h)
    for _ in 1:iterations
        seed = round(search_square_size/2)*rand(Float64,(4)) + [-70.0,-70.0,-70.0,-70.0]
        soli = nlsolve(voltage_eq!,seed,autodiff=:forward,method=:newton,iterations=1000)
        if converged(soli)
            if length(voltage_sols) == 0
                push!(voltage_sols,soli.zero)
            else
                if all([norm(s - soli.zero) > 0.0001 for s in voltage_sols])
                    push!(voltage_sols,soli.zero)
                end
            end
        end
    end
    stability_eigs = [crit_point_stability(fixed_point,model_parameters) for fixed_point in voltage_sols]
    is_stable_vec = [all([real.(si)<0 for si in s]) for s in stability_eigs]
    rate_sols = [rectified_power_rate(v,model_parameters["voltage_thresholds"],model_parameters["k"],model_parameters["power"]) for v in voltage_sols]
    return voltage_sols,stability_eigs,is_stable_vec,rate_sols
end

function crit_point_stability(fixed_point,model_parameters)
    weights = model_parameters["weights"]
    voltage_thresholds = model_parameters["voltage_thresholds"]
    power = model_parameters["power"]
    k = model_parameters["k"]
    voltage_time_constants = model_parameters["voltage_time_constants"]
    derivative = [power_rate_derivative(fixed_point[i],k,power,voltage_thresholds[i]) for i in eachindex(fixed_point)]
    matrix = inv(diagm(voltage_time_constants))*(-1*I + weights*diagm(derivative))
    stability = eigvals(matrix)
    return stability
end
export find_fixed_points
end
#fixed_points,stability,is_stable = find_fixed_points(α=0.5,β=0.4,h=1.0)
#p = scatter(Dict(1.1=>[0.1,1.1],1.5 => [0.2,1.5,1.1]))
#savefig(p,"./test.png")
#initial_condition = [-70.0,-70.0,-70.0,-70.0]
#network_eq! = deterministic_network_func()
#display(network_eq!(initial_condition))
#display(initial_condition)
#nlsolve()