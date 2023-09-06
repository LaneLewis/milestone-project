include("./network_sim.jl")
using .SSNNetwork
using LinearAlgebra
using ProgressBars
using Plots
using Statistics
function find_fixed_point(;input_current=1.0,time_delta=0.1,divisions=100000,
                          w_r1er2i::Float64 = 0.0, w_r1ir2i::Float64 = 0.0, w_r1ir2e::Float64 = 1.1, w_r1er2e::Float64 = 1.1,
                          w_r2er1e::Float64 = 1.1, w_r2er1i::Float64 = 0.0,
                          w_r2ir1i::Float64 = 0.0, w_r2ir1e::Float64 = 1.1,
                          v₀_r1e::Float64 = -70.0, v₀_r1i::Float64 = -70.0, v₀_r2e::Float64 = -70.0, v₀_r2i::Float64 = -70.0,plt=false,show_stability=false)
    input_current = ones(1,divisions)*input_current
    zero_current = zeros(1,divisions)
    sim_voltages,sim_rates,param_dict = simulate_deterministic_network(;h1_vec=input_current,h2_vec=zero_current,
                                                    w_r1er2i = w_r1er2i, w_r1ir2i = w_r1ir2i, w_r1ir2e=w_r1ir2e, w_r1er2e=w_r1er2e,
                                                    w_r2er1e = w_r2er1e, w_r2er1i= w_r2er1i,
                                                    w_r2ir1i= w_r2ir1i, w_r2ir1e= w_r2ir1e,
                                                    v₀_r1e=v₀_r1e,v₀_r1i=v₀_r1i,v₀_r2e=v₀_r2e,v₀_r2i=v₀_r2i)
    if show_stability
        voltage_thresholds = param_dict["voltage_thresholds"]
        voltage_time_constants = param_dict["voltage_time_constants"]
        weights = param_dict["weights"]
        fixed_point_estimate = mean(sim_voltages,dims=2)
        derivative = [power_rate_derivative(fixed_point_estimate[i],0.3,2,voltage_thresholds[i,1]) for i in eachindex(fixed_point_estimate)]
        matrix = inv(diagm(voltage_time_constants))*(-1*I + weights*diagm(derivative))
        jacobian_stability = eigvals(matrix)
    end
    if plt
        p = plot(sim_voltages[1,begin:end])
        savefig(p,"./test_plot.png")
    end
    #p = plot(sim_voltages[1,begin:end])
    if !(norm(sim_voltages[begin:end,end-200] - sim_voltages[begin:end,end]) < 0.0001)
        display("did not converge")
    end
    return sim_voltages[begin:end,end]#,jacobian_stability
end
#1.263
function test_fixed(iterations,w_r1ir2e::Float64 = 1.1, w_r1er2e::Float64 = 1.1,w_r2er1e::Float64 = 1.1,w_r2ir1e::Float64 = 1.1,multiplicative_factor=10.0,norm_cutoff=0.0001,divisions=100000)
    v₀_r1e,v₀_r1i,v₀_r2e,v₀_r2i = -1*[70.0,70.0,70.0,70.0] + 2*multiplicative_factor*(rand(Float64,4) - 0.5*ones(4))
    fixed_point_0 = find_fixed_point(v₀_r1e=v₀_r1e,v₀_r1i=v₀_r1i,v₀_r2e=v₀_r2e,v₀_r2i=v₀_r2i,w_r1ir2e=w_r1ir2e, w_r1er2e=w_r1er2e,w_r2er1e=w_r2er1e,w_r2ir1e=w_r2ir1e,divisions=divisions)
    fixed_point_arrs = [fixed_point_0]
    #stability_arr = [stability]
    for i in ProgressBar(1:iterations)
        v₀_r1e,v₀_r1i,v₀_r2e,v₀_r2i = 2*multiplicative_factor*(rand(Float64,4) - 0.5*ones(4)) -1*[70.0,70.0,70.0,70.0]
        fixed_point = find_fixed_point(v₀_r1e=v₀_r1e,v₀_r1i=v₀_r1i,v₀_r2e=v₀_r2e,v₀_r2i=v₀_r2i,w_r1ir2e=w_r1ir2e, w_r1er2e=w_r1er2e,w_r2er1e=w_r2er1e,w_r2ir1e=w_r2ir1e,divisions=divisions)
        display(fixed_point)
        if !any([sum((fixed_point-point).^2) < norm_cutoff for point in fixed_point_arrs])
            push!(fixed_point_arrs,fixed_point)
        end
    end
    return fixed_point_arrs
end
function power_rate_derivative(v,k,power,v₀)
    if (v-v₀)<0
        return 0.0
    else
        return k*power*(v-v₀)^(power-1)
    end
end
#α = 1.1
#β = 1.1
#display(test_fixed(10,β,α))
#fixed_point_list = []
#v₀_r1e,v₀_r1i,v₀_r2e,v₀_r2i = [-68.1523, -68.0286, -68.5274, -68.4521]#2*100*(rand(Float64,4) - 0.5*ones(4)) + [-70.0,-70.0,-70.0,-70.0] #2*multiplicative_factor*(rand(Float64,4) - 0.5*ones(4)) - 
#fixed_point_0 = find_fixed_point(v₀_r1e=v₀_r1e,v₀_r1i=v₀_r1i,v₀_r2e=v₀_r2e,v₀_r2i=v₀_r2i,w_r1ir2e=β, w_r1er2e=α,w_r2er1e=1.1,w_r2ir1e=1.1,plt=true,show_stability=true)
#display(fixed_point_0)
#display(fixed_point_0)
#weights = [1.25 -0.65  1.0   -0.0;#
#           1.2   -0.5   0.5   -0.0;
#          1.1   -0.0   1.25  -0.65;
#           1.1   -0.0   1.2   -0.5]
#steady_state_voltages = [-70.0,-70.0,-70.0,-70.0]
#fixed_point_voltages = [-65.89272519387337,-64.87175087494137,-68.59835824914107,-65.16098370308524]
#derivative = [power_rate_derivative(v,0.3,2.0,-70.0) for v in fixed_point_voltages]
#display(derivative)
#matrix = (-1*I + weights*diagm(derivative))*inv(diagm([20.0,10.0,20.0,10.0]))
#display("eigs")
#display(eigvals(matrix))