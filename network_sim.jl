module SSNNetwork
using Random,LinearAlgebra,HDF5
using DataFrames
using CSV
using StatsPlots
using Statistics
using Pkg
using ProgressBars
import Base.Threads
rng = MersenneTwister(780)
    function zero_or_one_row_matrix(rows::Int64,columns::Int64,non_zero_row_indicies::Array{Int64})::Array{Float64,2}
        zero_inputs = zeros(Float64,(rows,columns))
        for i ∈ non_zero_row_indicies
            zero_inputs[i,begin:end] = ones(Float64,columns)
        end
        return zero_inputs
    end

    function matrix_root(matrix::Array{Float64})::Array{Float64,2}
        S = eigvecs(matrix)
        root_eigs = eigvals(matrix).^(0.5)
        return S*diagm(root_eigs)*transpose(S)
    end

    function rectified_power_rate(voltage_vec::Array{Float64},voltage_threshold_vec::Array{Float64},k::Float64,power::Float64)::Array{Float64}
        rectified_voltage = map(x -> max(x,0),voltage_vec.-voltage_threshold_vec)
        return k*rectified_voltage.^power
    end

    function ornstein_uhlenbeck_process(noise_cov::Array{Float64,2},time_const::Float64,initial_cond_vec::Array{Float64,1},
                                        time_start::Float64,time_end::Float64,divisions::Int64)::Array{Float64,2}
        time_delta = (time_end-time_start)/divisions
        noise_dimension = size(noise_cov,1)
        out_arr = Array{Float64, 2}(undef,noise_dimension,divisions)
        noise_mult = matrix_root(noise_cov)*sqrt(2*time_const)
        for i ∈ 1:divisions
            normal_sample = randn(rng,Float64,(noise_dimension,1))
            noise = noise_mult*normal_sample
            if i == 1
                out_arr[begin:end,i] = initial_cond_vec  + -1*initial_cond_vec*(time_delta/time_const) + noise*(sqrt(time_delta)/time_const)
            else
                out_arr[begin:end,i] = out_arr[begin:end,i-1] + -1*out_arr[begin:end,i-1]*(time_delta/time_const) + noise*(sqrt(time_delta)/time_const)
            end
        end
        return out_arr
    end

    function ssn_diff_eq(weight_matrix::Array{Float64,2},time_const_vec::Array{Float64,1},previous_voltage_vec::Array{Float64,1},
                        resting_voltage_vec::Array{Float64,1},previous_rate_vec::Array{Float64,1},external_input::Array{Float64,1})::Array{Float64,1}
        rate_input = weight_matrix*previous_rate_vec
        dv_dt = (-1*previous_voltage_vec + resting_voltage_vec + external_input + rate_input)./time_const_vec
        return dv_dt
    end

    function euler_simulate_ssn(weight_matrix::Array{Float64,2},voltage_time_const_vec::Array{Float64,1},initial_voltage_vec::Array{Float64,1},
                            resting_voltage_vec::Array{Float64,1},input_matrix::Array{Float64,2},time_delta::Float64,rate_function::Function)
        #input matrix has size (network_size,divisions)
        divisions = size(input_matrix,2)
        network_size = length(initial_voltage_vec)
        voltage_out = Array{Float64,2}(undef,network_size,divisions)
        rate_out = Array{Float64,2}(undef,network_size,divisions)
        initial_rate_vec = rate_function(initial_voltage_vec)
        voltage_out[begin:end,1] = initial_voltage_vec + ssn_diff_eq(weight_matrix,voltage_time_const_vec,initial_voltage_vec,
                                            resting_voltage_vec,initial_rate_vec,input_matrix[begin:end,1])*time_delta
        rate_out[begin:end,1] = rate_function(voltage_out[begin:end,1])
        for i ∈ 2:size(input_matrix,2)
            voltage_out[begin:end,i] = voltage_out[begin:end,i-1] + ssn_diff_eq(weight_matrix,voltage_time_const_vec,voltage_out[begin:end,i-1],
                                                resting_voltage_vec,rate_out[begin:end,i-1],input_matrix[begin:end,i])*time_delta
            rate_out[begin:end,i] = rate_function(voltage_out[begin:end,i])
        end
        return voltage_out,rate_out
    end

    function simulate_network(;
        #Weights
        w_r1er1e::Float64=1.25,w_r1er1i::Float64=0.65,w_r1ir1i::Float64=0.5,w_r1ir1e::Float64=1.2,
        w_r2er2e::Float64=1.25,w_r2er2i::Float64=0.65,w_r2ir2i::Float64=0.5,w_r2ir2e::Float64=1.2,
        w_r1er2e::Float64=0.0,w_r1er2i::Float64=0.0,w_r1ir2i::Float64=0.0,w_r1ir2e::Float64=0.0,
        w_r2er1e::Float64=0.0,w_r2er1i::Float64=0.0,w_r2ir1i::Float64=0.0,w_r2ir1e::Float64=0.0,
        #Initial voltages
        v₀_r1e::Float64=-70.0,v₀_r1i::Float64=-70.0,v₀_r2e::Float64=-70.0,v₀_r2i::Float64=-70.0,
        #resting voltages
        vᵣ_r1e::Float64=-70.0,vᵣ_r1i::Float64=-70.0,vᵣ_r2e::Float64=-70.0,vᵣ_r2i::Float64=-70.0,
        #voltage time constants
        τ_r1e::Float64 = 20.0,τ_r1i::Float64 = 10.0,τ_r2e::Float64 = 20.0,τ_r2i::Float64 = 10.0,
        #voltage thresholds
        vₜ_r1e::Float64=-70.0,vₜ_r1i::Float64=-70.0,vₜ_r2e::Float64=-70.0,vₜ_r2i::Float64=-70.0,
        #rate variables
        k::Float64=0.3,power::Float64=2.0,
        #noise initial conditions
        η_r1e::Float64 = 0.0, η_r1i::Float64 = 0.0,η_r2e::Float64 = 0.0,η_r2i::Float64 = 0.0,
        #noise time constant
        τₙ::Float64 = 50.0,
        #noise stds
        noise_std_r1e::Float64 = 0.2,noise_std_r1i::Float64 = 0.1,noise_std_r2e::Float64 = 0.2,noise_std_r2i::Float64 = 0.1,
        #input amplitudes as a function of time. dim(h1_matrix) === dim(h2_matrix)
        h1_vec::Array{Float64}, h2_vec::Array{Float64},
        #euler integration variables
        time_delta::Float64=0.1)
        #builds the time components
        divisions = length(h1_vec)
        time_end::Float64 = divisions*time_delta
        time_start = 0.0
        #builts the stochastic part of the inputs
        noise_initial_conditions = Array([η_r1e;
                                        η_r1i;
                                        η_r2e;
                                        η_r2i;])
        #*np.sqrt(1+r1ETau/noiseTau) add this in!!!!
        noise_cov = diagm([noise_std_r1e*sqrt(1+τ_r1e/τₙ),
                           noise_std_r1i*sqrt(1+τ_r1i/τₙ),
                           noise_std_r2e*sqrt(1+τ_r2e/τₙ),
                           noise_std_r2i*sqrt(1+τ_r2i/τₙ)]).^2
            
        variable_input_currents = ornstein_uhlenbeck_process(noise_cov,τₙ,noise_initial_conditions,time_start,time_end,divisions)
        #builds the input matrix
        fixed_input_currents = [h1_vec;
                                h1_vec;
                                h2_vec;
                                h2_vec]
        total_input_currents = fixed_input_currents + variable_input_currents
        #builds the weights
        weight_matrix = Array([w_r1er1e -1.0*w_r1er1i w_r1er2e -1.0*w_r1er2i;
                            w_r1ir1e -1.0*w_r1ir1i w_r1ir2e -1.0*w_r1ir2i;
                            w_r2er1e -1.0*w_r2er1i w_r2er2e -1.0*w_r2er2i;
                            w_r2ir1e -1.0*w_r2ir1i w_r2ir2e -1.0*w_r2ir2i])
        #builds the initial voltages
        initial_voltages = Array([v₀_r1e;
                                v₀_r1i;
                                v₀_r2e;
                                v₀_r2i])
        #builds the resting voltages
        resting_voltages = Array([vᵣ_r1e;
                                vᵣ_r1i;
                                vᵣ_r2e;
                                vᵣ_r2i])
        #builds the voltage time constants
        voltage_time_constants = Array([τ_r1e;
                                        τ_r1i;
                                        τ_r2e;
                                        τ_r2i])
        #builds the voltage thresholds
        voltage_thresholds = Array([vₜ_r1e;
                                    vₜ_r1i;
                                    vₜ_r2e;
                                    vₜ_r2i])
        #makes the voltage to rate function 
        sim_voltages,sim_rates = euler_simulate_ssn(weight_matrix,voltage_time_constants,initial_voltages,
                                                    resting_voltages,total_input_currents,time_delta,
                                                    x->rectified_power_rate(x,voltage_thresholds,k,power))
        param_dict=Dict("noise_initial_conditions"=>noise_initial_conditions,"weights"=>weight_matrix,"initial_voltages"=>initial_voltages,
                        "resting_voltages"=>resting_voltages,"voltage_time_constants"=>voltage_time_constants,"voltage_thresholds"=>voltage_thresholds,
                        "time_delta"=>time_delta,"noise_cov"=>noise_cov)
        return sim_voltages,sim_rates,total_input_currents,param_dict
    end

    function simulate_deterministic_network(;
        #Weights
        w_r1er1e::Float64=1.25,w_r1er1i::Float64=0.65,w_r1ir1i::Float64=0.5,w_r1ir1e::Float64=1.2,
        w_r2er2e::Float64=1.25,w_r2er2i::Float64=0.65,w_r2ir2i::Float64=0.5,w_r2ir2e::Float64=1.2,
        w_r1er2e::Float64=0.0,w_r1er2i::Float64=0.0,w_r1ir2i::Float64=0.0,w_r1ir2e::Float64=0.0,
        w_r2er1e::Float64=0.0,w_r2er1i::Float64=0.0,w_r2ir1i::Float64=0.0,w_r2ir1e::Float64=0.0,
        #Initial voltages
        v₀_r1e::Float64=-70.0,v₀_r1i::Float64=-70.0,v₀_r2e::Float64=-70.0,v₀_r2i::Float64=-70.0,
        #resting voltages
        vᵣ_r1e::Float64=-70.0,vᵣ_r1i::Float64=-70.0,vᵣ_r2e::Float64=-70.0,vᵣ_r2i::Float64=-70.0,
        #voltage time constants
        τ_r1e::Float64 = 20.0,τ_r1i::Float64 = 10.0,τ_r2e::Float64 = 20.0,τ_r2i::Float64 = 10.0,
        #voltage thresholds
        vₜ_r1e::Float64=-70.0,vₜ_r1i::Float64=-70.0,vₜ_r2e::Float64=-70.0,vₜ_r2i::Float64=-70.0,
        #rate variables
        k::Float64=0.3,power::Float64=2.0,
        #input amplitudes as a function of time. dim(h1_matrix) === dim(h2_matrix)
        h1_vec::Array{Float64}, h2_vec::Array{Float64},
        #euler integration variables
        time_delta::Float64=0.1)
        #builds the time components
        divisions = length(h1_vec)
        time_end::Float64 = divisions*time_delta
        time_start = 0.0

        #builds the input matrix
        total_input_currents = [h1_vec;
                                h1_vec;
                                h2_vec;
                                h2_vec]
        #builds the weights
        weight_matrix = Array([w_r1er1e -1.0*w_r1er1i w_r1er2e -1.0*w_r1er2i;
                            w_r1ir1e -1.0*w_r1ir1i w_r1ir2e -1.0*w_r1ir2i;
                            w_r2er1e -1.0*w_r2er1i w_r2er2e -1.0*w_r2er2i;
                            w_r2ir1e -1.0*w_r2ir1i w_r2ir2e -1.0*w_r2ir2i])
        #builds the initial voltages
        initial_voltages = Array([v₀_r1e;
                                v₀_r1i;
                                v₀_r2e;
                                v₀_r2i])
        #builds the resting voltages
        resting_voltages = Array([vᵣ_r1e;
                                vᵣ_r1i;
                                vᵣ_r2e;
                                vᵣ_r2i])
        #builds the voltage time constants
        voltage_time_constants = Array([τ_r1e;
                                        τ_r1i;
                                        τ_r2e;
                                        τ_r2i])
        #builds the voltage thresholds
        voltage_thresholds = Array([vₜ_r1e;
                                    vₜ_r1i;
                                    vₜ_r2e;
                                    vₜ_r2i])
        #makes the voltage to rate function 
        sim_voltages,sim_rates = euler_simulate_ssn(weight_matrix,voltage_time_constants,initial_voltages,
                                                    resting_voltages,total_input_currents,time_delta,
                                                    x->rectified_power_rate(x,voltage_thresholds,k,power))
        param_dict=Dict("weights"=>weight_matrix,"initial_voltages"=>initial_voltages,
                        "resting_voltages"=>resting_voltages,"voltage_time_constants"=>voltage_time_constants,"voltage_thresholds"=>voltage_thresholds,
                        "time_delta"=>time_delta)
        return sim_voltages,sim_rates,param_dict
    end
    export simulate_network,simulate_deterministic_network
end
