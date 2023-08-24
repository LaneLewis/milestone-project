module SimulateOverTime
    include("./network_sim.jl")
    using .SSNNetwork
    using Random,LinearAlgebra,HDF5
    using Statistics
    using DataFrames
    using ProgressBars
#this function finds the mean rate and standard deviation across an input given for a time length
#then ramped to the next input, etc. It outputs the means and standard deviations in the format
#across inputs across neurons.
function simulate_over_time(;time_length_s=10.0,divisions_per_ms=10,input_vec=LinRange(0.0,5.0,21),
    w_r2er1e=1.1,w_r2ir1e=1.1,w_r1er2e=0.0,w_r1ir2e=0.0)
    ms = round(Int,time_length_s*1000)
    single_input_divisions = ms*divisions_per_ms
    h1s = hcat([amplitude*ones(Float64,(1,single_input_divisions)) for amplitude in input_vec]...)
    h2s = zeros(Float64,(1,single_input_divisions*length(input_vec)))
    sim_voltages,sim_rates,sim_inputs, params = simulate_network(h1_vec=h1s,h2_vec=h2s,w_r1er2e=w_r1er2e,w_r2er1e=w_r2er1e,w_r2ir1e=w_r2ir1e,w_r1ir2e=w_r1ir2e)
    rate_means = Array{Float64,2}(undef,4,length(input_vec))
    rate_stds = Array{Float64,2}(undef,4,length(input_vec))
    voltage_means = Array{Float64,2}(undef,4,length(input_vec))
    voltage_stds = Array{Float64,2}(undef,4,length(input_vec))
    for i in 1:length(input_vec)
        lower_bound = (i-1)*single_input_divisions+1
        upper_bound = i*single_input_divisions
        #gets the values from the rates
        rates = sim_rates[begin:end,lower_bound:upper_bound]
        rate_means[begin:end,i] = mean(rates,dims=2)
        rate_stds[begin:end,i] = std(rates,dims=2)
        #gets the values from the sim_voltages
        voltages = sim_voltages[begin:end,lower_bound:upper_bound]
        voltage_means[begin:end,i] = mean(voltages,dims=2)
        voltage_stds[begin:end,i] = std(voltages,dims=2)
    end
    linear_fisher_info,rate_derivative_matrix = linear_fisher(rate_means,rate_stds,input_vec)
    rate_mean_derivative_frame = DataFrame("stim" => input_vec[begin:end-1], "r1e" => rate_derivative_matrix[begin:end,1],"r1i" => rate_derivative_matrix[begin:end,2],
                                "r2e" => rate_derivative_matrix[begin:end,3],"r2i" => rate_derivative_matrix[begin:end,4])
    rate_mean_frame = DataFrame("stim" => input_vec,"r1e" => rate_means[1,begin:end],"r1i" => rate_means[2,begin:end],
                                "r2e" => rate_means[3,begin:end],"r2i" => rate_means[4,begin:end])
    rate_std_frame = DataFrame("stim" => input_vec,"r1e" => rate_stds[1,begin:end],"r1i" => rate_stds[2,begin:end],
                               "r2e" => rate_stds[3,begin:end],"r2i" => rate_stds[4,begin:end])
    rate_linear_fisher_frame = DataFrame("stim" => input_vec[begin:end-1],"r1e" => linear_fisher_info[begin:end,1],"r1i" => linear_fisher_info[begin:end,2],
                                         "r2e" => linear_fisher_info[begin:end,3],"r2i" => linear_fisher_info[begin:end,4])
    voltage_mean_frame = DataFrame("stim" => input_vec,"r1e" => voltage_means[1,begin:end],"r1i" => voltage_means[2,begin:end],
                                   "r2e" => voltage_means[3,begin:end],"r2i" => voltage_means[4,begin:end])
    voltage_std_frame = DataFrame("stim" => input_vec,"r1e" => voltage_stds[1,begin:end],"r1i" => voltage_stds[2,begin:end],
                                "r2e" => voltage_stds[3,begin:end],"r2i" => voltage_stds[4,begin:end])
    return voltage_mean_frame,voltage_std_frame,rate_mean_frame,rate_std_frame,rate_mean_derivative_frame,rate_linear_fisher_frame,params
end

#This function builds the linear fisher information across inputs and neurons
function linear_fisher(means::Array{Float64},stds::Array{Float64},input_vec)
    means_diff = means[begin:end,2:end].-means[begin:end,begin:end-1]
    time_deltas = input_vec[2:end] - input_vec[begin:end-1]
    derivative_matrix = transpose(means_diff)./time_deltas
    linear_fisher = (derivative_matrix.^2)./(transpose(stds).^2)[begin:end-1,begin:end]
    return linear_fisher,derivative_matrix
end

#this function returns a dictionary of the ei ratio to corresponding data
function simulate_over_feedback_ei_ratio(;ei_ratios=LinRange(0.5,1.5,6),time_length_s=10.0,divisions_per_ms=10,input_vec=LinRange(0.0,5.0,21),
    w_r2er1e=1.1,w_r2ir1e=1.1,w_r1ir2e=0.5)
    sim_dict_arr = []
    for ei_ratio in ProgressBar(ei_ratios)
        w_r1er2e = ei_ratio*w_r1ir2e
        voltage_mean_frame,voltage_std_frame,rate_mean_frame,rate_std_frame,rate_derivative_frame,rate_linear_fisher_frame,params = simulate_over_time(;time_length_s=time_length_s,
                        divisions_per_ms=divisions_per_ms,input_vec=input_vec,w_r2er1e=w_r2er1e,w_r2ir1e=w_r2ir1e,w_r1ir2e=w_r1ir2e,w_r1er2e=w_r1er2e)
        sim_dict = Dict("voltage_mean_frame"=>voltage_mean_frame,"voltage_std_frame"=>voltage_std_frame,"rate_mean_frame"=>rate_mean_frame,"rate_std_frame"=>rate_std_frame,
                        "rate_linear_fisher_frame"=>rate_linear_fisher_frame,"rate_derivative_frame"=>rate_derivative_frame,"params"=>params,"ei_ratio"=>ei_ratios)
        push!(sim_dict_arr,sim_dict)
        Base.GC.gc()
    end
    return sim_dict_arr
end
export simulate_over_time,simulate_over_feedback_ei_ratio
end
#simulate_over_time()