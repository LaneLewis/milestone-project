include("./multi_network_sim.jl")
include("./network_sim.jl")
using .SimulateOverTime
using LinearAlgebra
using DataFrames
using Plots

function power_rate_derivative(v,k,power,v₀)
    if (v-v₀)<0
        return 0.0
    else
        return k*power*(v-v₀)^(power-1)
    end
end

function fixed_point_derivative_frame(voltage_mean_frame,voltage_thresholds;k=0.3,power=2)
    derivative_frame = copy(voltage_mean_frame)
    columns = [:r1e,:r1i,:r2e,:r2i]
    for col in eachindex(columns)
        transform!(derivative_frame,columns[col]=>ByRow(x->power_rate_derivative(x,k,power,voltage_thresholds[col]))=>columns[col])
    end
    return derivative_frame
end

function fixed_point_linearized_response(weight_matrix,derivative_frame,rate_mean_frame)
    inputs_delta = derivative_frame[begin+1:end,:stim] - derivative_frame[begin:end-1,:stim]
    derivative_matrix = Matrix(derivative_frame[begin:end,[:r1e,:r1i,:r2e,:r2i]])
    linearized_rate_matrix = Matrix{Float64}(undef,length(inputs_delta),4)
    linearized_rate_derivative_matrix = Matrix{Float64}(undef,length(inputs_delta),4)
    for i in eachindex(inputs_delta)
        response_double_vec = inputs_delta[i]*diagm(derivative_matrix[i,begin:end])*inv(I - weight_matrix*diagm(derivative_matrix[i,begin:end]))*[[1],[1],[0],[0]]
        linearized_rate_matrix[i,begin:end] = [(response_double_vec...)...] + [values(rate_mean_frame[i,[:r1e,:r1i,:r2e,:r2i]])...]
        linearized_rate_derivative_matrix[i,begin:end] = [(response_double_vec...)...]
    end
    linearized_rate_frame = DataFrame(:stim=>derivative_frame[begin+1:end,:stim],:r1e=>linearized_rate_matrix[begin:end,1],:r1i=>linearized_rate_matrix[begin:end,2],
                                    :r2e=>linearized_rate_matrix[begin:end,3],:r2i=>linearized_rate_matrix[begin:end,4])
    linearized_rate_derivative_frame = DataFrame(:stim=>derivative_frame[begin+1:end,:stim],:r1e=>linearized_rate_derivative_matrix[begin:end,1],
                                    :r1i=>linearized_rate_derivative_matrix[begin:end,2],:r2e=>linearized_rate_derivative_matrix[begin:end,3],
                                    :r2i=>linearized_rate_derivative_matrix[begin:end,4])
    return linearized_rate_frame, linearized_rate_derivative_frame
end

function fixed_point_linearized_response_frame(voltage_mean_frame,rate_mean_frame, params)
    linearized_rate_frame, linearized_rate_derivative_frame = fixed_point_linearized_response(params["weights"],fixed_point_derivative_frame(voltage_mean_frame,params["voltage_thresholds"]),rate_mean_frame)
    return linearized_rate_frame, linearized_rate_derivative_frame
end
function simulate_and_plot_linearized(;saveloc="./r1_plots/linearized_approx.png",time_length_s=150.0,input_vec=LinRange(0.0, 16.0, 33),w_r1er2e = 0.0, w_r1ir2e = 0.0)
    voltage_mean_frame,voltage_std_frame,rate_mean_frame,rate_std_frame,rate_mean_derivative_frame,rate_linear_fisher_frame,params = simulate_over_time(time_length_s=time_length_s,input_vec = input_vec, w_r1er2e = w_r1er2e, w_r1ir2e = w_r1ir2e)
    linearized_frame,linearized_rate_derivative_frame = fixed_point_linearized_response_frame(voltage_mean_frame, rate_mean_frame,params)
    rate_diff_matrix = Matrix(rate_mean_frame[begin+1:end,[:r1e,:r1i,:r2e,:r2i]]) - Matrix(rate_mean_frame[begin:end-1,[:r1e,:r1i,:r2e,:r2i]])
    rate_diff_frame = DataFrame(:stim => rate_mean_frame[begin+1:end,:stim], :r1e => rate_diff_matrix[begin:end,1],
                                :r1i => rate_diff_matrix[begin:end,2],:r2e => rate_diff_matrix[begin:end,3],:r2i => rate_diff_matrix[begin:end,4])
    p = plot(linearized_rate_derivative_frame[!,:stim],linearized_rate_derivative_frame[!,:r1e],title="Linearized Fixed Point Derivative Vs Simulated",color=:green,label="Linearized",xlabel="Input Strength")
    plot!(p,rate_diff_frame[!,:stim],rate_diff_frame[!,:r1e],color=:black,label="Simulated",ylabel="Fixed Point Derivative")
    savefig(p,saveloc)
end

simulate_and_plot_linearized(;saveloc="./r1_plots/linearized_approx.png",time_length_s=150.0,input_vec=LinRange(0.0, 20.0, 21),w_r1er2e = 0.0, w_r1ir2e = 0.0)
Base.GC.gc()
simulate_and_plot_linearized(;saveloc="./plots/feedback_linearized_approx_1.1_0.5.png",time_length_s=150.0,input_vec=LinRange(0.0, 20.0, 21),w_r1er2e = 1.1, w_r1ir2e = 0.5)
Base.GC.gc()
simulate_and_plot_linearized(;saveloc="./plots/feedback_linearized_approx_1.1_1.1.png",time_length_s=150.0,input_vec=LinRange(0.0, 20.0, 21),w_r1er2e = 1.1, w_r1ir2e = 1.1)
Base.GC.gc()
simulate_and_plot_linearized(;saveloc="./plots/feedback_linearized_approx_1.1_1.5.png",time_length_s=150.0,input_vec=LinRange(0.0, 20.0, 21),w_r1er2e = 1.1, w_r1ir2e = 1.5)
Base.GC.gc()

#p = plot(linearized_frame[!,:stim],linearized_frame[!,:r1e],color=:red)
#plot!(p,rate_mean_frame[begin+1:end,:stim],rate_mean_frame[begin+1:end,:r1e],color=:blue)
#p1 = plot(linearized_rate_derivative_frame[!,:stim],linearized_rate_derivative_frame[!,:r1e],color=:red)
#plot!(p1,rate_mean_derivative_frame[!,:stim],rate_mean_derivative_frame[!,:r1e],color=:blue)
#savefig(p,"./plots/linearized.png")
#savefig(p1,"./plots/derivative.png")
#voltage_mean_frame,rate_mean_frame,rate_std_frame,rate_linear_fisher_frame,params = simulate_over_time()
#voltage_thresholds = params["voltage_thresholds"]
#res = fixed_point_linearized_response(params["weights"],fixed_point_derivative_frame(voltage_mean_frame,voltage_thresholds))
#display(res)
