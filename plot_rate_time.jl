include("./network_sim.jl")
using .SSNNetwork
using Plots
using Statistics
using ProgressBars
function simulate_experiment_run(;pre_stim_s=2.0,stim_presentation_s=6.0,divisions_per_ms=10,h=0.5,ei=1.0,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],
        saveloc="./graph.png",graph=false,stim_average_window_start=3.0,stim_average_window_end=3.1,show_window=false)
    stim_presentation_ms = round(Int,stim_presentation_s*1000)
    pre_stim_ms = round(Int,pre_stim_s*1000)
    average_window_start_index = round(Int,stim_average_window_start*1000+pre_stim_ms)*divisions_per_ms
    average_window_end_index = round(Int,stim_average_window_end*1000+pre_stim_ms)*divisions_per_ms
    stim_presentation_divisions = stim_presentation_ms*divisions_per_ms
    pre_stim_divisions = pre_stim_ms*divisions_per_ms
    h1s = hcat([zeros(Float64,(1,pre_stim_divisions)),h*ones(Float64,(1,stim_presentation_divisions))]...)
    h2s = zeros(Float64,(1,pre_stim_divisions + stim_presentation_divisions))
    sim_voltages,sim_rates,sim_inputs, params = simulate_network(h1_vec=h1s,h2_vec=h2s,w_r1er2e=β*ei,w_r2er1e=1.1,w_r2ir1e=1.1,w_r1ir2e=β,time_delta=1/divisions_per_ms,
                                                                v₀_r1e=initial_voltages[1],v₀_r1i=initial_voltages[2],v₀_r2e=initial_voltages[3],v₀_r2i=initial_voltages[4])
    volt_mean = mean(sim_voltages[begin:end,average_window_start_index:average_window_end_index],dims=2)
    rate_mean = mean(sim_rates[begin:end,average_window_start_index:average_window_end_index],dims=2)
    if graph
        rate_plot = plot(LinRange(0.0,stim_presentation_s+pre_stim_s,pre_stim_divisions + stim_presentation_divisions),sim_rates[1,begin:end],xaxis="Time (s)",
        yaxis="Rate",title="R1e Rate With Params: h=$h,β=$β,ei=$ei",color=:grey,label="Rate Trace")
        if show_window
            hline!(rate_plot,[rate_mean[1]],color=:red,linewidth=1,label="Window Average")
            vline!(rate_plot,[stim_average_window_start+pre_stim_s,stim_average_window_end+pre_stim_s],linestyle=:dash,color=:red,label=false,linewidth=1)
        end
        vline!(rate_plot,[pre_stim_s],linestyle=:dash,color=:blue,label=false,linewidth=5)
        savefig(rate_plot,saveloc)
    end
    return volt_mean,rate_mean
end

function find_mean_sd_for_h(;iterations=1000,pre_stim_s=1.0,stim_presentation_s=2.0,divisions_per_ms=10,h=0.5,ei=1.0,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],
    stim_average_window_start=1.5,stim_average_window_end=1.6)
    average_volt_arr = Array{Float64}(undef,iterations)
    average_rate_arr = Array{Float64}(undef,iterations)
    for i in 1:iterations
        average_volts,average_rates = simulate_experiment_run(pre_stim_s=pre_stim_s,stim_presentation_s=stim_presentation_s,divisions_per_ms=divisions_per_ms,h=h,ei=ei,β=β,initial_voltages=initial_voltages,
        stim_average_window_start=stim_average_window_start,stim_average_window_end=stim_average_window_end)
        average_volt_arr[i] = average_volts[1]
        average_rate_arr[i] = average_rates[1]
    end
    volt_mean_average = mean(average_volt_arr)
    rate_mean_average = mean(average_rate_arr)
    volt_mean_sd = std(average_volt_arr)
    rate_mean_sd = std(average_rate_arr)
    return volt_mean_average,volt_mean_sd,rate_mean_average,rate_mean_sd
end

function plot_average_mean_sd_lf(;iterations=5000,pre_stim_s=0.0,stim_presentation_s=1.5,divisions_per_ms=10,hs=LinRange(0.0,15.0,15),ei=0.5,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],
    stim_average_window_start=1.3,stim_average_window_end=1.4)
    volt_means = Array{Float64}(undef,length(hs))
    rate_means = Array{Float64}(undef,length(hs))
    volt_stds= Array{Float64}(undef,length(hs))
    rate_stds = Array{Float64}(undef,length(hs))
    for i in ProgressBars.ProgressBar(eachindex(hs))
        volt_mean,volt_std,rate_mean,rate_std = find_mean_sd_for_h(;iterations=iterations,pre_stim_s=pre_stim_s,stim_presentation_s=stim_presentation_s,divisions_per_ms=divisions_per_ms,h=hs[i],ei=ei,β=β,initial_voltages=initial_voltages,
        stim_average_window_start=stim_average_window_start,stim_average_window_end=stim_average_window_end)
        volt_means[i] = volt_mean
        rate_means[i] = rate_mean
        volt_stds[i] = volt_std
        rate_stds[i] = rate_std
    end
    voltage_mean_plot = plot(hs,volt_means,color=:red,title="Window Average Voltage Mean For h ∈[$(hs[1]),$(hs[end])],β=$β,E/I=$ei",xaxis="h",yaxis="Window Average Voltage Mean")
    voltage_std_plot = plot(hs,volt_stds,color=:red,title="Window Average Voltage Std For h ∈[$(hs[1]),$(hs[end])],β=$β,E/I=$ei",xaxis="h",yaxis="Window Average Voltage Std")
    rate_mean_plot = plot(hs,rate_means,color=:red,title="Window Average Rate Mean For h ∈[$(hs[1]),$(hs[end])],β=$β,E/I=$ei",xaxis="h",yaxis="Window Average Rate Mean")
    rate_std_plot = plot(hs,rate_stds,color=:red,title="Window Average Rate Mean For h ∈[$(hs[1]),$(hs[end])],β=$β,E/I=$ei",xaxis="h",yaxis="Window Average Rate Std")

    savefig(voltage_mean_plot,"./voltage_mean_plot.png")
    savefig(voltage_std_plot,"./voltage_std_plot.png")
    savefig(rate_mean_plot,"./rate_mean_plot.png")
    savefig(rate_std_plot,"./rate_std_plot.png")
end

function simulate_deterministic_experiment_run(;pre_stim_s=2.0,stim_presentation_s=6.0,divisions_per_ms=10,h=0.5,ei=1.0,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],
    saveloc="./graph.png",graph=false,show_window=false,stim_average_window_start=3.0,stim_average_window_end=3.1)
stim_presentation_ms = round(Int,stim_presentation_s*1000)
pre_stim_ms = round(Int,pre_stim_s*1000)
average_window_start_index = round(Int,stim_average_window_start*1000+pre_stim_ms)*divisions_per_ms
average_window_end_index = round(Int,stim_average_window_end*1000+pre_stim_ms)*divisions_per_ms
stim_presentation_divisions = stim_presentation_ms*divisions_per_ms
pre_stim_divisions = pre_stim_ms*divisions_per_ms
h1s = hcat([zeros(Float64,(1,pre_stim_divisions)),h*ones(Float64,(1,stim_presentation_divisions))]...)
h2s = zeros(Float64,(1,pre_stim_divisions + stim_presentation_divisions))
sim_voltages,sim_rates,params = simulate_deterministic_network(h1_vec=h1s,h2_vec=h2s,w_r1er2e=β*ei,w_r2er1e=1.1,w_r2ir1e=1.1,w_r1ir2e=β,time_delta=1/divisions_per_ms,
                                                            v₀_r1e=initial_voltages[1],v₀_r1i=initial_voltages[2],v₀_r2e=initial_voltages[3],v₀_r2i=initial_voltages[4])
volt_mean = mean(sim_voltages[begin:end,average_window_start_index:average_window_end_index],dims=2)
rate_mean = mean(sim_rates[begin:end,average_window_start_index:average_window_end_index],dims=2)
if graph
    rate_plot = plot(LinRange(0.0,stim_presentation_s+pre_stim_s,pre_stim_divisions + stim_presentation_divisions),sim_rates[1,begin:end],xaxis="Time (s)",
    yaxis="Rate",title="R1e Rate With Params: h=$h,β=$β,ei=$ei",color=:grey,label="Rate Trace")
    if show_window
        hline!(rate_plot,[rate_mean[1]],color=:red,linewidth=1,label="Window Average")
        vline!(rate_plot,[stim_average_window_start+pre_stim_s,stim_average_window_end+pre_stim_s],linestyle=:dash,color=:red,label=false,linewidth=1)
    end
    vline!(rate_plot,[pre_stim_s],linestyle=:dash,color=:blue,label=false,linewidth=5)
    savefig(rate_plot,saveloc)
end
return volt_mean,rate_mean
end
#loss of stability through ei
#simulate_deterministic_experiment_run(graph=true,h=5.0,ei=0.5,β=1.1,saveloc="./milestone_plots/transition/early_transition.png")
#simulate_deterministic_experiment_run(graph=true,h=5.0,ei=1.2,β=1.1,saveloc="./milestone_plots/transition/mid_transition.png")
#simulate_deterministic_experiment_run(graph=true,h=5.0,ei=1.5,β=1.1,saveloc="./milestone_plots/transition/late_transition.png")
#loss of stability through h
#simulate_deterministic_experiment_run(graph=true,h=0.5,ei=1.5,β=1.1,saveloc="./milestone_plots/unstable/early_transition.png")
#simulate_deterministic_experiment_run(graph=true,h=1.0,ei=1.5,β=1.1,saveloc="./milestone_plots/unstable/mid_transition.png")
#simulate_deterministic_experiment_run(graph=true,h=3.0,ei=1.5,β=1.1,saveloc="./milestone_plots/unstable/late_transition.png")
#make sure i know unbiased linear estimator stuff
#single_stability
#simulate_experiment_run(h=1.0,ei=0.5,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],saveloc="./milestone_plots/single_stable/trial_small_h.png",graph=true)
#simulate_experiment_run(h=5.0,ei=0.5,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],saveloc="./milestone_plots/single_stable/trial_large_h.png",graph=true)
#bistability 
#simulate_experiment_run(h=1.0,ei=1.0,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],saveloc="./milestone_plots/multi_stable/one_stable_point.png",graph=true)
#simulate_experiment_run(h=1.0,ei=1.0,β=1.1,initial_voltages=[-64.0;-60.0;-64.0;-60.0],saveloc="./milestone_plots/multi_stable/other_stable_point.png",graph=true)
#partial stability
simulate_experiment_run(h=8.0,ei=1.3,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],saveloc="./milestone_plots/partial_stable/low_h.png",graph=true)
simulate_experiment_run(h=15.0,ei=1.3,β=1.1,initial_voltages=[-70.0;-70.0;-70.0;-70.0],saveloc="./milestone_plots/partial_stable/high_h.png",graph=true)

#simulate_experiment_run(h=0.0,ei=0.85,β=1.1,initial_voltages=[-64.0;-60.0;-64.0;-60.0],saveloc="./graph2.png")