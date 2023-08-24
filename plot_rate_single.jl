using Plots
include("./multi_network_sim.jl")
using .SimulateOverTime
voltage_mean_frame,voltage_std_frame,rate_mean_frame,rate_std_frame,rate_mean_derivative_frame,rate_linear_fisher_frame,params = simulate_over_time(time_length_s=200.0,w_r2er1e=0.0,
w_r2ir1e=0.0,w_r1er2e=0.0,w_r1ir2e=0.0,input_vec=LinRange(0.0,20.0,21))

r1_v_std_plot = plot(voltage_std_frame[begin:end,:stim],voltage_std_frame[begin:end,:r1e],color=:red,marker=(:circle,4),xlabel="Input Strength",ylabel="Voltage Fixed Point Std",label="E",title="R1 Voltage Fixed Point Std Across Input Strengths")
plot!(r1_v_std_plot,voltage_std_frame[begin:end,:stim],voltage_std_frame[begin:end,:r1i],color=:blue,marker=(:circle,4),label="I")
savefig(r1_v_std_plot,"./r1_plots/voltage_std.png")

r1_v_mean_plot = plot(voltage_mean_frame[begin:end,:stim],voltage_mean_frame[begin:end,:r1e],color=:red,marker=(:circle,4),xlabel="Input Strength",ylabel="Voltage Fixed Point Mean",label="E",title="R1 Voltage Fixed Point Mean Across Input Strengths")
plot!(r1_v_mean_plot,voltage_mean_frame[begin:end,:stim],voltage_mean_frame[begin:end,:r1i],color=:blue,marker=(:circle,4),label="I")
savefig(r1_v_mean_plot,"./r1_plots/voltage_mean.png")

r1_r_mean_plot = plot(rate_mean_frame[begin:end,:stim],rate_mean_frame[begin:end,:r1e],color=:red,marker=(:circle,4),xlabel="Input Strength",ylabel="Rate Fixed Point Mean",label="E",title="R1 Rate Fixed Point Mean Across Input Strengths")
plot!(r1_r_mean_plot,rate_mean_frame[begin:end,:stim],rate_mean_frame[begin:end,:r1i],color=:blue,marker=(:circle,4),label="I")
savefig(r1_r_mean_plot,"./r1_plots/rate_mean.png")

r1_r_std_plot = plot(rate_std_frame[begin:end,:stim],rate_std_frame[begin:end,:r1e],color=:red,marker=(:circle,4),xlabel="Input Strength",ylabel="Rate Fixed Point Std",label="E",title="R1 Rate Fixed Point Std Across Input Strengths")
plot!(r1_r_std_plot,rate_std_frame[begin:end,:stim],rate_std_frame[begin:end,:r1i],color=:blue,marker=(:circle,4),label="I")
savefig(r1_r_std_plot,"./r1_plots/rate_std.png")

r1_rate_fisher_plot = plot(rate_linear_fisher_frame[begin:end,:stim],rate_linear_fisher_frame[begin:end,:r1e],color=:red,marker=(:circle,4),yaxis=:log,xlabel="Input Strength",ylabel="LFI",label="E",title="R1 Fixed Point LFI Across Input Strengths")
plot!(r1_rate_fisher_plot,rate_linear_fisher_frame[begin:end,:stim],rate_linear_fisher_frame[begin:end,:r1i],color=:blue,marker=(:circle,4),label="I")
savefig(r1_rate_fisher_plot,"./r1_plots/rate_linear_fisher_info.png")

