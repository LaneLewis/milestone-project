using Plots,CSV
include("./multi_network_sim.jl")
using .SimulateOverTime


function region_color_gradient(;r1e_color = :red,r1i_color=:blue,r2e_color=:orange,r2i_color=:purple,divisions=6)
    if divisions == 1
        r1e_colors = [r1e_color]
        r1i_colors = [r1i_color]
        r2e_colors = [r2e_color]
        r2i_colors = [r2i_color]
    else
        r1e_colors = palette([:black,r1e_color],divisions)
        r1i_colors = palette([:black,r1i_color],divisions)
        r2e_colors = palette([:black,r2e_color],divisions)
        r2i_colors = palette([:black,r2i_color],divisions)
    end
    return r1e_colors,r1i_colors,r2e_colors,r2i_colors
end

function plot_ei_rate_mean(ei_ratio_dict_arr;linewidth=2,combined_save_loc="./plots/ei_rate_mean.png",r1e_save_loc="./plots/r1e_ei_rate_mean.png")
    r1e_colors,r1i_colors,r2e_colors,r2i_colors = region_color_gradient(;divisions=length(ei_ratio_dict_arr))
    #generates the initial graphs which are then modified
    initial_mean_rate = ei_ratio_dict_arr[1]["rate_mean_frame"]
    r1e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1e],color=r1e_colors[1],legend=false,linewidth=linewidth,title="R1 Excitatory",ylabel="Fixed Point Rate")
    r1i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1i],color=r1i_colors[1],legend=false,linewidth=linewidth,title="R1 Inhibitory",ylabel="Fixed Point Rate",xlabel="Input Strength")
    r2e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2e],color=r2e_colors[1],legend=false,linewidth=linewidth,title="R2 Excitatory")
    r2i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2i],color=r2i_colors[1],legend=false,linewidth=linewidth,title="R2 Inhibitory",xlabel="Input Strength")
    for i in eachindex(ei_ratio_dict_arr)[begin+1:end]
        mean_rate = ei_ratio_dict_arr[i]["rate_mean_frame"]
        plot!(r1e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1e],color=r1e_colors[i],legend=false,linewidth=linewidth)
        plot!(r1i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1i],color=r1i_colors[i],legend=false,linewidth=linewidth)
        plot!(r2e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2e],color=r2e_colors[i],legend=false,linewidth=linewidth)
        plot!(r2i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2i],color=r2i_colors[i],legend=false,linewidth=linewidth)
    end
    subplots = plot(r1e_plot,r2e_plot,r1i_plot,r2i_plot,layout=(2,2))
    combined_plot = plot(subplots,plot_title="Network Fixed Point Rate With Modification Of Feedback EI Ratio Across Regions\n",thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plot,combined_save_loc)
    title!(r1e_plot,"Region 1 Excitatory Network Fixed Point Rate\n With Modification Of Feedback EI Ratio")
    xlabel!(r1e_plot,"Input Strength")
    savefig(r1e_plot,r1e_save_loc)
end

function plot_ei_rate_std(ei_ratio_dict_arr;linewidth=2,combined_save_loc="./plots/ei_rate_std.png",r1e_save_loc="./plots/r1e_ei_rate_std.png")
    r1e_colors,r1i_colors,r2e_colors,r2i_colors = region_color_gradient(;divisions=length(ei_ratio_dict_arr))
    #generates the initial graphs which are then modified
    initial_mean_rate = ei_ratio_dict_arr[1]["rate_std_frame"]
    r1e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1e],color=r1e_colors[1],legend=false,linewidth=linewidth,title="R1 Excitatory",ylabel="Fixed Point Std")
    r1i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1i],color=r1i_colors[1],legend=false,linewidth=linewidth,title="R1 Inhibitory",ylabel="Fixed Point Std",xlabel="Input Strength")
    r2e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2e],color=r2e_colors[1],legend=false,linewidth=linewidth,title="R2 Excitatory")
    r2i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2i],color=r2i_colors[1],legend=false,linewidth=linewidth,title="R2 Inhibitory",xlabel="Input Strength")
    for i in eachindex(ei_ratio_dict_arr)[begin+1:end]
        mean_rate = ei_ratio_dict_arr[i]["rate_std_frame"]
        plot!(r1e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1e],color=r1e_colors[i],legend=false,linewidth=linewidth)
        plot!(r1i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1i],color=r1i_colors[i],legend=false,linewidth=linewidth)
        plot!(r2e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2e],color=r2e_colors[i],legend=false,linewidth=linewidth)
        plot!(r2i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2i],color=r2i_colors[i],legend=false,linewidth=linewidth)
    end
    subplots = plot(r1e_plot,r2e_plot,r1i_plot,r2i_plot,layout=(2,2))
    combined_plot = plot(subplots,plot_title="Network Fixed Point Rate Std With Modification Of Feedback EI Ratio Across Regions\n",thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plot,combined_save_loc)
    title!(r1e_plot,"Region 1 Excitatory Fixed Point Rate Std \n With Modification Of Feedback EI Ratio")
    xlabel!(r1e_plot,"Input Strength")
    savefig(r1e_plot,r1e_save_loc)
end

function plot_ei_rate_fisher(ei_ratio_dict_arr;linewidth=2,combined_save_loc="./plots/ei_rate_fisher.png",r1e_save_loc="./plots/r1e_ei_rate_fisher.png")
    r1e_colors,r1i_colors,r2e_colors,r2i_colors = region_color_gradient(;divisions=length(ei_ratio_dict_arr))
    #generates the initial graphs which are then modified
    initial_mean_rate = ei_ratio_dict_arr[1]["rate_linear_fisher_frame"]
    r1e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1e],color=r1e_colors[1],legend=false,linewidth=linewidth,title="R1 Excitatory",ylabel="Fixed Point LFI",yaxis=:log10)
    r1i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1i],color=r1i_colors[1],legend=false,linewidth=linewidth,title="R1 Inhibitory",ylabel="Fixed Point LFI",yaxis=:log10,xlabel="Input Strength")
    r2e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2e],color=r2e_colors[1],legend=false,linewidth=linewidth,title="R2 Excitatory",yaxis=:log10)
    r2i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2i],color=r2i_colors[1],legend=false,linewidth=linewidth,title="R2 Inhibitory",xlabel="Input Strength",yaxis=:log10)
    for i in eachindex(ei_ratio_dict_arr)[begin+1:end]
        mean_rate = ei_ratio_dict_arr[i]["rate_linear_fisher_frame"]
        plot!(r1e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1e],color=r1e_colors[i],legend=false,linewidth=linewidth)
        plot!(r1i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1i],color=r1i_colors[i],legend=false,linewidth=linewidth)
        plot!(r2e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2e],color=r2e_colors[i],legend=false,linewidth=linewidth)
        plot!(r2i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2i],color=r2i_colors[i],legend=false,linewidth=linewidth)
    end
    subplots = plot(r1e_plot,r2e_plot,r1i_plot,r2i_plot,layout=(2,2))
    combined_plot = plot(subplots,plot_title="Network Fixed Point LFI With Modification Of Feedback EI Ratio Across Regions\n",thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plot,combined_save_loc)
    title!(r1e_plot,"Region 1 Excitatory Fixed Point LFI \n With Modification Of Feedback EI Ratio")
    xlabel!(r1e_plot,"Input Strength")
    savefig(r1e_plot,r1e_save_loc)
end

function plot_ei_rate_derivative(ei_ratio_dict_arr;linewidth=2,combined_save_loc="./plots/ei_rate_derivative_fisher.png",r1e_save_loc="./plots/r1e_ei_rate_derivative_fisher.png")
    r1e_colors,r1i_colors,r2e_colors,r2i_colors = region_color_gradient(;divisions=length(ei_ratio_dict_arr))
    #generates the initial graphs which are then modified
    initial_mean_rate = ei_ratio_dict_arr[1]["rate_derivative_frame"]
    r1e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1e],color=r1e_colors[1],legend=false,linewidth=linewidth,title="R1 Excitatory",ylabel="Fixed Point Derivative")
    r1i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r1i],color=r1i_colors[1],legend=false,linewidth=linewidth,title="R1 Inhibitory",ylabel="Fixed Point Derivative",xlabel="Input Strength")
    r2e_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2e],color=r2e_colors[1],legend=false,linewidth=linewidth,title="R2 Excitatory",)
    r2i_plot = plot(initial_mean_rate[begin:end,:stim],initial_mean_rate[begin:end,:r2i],color=r2i_colors[1],legend=false,linewidth=linewidth,title="R2 Inhibitory",xlabel="Input Strength")
    for i in eachindex(ei_ratio_dict_arr)[begin+1:end]
        mean_rate = ei_ratio_dict_arr[i]["rate_derivative_frame"]
        plot!(r1e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1e],color=r1e_colors[i],legend=false,linewidth=linewidth)
        plot!(r1i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r1i],color=r1i_colors[i],legend=false,linewidth=linewidth)
        plot!(r2e_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2e],color=r2e_colors[i],legend=false,linewidth=linewidth)
        plot!(r2i_plot,initial_mean_rate[begin:end,:stim],mean_rate[begin:end,:r2i],color=r2i_colors[i],legend=false,linewidth=linewidth)
    end
    subplots = plot(r1e_plot,r2e_plot,r1i_plot,r2i_plot,layout=(2,2))
    combined_plot = plot(subplots,plot_title="Network Fixed Point Derivative With Modification Of Feedback EI Ratio Across Regions\n",thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plot,combined_save_loc)
    title!(r1e_plot,"Region 1 Excitatory Fixed Point Derivative \n With Modification Of Feedback EI Ratio")
    xlabel!(r1e_plot,"Input Strength")
    savefig(r1e_plot,r1e_save_loc)
end

out_dict = simulate_over_feedback_ei_ratio(time_length_s=200.0; w_r2er1e=1.1, w_r2ir1e=1.1,
                            w_r1ir2e=1.1,input_vec=LinRange(0.0,10,21),ei_ratios = LinRange(0.1, 1.1, 10))
plot_ei_rate_mean(out_dict)
plot_ei_rate_std(out_dict)
plot_ei_rate_fisher(out_dict)
plot_ei_rate_derivative(out_dict)
#plot_rate_ei_mean(out_dict)
#plot_rate_ei_std(out_dict)
#plot_linear_fisher_ei(out_dict)