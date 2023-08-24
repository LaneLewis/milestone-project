include("./critical_points.jl")
using .Stability
using Plots
using ProgressBars
#plots an array of dicts mapping floats to floats. All these recieve a single label and the same arguments
function plot_arr_dict!(initial_plot,arr_dict;label="",kwargs...)
    seed_plot = initial_plot
    if length(arr_dict) >= 1
        for i in eachindex(arr_dict)
            if i == 1
                plot_label = label
            else
                plot_label = false
            end
            plot!(seed_plot,arr_dict[i];label=plot_label,kwargs...)
        end
    end
    return seed_plot
end

function plot_quadrant(stable_dicts,unstable_dicts;quadrant=1,color=:red,xlabel="α",ylabel="Critical Point Voltage")
    titles = ["R1e","R1i","R2e","R2i"]
    quadrant_stable_dict = [index_dictionary(dictionary,quadrant) for dictionary in stable_dicts]
    quadrant_unstable_dict = [index_dictionary(dictionary,quadrant) for dictionary in unstable_dicts]
    stability_plot = plot_arr_dict!(plot(),quadrant_stable_dict,seriestype =:scatter,color=color,label="Stable",title=titles[quadrant])
    plot_arr_dict!(stability_plot,quadrant_unstable_dict,seriestype =:scatter,color=:grey,label="Unstable")
    xaxis!(stability_plot,xlabel)
    yaxis!(stability_plot,ylabel)
end

function index_dictionary(dictionary,index)
    index_dict = Dict()
    for (key,value) in dictionary
        index_dict[key] = value[index]
    end
    return index_dict
end

function plot_critical_across_α(;αs = LinRange(0.0,2.0,30),β=0.5,h=1.0,iterations=1000,voltage_saveloc="./alpha_vary_critical_plot_b5_volt.png",rate_saveloc="./alpha_vary_critical_plot_b5_rate.png")
    #each dict will correspond to the index of true/false 
    voltage_fixed_point_arrs = []
    rate_fixed_point_arrs = []
    stability_arrs = []
    is_stable_arrs = []
    max_stability = 0
    max_instability = 0
    #finds max number of stable points and unstable points
    for α ∈ αs
        voltage_fixed_points,stability,is_stable,rate_fixed_points = find_fixed_points(α=α,β=β,h=h,iterations=iterations)
        stable_points = sum(is_stable)
        unstable_points = length(is_stable) - stable_points
        if stable_points > max_stability
            max_stability = stable_points
        end
        if unstable_points > max_instability
            max_instability = unstable_points
        end
        push!(voltage_fixed_point_arrs,voltage_fixed_points)
        push!(stability_arrs,stability)
        push!(is_stable_arrs,is_stable)
        push!(rate_fixed_point_arrs,rate_fixed_points)
    end
    #constructs point dict for scatter plots
    stable_voltage_dicts = [Dict() for _ in 1:max_stability]
    unstable_voltage_dicts = [Dict() for _ in 1:max_instability]
    stable_rate_dicts = [Dict() for _ in 1:max_stability]
    unstable_rate_dicts = [Dict() for _ in 1:max_instability]
    for i in eachindex(αs)
        stable_counter = 1
        unstable_counter = 1
        for j in eachindex(is_stable_arrs[i])
            if is_stable_arrs[i][j]
                stable_voltage_dicts[stable_counter][αs[i]] = voltage_fixed_point_arrs[i][j]
                stable_rate_dicts[stable_counter][αs[i]] = rate_fixed_point_arrs[i][j]
                stable_counter += 1
            else
                unstable_voltage_dicts[unstable_counter][αs[i]] = voltage_fixed_point_arrs[i][j]
                unstable_rate_dicts[unstable_counter][αs[i]] = rate_fixed_point_arrs[i][j]
                unstable_counter += 1
            end
        end
    end
    #plots voltage
    r1e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=1,color=:red,xlabel="α",ylabel="Critical Point Voltage")
    r1i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=2,color=:blue,xlabel="α",ylabel="Critical Point Voltage")
    r2e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=3,color=:orange,xlabel="α",ylabel="Critical Point Voltage")
    r2i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=4,color=:purple,xlabel="α",ylabel="Critical Point Voltage")
    v_sub_plots = plot(r1e_v_stability_plot,r1i_v_stability_plot,r2e_v_stability_plot,r2i_v_stability_plot,layout=(2,2))
    beginning_alpha = αs[begin]
    ending_alpha = αs[end]
    combined_plots = plot(v_sub_plots,plot_title="Critical Points Voltage Across α ∈ [$beginning_alpha,$ending_alpha] with fixed β=$β,h=$h",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,voltage_saveloc)
    #plots rate
    r1e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=1,color=:red,xlabel="α",ylabel="Critical Point Rate")
    r1i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=2,color=:blue,xlabel="α",ylabel="Critical Point Rate")
    r2e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=3,color=:orange,xlabel="α",ylabel="Critical Point Rate")
    r2i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=4,color=:purple,xlabel="α",ylabel="Critical Point Rate")
    r_sub_plots = plot(r1e_r_stability_plot,r1i_r_stability_plot,r2e_r_stability_plot,r2i_r_stability_plot,layout=(2,2))
    beginning_alpha = αs[begin]
    ending_alpha = αs[end]
    combined_plots = plot(r_sub_plots,plot_title="Critical Points Rate Across α ∈ [$beginning_alpha,$ending_alpha] with fixed β=$β,h=$h",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,rate_saveloc)
end

function plot_critical_across_ei(;eis = LinRange(0.1,1.5,30),β=0.5,h=1.0,iterations=1000,voltage_saveloc="./ei_vary_critical_plot_b5_volt.png",rate_saveloc="./ei_vary_critical_plot_b5_rate.png")
    #each dict will correspond to the index of true/false 
    voltage_fixed_point_arrs = []
    rate_fixed_point_arrs = []
    stability_arrs = []
    is_stable_arrs = []
    max_stability = 0
    max_instability = 0
    #finds max number of stable points and unstable points
    for ei ∈ eis
        voltage_fixed_points,stability,is_stable,rate_fixed_points = find_fixed_points(α=ei*β,β=β,h=h,iterations=iterations)
        stable_points = sum(is_stable)
        unstable_points = length(is_stable) - stable_points
        if stable_points > max_stability
            max_stability = stable_points
        end
        if unstable_points > max_instability
            max_instability = unstable_points
        end
        push!(voltage_fixed_point_arrs,voltage_fixed_points)
        push!(stability_arrs,stability)
        push!(is_stable_arrs,is_stable)
        push!(rate_fixed_point_arrs,rate_fixed_points)
    end
    #constructs point dict for scatter plots
    stable_voltage_dicts = [Dict() for _ in 1:max_stability]
    unstable_voltage_dicts = [Dict() for _ in 1:max_instability]
    stable_rate_dicts = [Dict() for _ in 1:max_stability]
    unstable_rate_dicts = [Dict() for _ in 1:max_instability]
    for i in eachindex(eis)
        stable_counter = 1
        unstable_counter = 1
        for j in eachindex(is_stable_arrs[i])
            if is_stable_arrs[i][j]
                stable_voltage_dicts[stable_counter][eis[i]] = voltage_fixed_point_arrs[i][j]
                stable_rate_dicts[stable_counter][eis[i]] = rate_fixed_point_arrs[i][j]
                stable_counter += 1
            else
                unstable_voltage_dicts[unstable_counter][eis[i]] = voltage_fixed_point_arrs[i][j]
                unstable_rate_dicts[unstable_counter][eis[i]] = rate_fixed_point_arrs[i][j]
                unstable_counter += 1
            end
        end
    end
    #plots voltage
    r1e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=1,color=:red,xlabel="E/I",ylabel="Critical Point Voltage")
    r1i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=2,color=:blue,xlabel="E/I",ylabel="Critical Point Voltage")
    r2e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=3,color=:orange,xlabel="E/I",ylabel="Critical Point Voltage")
    r2i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=4,color=:purple,xlabel="E/I",ylabel="Critical Point Voltage")
    v_sub_plots = plot(r1e_v_stability_plot,r1i_v_stability_plot,r2e_v_stability_plot,r2i_v_stability_plot,layout=(2,2))
    beginning_ei = eis[begin]
    ending_ei = eis[end]
    combined_plots = plot(v_sub_plots,plot_title="Critical Points Voltage Across E/I ∈ [$beginning_ei,$ending_ei] with fixed β=$β,h=$h",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,voltage_saveloc)
    #plots rate
    r1e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=1,color=:red,xlabel="E/I",ylabel="Critical Point Rate")
    r1i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=2,color=:blue,xlabel="E/I",ylabel="Critical Point Rate")
    r2e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=3,color=:orange,xlabel="E/I",ylabel="Critical Point Rate")
    r2i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=4,color=:purple,xlabel="E/I",ylabel="Critical Point Rate")
    r_sub_plots = plot(r1e_r_stability_plot,r1i_r_stability_plot,r2e_r_stability_plot,r2i_r_stability_plot,layout=(2,2))
    combined_plots = plot(r_sub_plots,plot_title="Critical Points Rate Across E/I ∈ [$beginning_ei,$ending_ei] with fixed β=$β,h=$h",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,rate_saveloc)
end

function plot_critical_across_h(;ei = 1.0,β=0.5,hs=LinRange(0.0,10.0,30),iterations=1000,voltage_saveloc="./h_vary_critical_plot_b5_volt.png",rate_saveloc="./h_vary_critical_plot_b5_rate.png")
    #each dict will correspond to the index of true/false 
    voltage_fixed_point_arrs = []
    rate_fixed_point_arrs = []
    stability_arrs = []
    is_stable_arrs = []
    max_stability = 0
    max_instability = 0
    #finds max number of stable points and unstable points
    for h ∈ hs
        voltage_fixed_points,stability,is_stable,rate_fixed_points = find_fixed_points(α=β*ei,β=β,h=h,iterations=iterations)
        stable_points = sum(is_stable)
        unstable_points = length(is_stable) - stable_points
        if stable_points > max_stability
            max_stability = stable_points
        end
        if unstable_points > max_instability
            max_instability = unstable_points
        end
        push!(voltage_fixed_point_arrs,voltage_fixed_points)
        push!(stability_arrs,stability)
        push!(is_stable_arrs,is_stable)
        push!(rate_fixed_point_arrs,rate_fixed_points)
    end
    #constructs point dict for scatter plots
    stable_voltage_dicts = [Dict() for _ in 1:max_stability]
    unstable_voltage_dicts = [Dict() for _ in 1:max_instability]
    stable_rate_dicts = [Dict() for _ in 1:max_stability]
    unstable_rate_dicts = [Dict() for _ in 1:max_instability]
    for i in eachindex(hs)
        stable_counter = 1
        unstable_counter = 1
        for j in eachindex(is_stable_arrs[i])
            if is_stable_arrs[i][j]
                stable_voltage_dicts[stable_counter][hs[i]] = voltage_fixed_point_arrs[i][j]
                stable_rate_dicts[stable_counter][hs[i]] = rate_fixed_point_arrs[i][j]
                stable_counter += 1
            else
                unstable_voltage_dicts[unstable_counter][hs[i]] = voltage_fixed_point_arrs[i][j]
                unstable_rate_dicts[unstable_counter][hs[i]] = rate_fixed_point_arrs[i][j]
                unstable_counter += 1
            end
        end
    end
    #plots voltage
    r1e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=1,color=:red,xlabel="h",ylabel="Critical Point Voltage")
    r1i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=2,color=:blue,xlabel="h",ylabel="Critical Point Voltage")
    r2e_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=3,color=:orange,xlabel="h",ylabel="Critical Point Voltage")
    r2i_v_stability_plot = plot_quadrant(stable_voltage_dicts,unstable_voltage_dicts;quadrant=4,color=:purple,xlabel="h",ylabel="Critical Point Voltage")
    v_sub_plots = plot(r1e_v_stability_plot,r1i_v_stability_plot,r2e_v_stability_plot,r2i_v_stability_plot,layout=(2,2))
    beginning_h = hs[begin]
    ending_h = hs[end]
    combined_plots = plot(v_sub_plots,plot_title="Critical Points Voltage Across h ∈ [$beginning_h,$ending_h] with fixed β=$β,E/I=$(α/β)",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,voltage_saveloc)
    #plots rate
    r1e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=1,color=:red,xlabel="h",ylabel="Critical Point Rate")
    r1i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=2,color=:blue,xlabel="h",ylabel="Critical Point Rate")
    r2e_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=3,color=:orange,xlabel="h",ylabel="Critical Point Rate")
    r2i_r_stability_plot = plot_quadrant(stable_rate_dicts,unstable_rate_dicts;quadrant=4,color=:purple,xlabel="h",ylabel="Critical Point Rate")
    r_sub_plots = plot(r1e_r_stability_plot,r1i_r_stability_plot,r2e_r_stability_plot,r2i_r_stability_plot,layout=(2,2))
    combined_plots = plot(r_sub_plots,plot_title="Critical Points Rate Across h ∈ [$beginning_h,$ending_h] with fixed β=$β,E/I=$(α/β)",
    thickness_scaling=1,titlefontsize=10,size=(700,600))
    savefig(combined_plots,rate_saveloc)
end

function plot_critical_count_across_α_h(;β=1.1,hs=LinRange(0.0,20.0,50),αs=LinRange(0.0,2.0,20),stable_saveloc="./stable_as_hs_plots.png",unstable_saveloc="./unstable_as_hs_plots.png")
    stable_count = Array{Int}(undef,length(hs),length(αs))
    unstable_count = Array{Int}(undef,length(hs),length(αs))
    for hi ∈ ProgressBar(eachindex(hs))
        for αi ∈ eachindex(αs)
            _,_,is_stable = find_fixed_points(α=αs[αi],β=β,h=hs[hi],iterations=10000)
            stable_equilibria = sum(is_stable)
            unstable_equilibria = length(is_stable) - stable_equilibria
            stable_count[hi,αi] = stable_equilibria
            unstable_count[hi,αi] = unstable_equilibria
        end
    end
    stable_color_pal = palette([:black,:cyan],max(stable_count...)+1)
    stable_plot = heatmap(αs,hs,stable_count;c=stable_color_pal,xaxis="α",yaxis="h",colorbar=false,title="Number of Stable Critical Points With Fixed β=$β")
    for i in 1:max(stable_count...)+1
        plot!(stable_plot,[], [], seriestype=:shape,label=i-1, color=stable_color_pal[i])
    end

    unstable_color_pal = palette([:black,:fuchsia],max(unstable_count...)+1)
    unstable_plot = heatmap(αs,hs,unstable_count;c=unstable_color_pal,xaxis="α",yaxis="h",colorbar=false,title="Number of Unstable Critical Points With Fixed β=$β")
    for i in 1:max(unstable_count...)+1
        plot!(unstable_plot,[], [], seriestype=:shape,label=i-1, color=unstable_color_pal[i])
    end
    savefig(stable_plot,stable_saveloc)
    savefig(unstable_plot,unstable_saveloc)
end

function plot_critical_number_across_ei_h(;β=1.1,hs=LinRange(0.0,20.0,200),ei_ratios=LinRange(0.1,2.0,200),max_critical_points=2,stable_saveloc="./stable_ei_hs_plots.png",unstable_saveloc="./unstable_ei_hs_plots.png")
    display("Running $(Threads.nthreads()) Threads")
    stable_count = Array{Int}(undef,length(hs),length(ei_ratios))
    unstable_count = Array{Int}(undef,length(hs),length(ei_ratios))
    Threads.@threads for hi ∈ ProgressBar(eachindex(hs))
        for ei_ratios_i ∈ eachindex(ei_ratios)
            _,_,is_stable = find_fixed_points(α=ei_ratios[ei_ratios_i]*β,β=β,h=hs[hi],iterations=1000)
            stable_equilibria = sum(is_stable)
            unstable_equilibria = length(is_stable) - stable_equilibria
            stable_count[hi,ei_ratios_i] = stable_equilibria
            unstable_count[hi,ei_ratios_i] = unstable_equilibria
        end
    end
    stable_color_pal = palette([:black,:cyan],max_critical_points+1)
    unstable_color_pal = palette([:black,:fuchsia],max_critical_points+1)
    numbers_in_stable = sort([x for x in Set([stable_count...])])
    numbers_in_unstable = sort([x for x in Set([unstable_count...])])
    display_stable_color_pal = [stable_color_pal[c+1] for c in numbers_in_stable]
    display_unstable_color_pal = [unstable_color_pal[c+1] for c in numbers_in_unstable]
    if max(stable_count...) > max_critical_points
        display("WARNING: Stable Critical Points Exceed Max")
    end
    if max(unstable_count...) > max_critical_points
        display("WARNING: Unstable Critical Points Exceed Max")
    end
    stable_plot = heatmap(ei_ratios,hs,stable_count;c=display_stable_color_pal,xaxis="E/I Ratio",yaxis="h",colorbar=false,title="Number of Stable Critical Points With Fixed β=$β")
    for i in 1:max_critical_points+1
        plot!(stable_plot,[], [], seriestype=:shape,label=i-1, color=stable_color_pal[i])
    end
    unstable_plot = heatmap(ei_ratios,hs,unstable_count;c=display_unstable_color_pal,xaxis="E/I Ratio",yaxis="h",colorbar=false,title="Number of Unstable Critical Points With Fixed β=$β")
    for i in 1:max_critical_points+1
        plot!(unstable_plot,[], [], seriestype=:shape,label=i-1, color=unstable_color_pal[i])
    end
    savefig(stable_plot,stable_saveloc)
    savefig(unstable_plot,unstable_saveloc)
end
#plots for β=0.5
#β = 0.5
#plot_critical_across_ei(β=β,h=0.0,voltage_saveloc="./critical_plots/b05/EI_critical_voltage_h10.png",rate_saveloc="./critical_plots/b05/EI_critical_rate_h10.png")
#plot_critical_across_h(β=β,α=1.1,voltage_saveloc="./critical_plots/b05/h_critical_voltage_α05.png",rate_saveloc="./critical_plots/b05/h_critical_rate_α05.png")
#plot_critical_number_across_ei_h(β=β,hs = LinRange(0.0, 20.0, 200), ei_ratios = LinRange(0.1, 2.0, 200),max_critical_points = 2, stable_saveloc = "./critical_plots/b05/stable_ei_hs_plots.png", unstable_saveloc = "./critical_plots/b05/unstable_ei_hs_plots.png")
#plots for β = 1.1
β = 1.1
#plot_critical_across_ei(β=β,h=0.0,voltage_saveloc="./critical_plots/EI_critical_voltage_h10.png",rate_saveloc="./critical_plots/EI_critical_rate_h10.png")
plot_critical_across_h(β=β,α=1.32,voltage_saveloc="./critical_plots/h_critical_voltage_α132.png",rate_saveloc="./critical_plots/h_critical_rate_α132.png")
plot_critical_across_h(β=β,α=1.65,voltage_saveloc="./critical_plots/h_critical_voltage_α165.png",rate_saveloc="./critical_plots/h_critical_rate_α165.png")

#plot_critical_number_across_ei_h(β=β,hs = LinRange(0.0, 20.0, 200), ei_ratios = LinRange(0.1, 2.0, 200),max_critical_points = 2, stable_saveloc = "./critical_plots/b11/stable_ei_hs_plots.png", unstable_saveloc = "./critical_plots/b11/unstable_ei_hs_plots.png")
#plots for β = 1.5
#β = 1.5
#plot_critical_across_ei(β=β,h=0.0,voltage_saveloc="./critical_plots/b15/EI_critical_voltage_h10.png",rate_saveloc="./critical_plots/b15/EI_critical_rate_h10.png")
#plot_critical_across_h(β=β,α=1.1,voltage_saveloc="./critical_plots/b15/h_critical_voltage_α11.png",rate_saveloc="./critical_plots/b15/h_critical_rate_α15.png")
#plot_critical_number_across_ei_h(β=β,hs = LinRange(0.0, 20.0, 200), ei_ratios = LinRange(0.1, 2.0, 200),max_critical_points = 2, stable_saveloc = "./critical_plots/b15/stable_ei_hs_plots.png", unstable_saveloc = "./critical_plots/b15/unstable_ei_hs_plots.png")
#plot_critical_across_ei()
#plot_across_ei_h(β=0.5,stable_saveloc="./stable_ei_hs_plots_b05.png",unstable_saveloc="./unstable_ei_hs_plots_b05.png")
#plot_across_ei_h(β=0.5,stable_saveloc="./critical_plots/stable_ei_hs_plots_b05.png",unstable_saveloc="./critical_plots/unstable_ei_hs_plots_b05.png")
#plot_across_ei_h(β=1.1,stable_saveloc="./critical_plots/stable_ei_hs_plots_b11.png",unstable_saveloc="./critical_plots/unstable_ei_hs_plots_b11.png")
#plot_across_ei_h(β=1.5,stable_saveloc="./critical_plots/stable_ei_hs_plots_b15.png",unstable_saveloc="./critical_plots/unstable_ei_hs_plots_b15.png")

#plot_across_h(α=0.5,β=1.1,saveloc="./critical_plots/h_vary_critical_plot_a05.png")
#plot_across_h(α=1.1,β=1.1,saveloc="./critical_plots/h_vary_critical_plot_a10.png")
#plot_across_h(α=1.5,β=1.1,saveloc="./critical_plots/h_vary_critical_plot_a15.png")

#plot_across_α(β=0.5,saveloc="./critical_plots/alpha_vary_critical_plot_b05.png")
#plot_across_α(β=1.1,h=1.0,saveloc="./critical_plots/alpha_vary_critical_plot_b11_h10.png")
#plot_across_α(β=1.1,h=0.5,saveloc="./critical_plots/alpha_vary_critical_plot_b11_h05.png")
#plot_across_α(β=1.1,h=5.0,saveloc="./critical_plots/alpha_vary_critical_plot_b11_h50.png")