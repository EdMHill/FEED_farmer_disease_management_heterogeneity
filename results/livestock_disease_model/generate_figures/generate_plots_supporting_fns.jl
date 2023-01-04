#=
Supporting functions to generate the plots

Julia version: 1.8.3
=#


#========================
Tick length modifier function
========================#
function ticks_length!(; tl=0.02)
    p = Plots.current()
    xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
    xl, yl = Plots.xlims(p), Plots.ylims(p)
    x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
    sz = p.attr[:size]
    r = sz[1] / sz[2]
    dx, dy = tl * (xl[2] - xl[1]), tl * r * (yl[2] - yl[1])
    plot!([xticks xticks]', [y1 y1 .+ dy]', c=:black, labels=false)
    plot!([x1 x1 .+ dx]', [yticks yticks]', c=:black, labels=false, xlims=xl, ylims=yl)
    return Plots.current()
end


#========================
Spatial map plotting functions
========================#

function plot_spatial_map_fn(input_event_array::Union{Array{Float64,3},Array{Int64,3}},
    quantile_vals_spatial_maps::Vector{Float64},
    n_quantile_vals::Int64,
    n_seed_region_scenarios::Int64,
    n_behav_group_configs::Int64,
    n_counties_UAs::Int64,
    region_seed_infection_scenario_ID::Vector{Int64},
    plot_title_vec::Vector{String},
    cmap::Union{Symbol,PlotUtils.ContinuousColorGradient},
    clim_vals::Tuple{Float64,Float64},
    colorbar_string::String,
    plot_fontsize::Int64,
    save_fig_flag::Bool,
    save_filename_prefix::String,
    fig_width::Int64,
    fig_height::Int64)

    # Initialise arrays to store percentile results & compute specified percentiles
    quantiles_array = zeros(Float64, n_quantile_vals, n_seed_region_scenarios, n_behav_group_configs)
    for behav_config_idx = 1:n_behav_group_configs
        for seed_infection_scen_idx = 1:n_seed_region_scenarios

            # Percentiles based on linked regions for seeding infections
            quantiles_array[:, seed_infection_scen_idx, behav_config_idx] =
                quantile(input_event_array[:, seed_infection_scen_idx, behav_config_idx], quantile_vals_spatial_maps)
        end
    end

    # Map seed region locations to individual county/unitary authority IDs
    # Uses region_seed_infection_scenario_ID
    quantiles_by_county_UA_array = zeros(Float64, n_quantile_vals, n_counties_UAs, n_behav_group_configs)
    for county_UA_ID = 1:n_counties_UAs

        # Check if county/UA is in landscape and has cattle holdings
        # If not, assign quantile values as -1
        if region_seed_infection_scenario_ID[county_UA_ID] == 0
            quantiles_by_county_UA_array[:, county_UA_ID, :] .= -1
        else
            quantiles_by_county_UA_array[:, county_UA_ID, :] .= quantiles_array[:, region_seed_infection_scenario_ID[county_UA_ID], :]
        end
    end

    println("quantiles_by_county_UA_array minimum value: $(minimum(quantiles_by_county_UA_array))")
    println("quantiles_by_county_UA_array maximum value: $(maximum(quantiles_by_county_UA_array))")

    # Initialise plot vector
    p_array = Matrix{Plots.Plot{Plots.GRBackend}}(undef, n_behav_group_configs, n_quantile_vals)

    # Generate the figures
    for behav_config_idx = 1:n_behav_group_configs

        # Set title for plot (if applicable)
        if length(plot_title_vec) == n_behav_group_configs
            plot_title = plot_title_vec[behav_config_idx]
        end

        for quantile_idx = 1:n_quantile_vals
            fill_z_data = quantiles_by_county_UA_array[quantile_idx, [1:152; 164:189; 191:217], behav_config_idx]
            println("fill_z_data minimum & maximum: $(minimum(fill_z_data)),$(maximum(fill_z_data))")

            # Set title for plot (if applicable)
            if length(plot_title_vec) == n_quantile_vals
                plot_title = plot_title_vec[quantile_idx]
            end

            p_array[behav_config_idx, quantile_idx] =
                plot(county_UA_polygons[[1:152; 164:189; 191:217]],
                    fill_z=fill_z_data', # Shade each polygon based on fill_z vector (input needs to be a row vector)
                    title=plot_title,
                    cbar=true,
                    clim=clim_vals,
                    colorbar_title=colorbar_string,
                    colorbar_titlefontsize=plot_fontsize,
                    color=cmap,
                    grid=false,
                    showaxis=false,
                    xticks=nothing,
                    yticks=nothing,
                    xtickfontsize=plot_fontsize,
                    ytickfontsize=plot_fontsize,
                    guidefontsize=plot_fontsize,
                    framestyle=:box,
                    left_margin=2Plots.mm,
                    right_margin=5Plots.mm,
                    aspect_ratio=1.0,
                    xlims=(0, 700000))

            #Specify the overall plot size
            plot!(size=(fig_width, fig_height))

            # If applicable, save the figures
            if save_fig_flag == true
                savefig(p_array[behav_config_idx, quantile_idx], string(save_filename_prefix, "_config_ID$(behav_config_idx)_quantile_$(quantile_vals_spatial_maps[quantile_idx]).pdf"))
            end
        end
    end

    return nothing
end

function plot_outbreak_threshold_spatial_map_fn(input_event_array::Union{Array{Float64,3},Array{Int64,3}},
    threshold_vals_vec::Union{Vector{Int64},Vector{Float64}},
    n_outbreak_threshold_sizes::Int64,
    n_seed_region_scenarios::Int64,
    n_behav_group_configs::Int64,
    n_counties_UAs::Int64,
    region_seed_infection_scenario_ID::Vector{Int64},
    plot_title_vec::Vector{String},
    cmap::Union{Symbol,PlotUtils.ContinuousColorGradient},
    clim_vals::Tuple{Float64,Float64},
    colourbar_string_prefix::String,
    colourbar_string_suffix::String,
    plot_fontsize::Int64,
    save_fig_flag::Bool,
    save_filename_prefix::String,
    fig_width::Int64,
    fig_height::Int64)

    # Map seed region locations to individual county/unitary authority IDs
    # Uses region_seed_infection_scenario_ID
    input_data_by_county_UA_array = zeros(Float64, n_counties_UAs, n_behav_group_configs, n_outbreak_threshold_sizes)
    for county_UA_ID = 1:n_counties_UAs

        # Check if county/UA is in landscape and has cattle holdings
        # If not, assign quantile values as -1
        if region_seed_infection_scenario_ID[county_UA_ID] == 0
            input_data_by_county_UA_array[county_UA_ID, :, :] .= -1
        else
            println(size(input_event_array))
            println(region_seed_infection_scenario_ID[county_UA_ID])
            println(size(input_event_array[region_seed_infection_scenario_ID[county_UA_ID], :, :]))
            input_data_by_county_UA_array[county_UA_ID, :, :] = 1 .* input_event_array[region_seed_infection_scenario_ID[county_UA_ID], :, :]
        end
    end

    println("input_data_by_county_UA_array minimum value: $(minimum(input_data_by_county_UA_array))")
    println("input_data_by_county_UA_array maximum value: $(maximum(input_data_by_county_UA_array))")

    # Initialise plot vector
    p_array = Matrix{Plots.Plot{Plots.GRBackend}}(undef, n_behav_group_configs, n_outbreak_threshold_sizes)

    # Generate the figures
    for behav_config_idx = 1:n_behav_group_configs

        # Set title for plot (if applicable)
        if length(plot_title_vec) == n_behav_group_configs
            plot_title = plot_title_vec[behav_config_idx]
        end

        for outbreak_threshold_idx = 1:n_outbreak_threshold_sizes
            fill_z_data = input_data_by_county_UA_array[[1:152; 164:189; 191:217], behav_config_idx, outbreak_threshold_idx]
            println("fill_z_data minimum & maximum: $(minimum(fill_z_data)),$(maximum(fill_z_data))")

            # Modify colorbar_string
            println(threshold_vals_vec)
            println(threshold_vals_vec[outbreak_threshold_idx])
            colourbar_string = string(colourbar_string_prefix, " > $(threshold_vals_vec[outbreak_threshold_idx])",colourbar_string_suffix, " (% of simns)")

            # Construct the plot
            p_array[behav_config_idx, outbreak_threshold_idx] =
                plot(county_UA_polygons[[1:152; 164:189; 191:217]],
                    fill_z=100 * fill_z_data', # Shade each polygon based on fill_z vector (input needs to be a row vector). Multiply by 100 to convert proportion to a %
                    title=plot_title,
                    cbar=true,
                    clim=clim_vals,
                    colorbar_title=colourbar_string,
                    colorbar_titlefontsize=plot_fontsize,
                    color=cmap,
                    grid=false,
                    showaxis=false,
                    xticks=nothing,
                    yticks=nothing,
                    xtickfontsize=plot_fontsize,
                    ytickfontsize=plot_fontsize,
                    guidefontsize=plot_fontsize,
                    framestyle=:box,
                    left_margin=2Plots.mm,
                    right_margin=5Plots.mm,
                    aspect_ratio=1.0,
                    xlims=(0, 700000))

            #Specify the overall plot size
            plot!(size=(fig_width, fig_height))

            # If applicable, save the figures
            if save_fig_flag == true
                savefig(p_array[behav_config_idx, outbreak_threshold_idx], string(save_filename_prefix, "_config_ID$(behav_config_idx)_threshold_value_$(threshold_vals_vec[outbreak_threshold_idx]).pdf"))
            end
        end
    end

    return nothing
end


#========================
Violin plot
========================#
function plot_violin_figs_fn(input_event_array::Union{Array{Float64,2},Array{Int64,2}},
    uncertainty_interval_quantile_vals::Vector{Float64},
    n_behav_group_configs::Int64,
    x_vals::Vector{Int64},
    config_xtick_label::Vector{String},
    y_axis_label::String,
    y_lim_max_val::Float64,
    plot_fontsize::Int64,
    save_fig_flag::Bool,
    save_filename_violin_plot::String,
    fig_width::Float64,
    fig_height::Float64)

    # Generate violin plots (one at a time, iterate over each column, enables
    # specifying the x-position of the violin plot)
    p_violins = violin([x_vals[1]],
        input_event_array[:, 1],
        xlabel="Behavioural configuration",
        xticks=(x_vals, config_xtick_label),
        xrotation=30,
        ylabel=y_axis_label,
        formatter=:plain,
        color=theme_palette(:auto).colors.colors[1],
        legend=false,
        grid=false,
        framestyle=:box
    )

    for ii = 2:length(x_vals)
        violin!(p_violins,
            [x_vals[ii]],
            input_event_array[:, ii],
            color=theme_palette(:auto).colors.colors[1])
    end

    # Add median point
    scatter!(p_violins, x_vals, vec(median(input_event_array,dims=1)), markersize=2, color=:red)

    # Add divider lines to separate out configs that are varying similar attributes
    line_divide_x_vals = [2; 6; 9]
    for add_dividing_line_itr = 1:length(line_divide_x_vals)
        plot!(p_violins,
            [line_divide_x_vals[add_dividing_line_itr]; line_divide_x_vals[add_dividing_line_itr]], # x position
            [-0.035 * y_lim_max_val; y_lim_max_val * 1.1],   # y position
            linestyle=:dash,
            linecolor=:gray,
            linewidth=1.5,
            ylims=(-0.035 * y_lim_max_val, y_lim_max_val)
        )
    end

    # Set up margins
    plot!(left_margin=2Plots.mm,
        right_margin=2Plots.mm,
        bottom_margin=5Plots.mm,
        dpi=600)

    # If applicable, save the figures
    if save_fig_flag == true
        savefig(p_violins, string(save_filename_violin_plot, ".pdf"))
        savefig(p_violins, string(save_filename_violin_plot, ".png"))
    end
end


#========================
Bar plot: Proportion of simulations that exceed a given threshold value
========================#
function generate_threshold_grouped_bar_fig(results_array::Union{Matrix{Int64},Matrix{Float64}},
    threshold_vals_vec::Union{Vector{Int64},Vector{Float64}},
    n_simns_per_config::Int64,
    x_vals::Vector{Int64},
    config_xtick_label::Vector{String},
    yaxis_label::String,
    legend_position_vals::Tuple{Float64,Float64},
    legend_labels::Matrix{String},
    bar_fillcolour::Matrix{RGB{Float64}},
    save_fig_flag::Bool,
    save_filename::String)

    # Get number of threshold values input
    n_thresholds = length(threshold_vals_vec)

    # Error checks
    if length(legend_labels) != n_thresholds
        error("Number of threshold values: $n_thresholds. Number of legend labels: $(length(legend_labels)). Number of entries in each should be equal.")
    end

    if length(bar_fillcolour) != n_thresholds
        error("Number of threshold values: $n_thresholds. Number of bar fill colours: $(length(bar_fillcolour)). Number of entries in each should be equal.")
    end

    # Initialise the arrays to store the proportion & quantile values
    n_configs = size(results_array, 2)
    propn_outbreaks_satisfying_threshold = zeros(Float64, n_configs, n_thresholds)

    # Iterate over each thresholdvalue
    # Get counts & proportion of simulations that satisfy the threshold condition
    for threshold_val_itr = 1:n_thresholds

        # Get count of simulations that were shorter than the specified outbreak duration
        threshold_check_result_array = results_array .> threshold_vals_vec[threshold_val_itr]
        count_outbreaks_satisfying_threshold = vec(sum(threshold_check_result_array, dims=1))

        # Get proportion of simulations that were shorter than the specified outbreak duration
        propn_outbreaks_satisfying_threshold[:, threshold_val_itr] = count_outbreaks_satisfying_threshold ./ n_simns_per_config

        # # Compute the binomial proportion confidence intervals
        # # Based on jeffreys interval:
        # # After observing x successes in n trials, the posterior distribution for p
        # # is a Beta distribution with parameters (x + 1/2, n – x + 1/2).
        # quantiles_propn_outbreaks_shorter_than_duration_threshold[:,:,outbreak_duration_itr] =
        #         quantile.(Beta.(count_outbreaks_shorter_than_duration_threshold .+ 0.5,
        #                         n_simns_per_config .- count_outbreaks_shorter_than_duration_threshold .+ 0.5), quantile_vals)
    end

    # Generate grouped bar plot
    p_threshold_bar = groupedbar(x_vals,
        propn_outbreaks_satisfying_threshold * 100,
        bar_width=0.7,
        color=bar_fillcolour,
        xlabel="Behavioural configuration",
        xticks=(x_vals, config_xtick_label),
        xrotation=30,
        ylabel=yaxis_label,
        formatter=:plain,
        legend=true,
        label=legend_labels,
        #legend_position = legend_position_vals,
        legend_column=-1, # Specify horizontal orientation with input to egend_column of -1
        grid=false,
        framestyle=:box,
        left_margin=2Plots.mm,
        right_margin=2Plots.mm,
        bottom_margin=5Plots.mm,
        dpi=600)

    # Add divider lines to separate out configs that are varying similar attributes
    line_divide_x_vals = [2; 6; 9]
    for add_dividing_line_itr = 1:length(line_divide_x_vals)
        plot!(p_threshold_bar,
            [line_divide_x_vals[add_dividing_line_itr]; line_divide_x_vals[add_dividing_line_itr]], # x position
            [-3.5; 105],   # y position
            ylims=(0, 101),
            linestyle=:dash,
            linecolor=:gray,
            linewidth=1.5,
            label="",
        )
    end

    # Set legend position
    plot!(legend_position=:outertop)

    # If applicable, save the figures
    if save_fig_flag == true
        savefig(p_threshold_bar, string(save_filename, ".pdf"))
        savefig(p_threshold_bar, string(save_filename, ".png"))
    end

    return nothing
end


#========================
Compute the proportion of simulations satisfying threshold criteria
========================#
function get_propn_simn_satisfying_threshold_summ_stats(reordered_input_array::Union{Array{Int64,3},Array{Float64,3}},
                                                            threshold_vals_vec::Union{Vector{Int64},Vector{Float64}},
                                                            n_threshold_vals::Int64,
                                                            n_seed_region_scenarios::Int64,
                                                            n_configs::Int64,
                                                            n_simn_replicates::Int64,
                                                            quantile_vals::Vector{Float64},
                                                            n_quantile_vals::Int64)

    # Initialise the array to store the proportion of simulation where outbreak size
    # did not exceed the stipulated threshold
    propn_outbreaks_satisfying_threshold = zeros(Float64, n_seed_region_scenarios, n_configs, n_threshold_vals)

    # Initialise the array to store the proportion & quantile values
    quantiles_propn_outbreaks_satisfying_threshold = zeros(Float64, n_seed_region_scenarios, n_configs, n_quantile_vals, n_threshold_vals)

    # Iterate over threshold values
    # Get counts & proportion of simulations that are within the stipulated threshold value
    for threshold_val_itr = 1:n_threshold_vals

        # Get count of simulations that were contained below the specified outbreak size
        outbreak_threshold_result = reordered_input_array .> threshold_vals_vec[threshold_val_itr]
        count_outbreaks_satisfying_threshold = dropdims(sum(outbreak_threshold_result, dims=1); dims=1)

        # Get proportion of simulations that were contained below the specified outbreak size
        propn_outbreaks_satisfying_threshold[:, :, threshold_val_itr] = count_outbreaks_satisfying_threshold ./ n_simn_replicates

        # Compute the binomial proportion confidence intervals
        # Based on jeffreys interval:
        # After observing x successes in n trials, the posterior distribution for p
        # is a Beta distribution with parameters (x + 1/2, n – x + 1/2).
        for quantile_itr = 1:n_quantile_vals
            quantiles_propn_outbreaks_satisfying_threshold[:, :, quantile_itr, threshold_val_itr] =
                quantile.(Beta.(count_outbreaks_satisfying_threshold .+ 0.5,
                        n_simn_replicates .- count_outbreaks_satisfying_threshold .+ 0.5), quantile_vals[quantile_itr])
        end
    end

    return propn_outbreaks_satisfying_threshold::Array{Float64,3},
    quantiles_propn_outbreaks_satisfying_threshold::Array{Float64,4}
end