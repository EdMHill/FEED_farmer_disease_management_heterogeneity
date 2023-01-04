#=
Script to generate plots for the paper

Useful Julia plotting resource:
https://nextjournal.com/leandromartinez98/tips-to-create-beautiful-publication-quality-plots-in-julia

Font family information:
https://gr-framework.org/fonts.html

Julia version: 1.8.3
=#

#========================
LOAD REQUIRED ENVIRONMENT
========================#

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../../")

#========================
LOAD REQUIRED PACKAGES
========================#
using DelimitedFiles
using StatsBase
using Plots, StatsPlots
using Shapefile
using ColorSchemes
using DataFrames
using Distributions
using JLD2

#========================
INCLUDE REQUIRED FILES
========================#
include("generate_plots_supporting_fns.jl")

#========================
SET GLOBAL FIGURE PROPERTIES
========================#

# Fontsize
plot_fontsize_default = 12
plot_fontsize_leg = 12
plot_fontsize_violin_plot = 12
plot_fontsize_bar_vals = 6

# Set figure width and height
fig_width = 1.2*550
fig_height = 1.2*450
fig_width_spatial_map = 450
fig_height_spatial_map = 600
fig_width_violin_plot = 1.2*550
fig_height_violin_plot = 1.2*450

# Save variables
save_fig_flag = true

# Shapefile data
shapefile_table =  Shapefile.Table("../../../data/shapefiles/Counties_and_Unitary_Authorities_(December_2021)_UK_BGC/CTYUA_DEC_2021_UK_BGC.shp")
county_UA_polygons = Shapefile.shapes(shapefile_table)
n_counties_UAs = length(county_UA_polygons)

#====================================================
INTERVIEW DATA: HERD SIZE & BEHAVIOURAL CLUSTER GROUPED BAR PLOT (FIGURE YYY)
====================================================#

# Set up the data array
behavioural_group_propns_by_herd_size_empirical_estimates =
        [0.52  0.66  0.20  0.00  0.00;
         0.48  0.29  0.47  0.50  0.00;
         0.00  0.05  0.33  0.50  1.00]'
         # For the plot, group of bars per row (herd size group),
         # behavioural group by column

n_herd_size_groups = size(behavioural_group_propns_by_herd_size_empirical_estimates,1)
n_behavioural_groups = size(behavioural_group_propns_by_herd_size_empirical_estimates,2)

# Set up plot legend labels and bar colours
bar_colour_scheme = palette(:Reds_3,3)
bar_fillcolor = [bar_colour_scheme[1] bar_colour_scheme[2] bar_colour_scheme[3]]
legend_labels = ["Behavioural group 1" "Behavioural group 2" "Behavioural group 3"]

# Generate the grouped bar plot
p_grouped_bar = groupedbar(behavioural_group_propns_by_herd_size_empirical_estimates,
                        xticks = (1:n_herd_size_groups, ["1-100","101-250","251-500","501-1000",">1000"]),
                        bar_width = 0.7,
                        xlabel = "Herd size",
                        ylabel = "Proportion",
                        ylim = (0.,1.05),
                        formatter = :plain,
                        color = bar_fillcolor,
                        label = legend_labels,
                        legend = :topleft,
                        grid = false,
                        framestyle = :box,
                        left_margin = 2Plots.mm,
                        right_margin = 5Plots.mm,
                        dpi = 600)

# Set fontsizes
plot!(xtickfontsize = plot_fontsize_default-2,
        ytickfontsize = plot_fontsize_default-2,
        guidefontsize = plot_fontsize_default)

# Add values to top of each bar
# First sets x-position for the labels, then uses the annotate! function
# to add to the plot.
x_pos = zeros(Float64,n_herd_size_groups,n_behavioural_groups)
for herd_size_grp_itr = 1:n_herd_size_groups
    x_pos[herd_size_grp_itr,:] = (herd_size_grp_itr-1) .+ [0.77 1.0 1.23]
end
for behavioural_grp_itr = 1:n_behavioural_groups
    for herd_size_grp_itr = 1:n_herd_size_groups
        annotate!(x_pos[herd_size_grp_itr,behavioural_grp_itr], # x position
                  behavioural_group_propns_by_herd_size_empirical_estimates[herd_size_grp_itr,behavioural_grp_itr]+0.02, # y position
                  text(string(behavioural_group_propns_by_herd_size_empirical_estimates[herd_size_grp_itr,behavioural_grp_itr]), :black, :centre, plot_fontsize_bar_vals)) # Text to add to plot
    end
end

# Save the figure
if save_fig_flag == true
    savefig(p_grouped_bar, "../../../data/datavis/elicitation_interview_data_figs/grouped_bar_plot_interview_herd_size_behaviour_group_data.pdf")
    savefig(p_grouped_bar, "../../../data/datavis/elicitation_interview_data_figs/grouped_bar_plot_interview_herd_size_behaviour_group_data.png")
end


#========================
IMPORT RESULTS DATA
========================#

# Set up filename prefixes
behav_group_config_disease_state_filename = "../GB_model_epi_outputs_aggregated/holdings_per_disease_state_"
behav_group_config_outbreak_duration_filename = "../GB_model_epi_outputs_aggregated/outbreak_duration_"

# Set up array construction variables
n_simn_replicates = 500
n_seed_region_scenarios = 89
n_behav_group_configs = 7

# Initialise arrays to consolidate results from different scenarios
total_infected_holdings_results_array = zeros(Int64,n_simn_replicates,n_seed_region_scenarios,n_behav_group_configs)
outbreak_duration_results_array = zeros(Int64,n_simn_replicates,n_seed_region_scenarios,n_behav_group_configs)

# # Specify initial batch ID value that batch_ID_offset is then applied to
# batch_baseline_val = 2000

# Specify values for batch_ID_offset
batch_ID_offset_vec = [2000;4100;4200;4300;4400;4500;4600]

# Iterate over behavioural group configs & specified batch IDs.
for behav_config_idx = 1:n_behav_group_configs

    # Specify batch ID offset for this behavioural config
    batch_ID_offset = batch_ID_offset_vec[behav_config_idx]
    # batch_ID_offset = batch_baseline_val + 100*(behav_config_idx-1)

    # For each batch, read data and assign to results arrays
    for seed_infection_scen_idx = 1:n_seed_region_scenarios
        current_disease_state_filename = string(behav_group_config_disease_state_filename,"batchID$(batch_ID_offset+seed_infection_scen_idx).txt")
        current_outbreak_duration_filename = string(behav_group_config_outbreak_duration_filename,"batchID$(batch_ID_offset+seed_infection_scen_idx).txt")

        total_infected_holdings_results_array[:,seed_infection_scen_idx,behav_config_idx] = readdlm(current_disease_state_filename)[:,6] # Sixth column corresponds to cumulative infected holdings
        outbreak_duration_results_array[:,seed_infection_scen_idx,behav_config_idx] = readdlm(current_outbreak_duration_filename)
    end
end

# Get percentage results from absolute counts
n_holdings = 59774
percentage_infected_holdings_results_array = 100. .*(total_infected_holdings_results_array./n_holdings)

#========================
VIOLIN PLOTS
========================#

# Load intervention unit threshold cost data
intervention_unit_threshold_cost_array = load("JLD2_files/cost_related_statistics.jld2", "intervention_unit_threshold_cost_array")

# Replace NaN values with zeros
intervention_unit_threshold_cost_array[:,:,1] = zeros(Float64,n_simn_replicates,n_seed_region_scenarios)

# To outbreak size and duration result arrays, add a zero slice to correspond
# to all precautionary behaviour configuration (no outbreak scenario)
total_infected_holdings_results_incl_no_outbreak = cat(total_infected_holdings_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)
percentage_infected_holdings_results_incl_no_outbreak = cat(percentage_infected_holdings_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)
outbreak_duration_results_temp_incl_no_outbreak = cat(outbreak_duration_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)

# Assign data independent of seed location to 2d array
total_infected_holdings_results_all_seed_locs = reshape(total_infected_holdings_results_incl_no_outbreak,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs+1))
percentage_infected_holdings_results_all_seed_locs = reshape(percentage_infected_holdings_results_incl_no_outbreak,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs+1))
outbreak_duration_results_all_seed_locs = reshape(outbreak_duration_results_temp_incl_no_outbreak,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs+1))
intervention_unit_threshold_cost_array_all_seed_locs = reshape(intervention_unit_threshold_cost_array,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs+1))

# Configuration order for whisker plots
# Configuration 1: Infected premises control only.
# Configuration 2: All holdings apply both interventions upon confirmed infection within 50km (30 miles).
# Configuration 3: All holdings apply both interventions upon confirmed infection within 320km (200 miles).
# Configuration 4: All holdings apply both interventions prior to pathogen emergence (no outbreak occurs).
# Configuration 5: Uniform split of holdings across four intervention stance groups
# Configuration 6: Uniform split of holdings across three intervention timing groups (none assigned to non-user group)
# Configuration 7: Holdings partitioned into four behavioural groups using empirical data (model with two most stable variables).
# Configuration 8: Holdings partitioned into three behavioural groups using empirical data, includes dependence on herd size (model with five most stable variables).

# Mapping from result arrays slices to above configuration order
map_result_array_slice_to_whisker_plot_order = [1;2;3;8;4;5;6;7]

# Set up tuple of arrays with data inputs
# Use slice mapping to reorder columns for the whisker plots
data_tuple = [total_infected_holdings_results_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order],
              percentage_infected_holdings_results_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order],
              outbreak_duration_results_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order],
              intervention_unit_threshold_cost_array_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order]]


# Set up x-axis plot values
x_vals = [1;3:1:5;7:1:8;10:1:11]

# Set tick label for each configuration
config_xtick_label = ["Uncooperative";
                    "Strong parasitism";
                    "Weak parasitism";
                    "Mutual cooperation";
                    "Coop-Parasitism-FR";
                    "Coop-Parasitism";
                    "Trust-Expectancy";
                    "Herd size dependent"]

# Set up yaxis labels
y_axis_label_vec = ["Number of holdings infected";
                    "Holdings infected (%)";
                    "Outbreak duration (days)";
                    "Vaccine dose threshold cost\n"]

# Set up yaxis limit values
y_lim_max_val = [6.1e4;
                 103.;
                 500.;
                 12.]

y_lim_vals = [ (0. , 6.1e4),
               (0. , 103.),
               (0. , 500.),
               (0. , 12.) ]

# Set up save filenames
violin_plot_save_filename_vec = ["violin_plots/500_replicate_runs/total_infected_violin_plot";
                                "violin_plots/500_replicate_runs/percentage_infected_violin_plot";
                                "violin_plots/500_replicate_runs/outbreak_duration_violin_plot";
                                "violin_plots/500_replicate_runs/threshold_intervention_cost_violin_plot"]

# Specify the quantile values to be computed for this analysis
uncertainty_interval_quantile_vals = [0.025, 0.975]

# Generate the plots
for plot_idx = 1:length(data_tuple)
    plot_violin_figs_fn(data_tuple[plot_idx],
                                 uncertainty_interval_quantile_vals,
                                 n_behav_group_configs,
                                 x_vals,
                                 config_xtick_label,
                                 y_axis_label_vec[plot_idx],
                                 y_lim_max_val[plot_idx],
                                 plot_fontsize_violin_plot,
                                 save_fig_flag,
                                 violin_plot_save_filename_vec[plot_idx],
                                 fig_width_violin_plot,
                                 fig_height_violin_plot)
end


#========================
OUTBREAK THRESHOLD STATISTICS
GROUPED BAR PLOTS: Aggregating data over all seed infection locations
(Figure ZZZ)
========================#

### Set shared plotting variables ###

# Mapping from result arrays slices to this configuration order
# Set tick label for each configuration
map_result_array_slice_to_plot_order = [1;2;3;8;4;5;6;7]

config_xtick_label = ["Uncooperative";
                    "Strong parasitism";
                    "Weak parasitism";
                    "Mutual cooperation";
                    "Coop-Parasitism-FR";
                    "Coop-Parasitism";
                    "Trust-Expectancy";
                    "Herd size dependent"]

# Set number of simulations performed per behavioural group 
n_simns_per_config = n_simn_replicates*n_seed_region_scenarios

# Set up x-axis plot values
x_tick_vals = [1;3:1:5;7:1:8;10:1:11]

### Outbreak size (percentage) ###

# Set threshold values to assess
outbreak_threshold_size_vec = [1;10;20]

# Map from result arrays slices to the specified configuration order
reordered_percentage_infected_holdings_results_all_seed_locs = total_infected_holdings_results_all_seed_locs[:,map_result_array_slice_to_plot_order].*(100/n_holdings)

# Set y-axis label
yaxis_label_threshold_outbreak_size = "% holdings infected > X (% of simns) "

# Legend input values
legend_position_vals_threshold_outbreak_size = (0.86,0.92)
legend_labels_threshold_outbreak_size = ["X = 1%" "X = 10%" "X = 20%"]

# Set up plot bar colours
bar_colour_scheme_size = palette(:Blues_3,3)
bar_fillcolor_threshold_outbreak_size = [bar_colour_scheme_size[1] bar_colour_scheme_size[2] bar_colour_scheme_size[3]]

# Specify the save filename
threshold_outbreak_size_bar_plot_save_filename = "bar_plots/500_replicate_runs/threshold_outbreak_percentage_bar_plot"

# Generate the grouped bar plot
generate_threshold_grouped_bar_fig(reordered_percentage_infected_holdings_results_all_seed_locs,
                                    outbreak_threshold_size_vec,
                                    n_simns_per_config,
                                    x_tick_vals,
                                    config_xtick_label,
                                    yaxis_label_threshold_outbreak_size,
                                    legend_position_vals_threshold_outbreak_size,
                                    legend_labels_threshold_outbreak_size,
                                    bar_fillcolor_threshold_outbreak_size,
                                    save_fig_flag,
                                    threshold_outbreak_size_bar_plot_save_filename)


### Outbreak duration ###

# Set threshold values to assess
outbreak_threshold_durations_vec = [30;100;180]

# Map from result arrays slices to the specified configuration order
reordered_outbreak_duration_results_all_seed_locs = outbreak_duration_results_all_seed_locs[:,map_result_array_slice_to_plot_order]

# Set y-axis label
yaxis_label_threshold_outbreak_duration = "Outbreak duration > t (% of simns) "

# Legend input values
legend_position_vals_threshold_outbreak_duration = (0.82,0.92)
legend_labels_threshold_outbreak_duration = ["t = 30 days" "t = 100 days" "t = 180 days"]

# Set up plot bar colours
bar_colour_scheme_duration = palette(:Reds_3,3)
bar_fillcolor_threshold_outbreak_duration = [bar_colour_scheme_duration[1] bar_colour_scheme_duration[2] bar_colour_scheme_duration[3]]

# Specify the save filename
threshold_outbreak_duration_bar_plot_save_filename = "bar_plots/500_replicate_runs/threshold_outbreak_duration_bar_plot"

# Generate the grouped bar plot
generate_threshold_grouped_bar_fig(reordered_outbreak_duration_results_all_seed_locs,
                                    outbreak_threshold_durations_vec,
                                    n_simns_per_config,
                                    x_tick_vals,
                                    config_xtick_label,
                                    yaxis_label_threshold_outbreak_duration,
                                    legend_position_vals_threshold_outbreak_duration,
                                    legend_labels_threshold_outbreak_duration,
                                    bar_fillcolor_threshold_outbreak_duration,
                                    save_fig_flag,
                                    threshold_outbreak_duration_bar_plot_save_filename)

### Intervention unit cost ###

# Set threshold values to assess
intervention_unit_cost_threshold_vec = [0.5;1;2]

# Map from result arrays slices to the specified configuration order
reordered_intervention_unit_threshold_cost_results_all_seed_locs = intervention_unit_threshold_cost_array_all_seed_locs[:,map_result_array_slice_to_plot_order]

# Set y-axis label
yaxis_label_threshold_intervention_cost = "Threshold vaccine cost > C (% of simns) "

# Legend input values
legend_position_vals_threshold_intervention_cost= (0.10,0.92)
legend_labels_threshold_intervention_cost = ["C = 0.5" "C = 1.0" "C = 2.0"]

# Set up plot bar colours
bar_colour_scheme_intervention_cost = palette(:Purples_3,3)
bar_fillcolor_threshold_intervention_cost = [bar_colour_scheme_intervention_cost[1] bar_colour_scheme_intervention_cost[2] bar_colour_scheme_intervention_cost[3]]

# Specify the save filename
threshold_intervention_cost_bar_plot_save_filename = "bar_plots/500_replicate_runs/threshold_intervention_cost_bar_plot"

# Generate the grouped bar plot
generate_threshold_grouped_bar_fig(reordered_intervention_unit_threshold_cost_results_all_seed_locs,
                                    intervention_unit_cost_threshold_vec,
                                    n_simns_per_config,
                                    x_tick_vals,
                                    config_xtick_label,
                                    yaxis_label_threshold_intervention_cost,
                                    legend_position_vals_threshold_intervention_cost,
                                    legend_labels_threshold_intervention_cost,
                                    bar_fillcolor_threshold_intervention_cost,
                                    save_fig_flag,
                                    threshold_intervention_cost_bar_plot_save_filename)



#========================
SPATIAL PLOTS PRE-PROCESSING 
Set up linked regions & quantile parameters
========================#

# Import linked region array
# Row for valid seed regions per scenario
seed_region_scenarios = readdlm("../../../src/livestock_disease_model/GB_model_with_elicitation_data/region_seed_infection_scenario_files/county_ID_groups_for_region_seed_infection.txt") # Load the parameter combinations from file

# Assign number of scenarios to variable
n_scenarios = size(seed_region_scenarios,1)

# Initialise vector to store scenario ID for each region
region_seed_infection_scenario_ID = zeros(Int64,n_counties_UAs)

# Iterate over each scenario.
# Assign that scenario ID to the regions that were
# seed infection locations (in that given scenario)
for scenario_ID = 1:n_scenarios
    # Extract regions that were valid seed infection locations.
    # Convert to a vector
    selected_seed_region_data = seed_region_scenarios[scenario_ID,:]
    selected_seed_region_IDs = convert(Vector{Int64},selected_seed_region_data[selected_seed_region_data .!= ""])

    # For the seed infection regions in this scenarios, add scenario_ID to vector
    # storing scenario IDs where each region was a valid seed infection location
    region_seed_infection_scenario_ID[selected_seed_region_IDs] .= scenario_ID
end

# Specify the quantile values to be computed
quantile_vals_spatial_maps = [0.025, 0.5, 0.975]
n_quantile_vals = length(quantile_vals_spatial_maps)

# Set up plot title vector
plot_quantile_title_vec = Vector{String}(undef,3)
for quantile_idx = 1:n_quantile_vals
    # Set up titles, with median singled out
    if quantile_vals_spatial_maps[quantile_idx] == 0.5
        plot_quantile_title_vec[quantile_idx] = "Median"
    else
        plot_quantile_title_vec[quantile_idx] = "$(quantile_vals_spatial_maps[quantile_idx]*100)th percentile"
    end
end

#========================
SPATIAL MAPS BY SEED INFECTION LOCATION
INFECTION & OUTBREAK DURATION
(FIGURES YYY)
========================#

# To outbreak size and duration result arrays, add a zero slice to correspond
# to all precautionary behaviour configuration (no outbreak scenario)
total_infected_holdings_results_incl_no_outbreak = cat(total_infected_holdings_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)
percentage_infected_holdings_results_incl_no_outbreak = cat(percentage_infected_holdings_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)
outbreak_duration_results_temp_incl_no_outbreak = cat(outbreak_duration_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)

# Reset number of behavioural configurations values
n_behav_group_configs = size(percentage_infected_holdings_results_incl_no_outbreak,3)

# Set the figure panel titles
plot_config_title_vec = ["Uncooperative";
                        "Strong parasitism";
                        "Weak parasitism";
                        "Mutual cooperation";
                        "Coop-Parasitism-FR";
                        "Coop-Parasitism";
                        "Trust-Expectancy";
                        "Herd size dependent"
                        ]

# Set mapping from initial results array to match planned plot order for beahvioural configurations
map_result_array_slice_to_plot_order = [1;2;3;8;4;5;6;7]

 # Percentage of holdings infected (for given seed infection region)
 percentage_infected_holdings_cmap = :Blues
 percentage_infected_holdings_clim_vals = (0.,100.)
 percentage_infected_holdings_colorbar_string = "  \nHoldings infected (%)"
 percentage_infected_holdings_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/percentage_outbreak_size/holding_percentage_outbreak_size_seed_region_spatial_map_behav"
 plot_spatial_map_fn(percentage_infected_holdings_results_incl_no_outbreak[:,:,map_result_array_slice_to_plot_order],
                   quantile_vals_spatial_maps,
                   n_quantile_vals,
                   n_seed_region_scenarios,
                   n_behav_group_configs,
                   n_counties_UAs,
                   region_seed_infection_scenario_ID,
                   plot_config_title_vec,
                   percentage_infected_holdings_cmap,
                   percentage_infected_holdings_clim_vals,
                   percentage_infected_holdings_colorbar_string,
                   plot_fontsize_default,
                   save_fig_flag,
                   percentage_infected_holdings_save_filename_prefix,
                   fig_width_spatial_map,
                   fig_height_spatial_map)

 # Outbreak duration (for given seed infection region)
 outbreak_duration_cmap = :Reds
 outbreak_duration_clim_vals = (0.,500.)
 outbreak_duration_colorbar_string = "  \nOutbreak duration (days)"
 outbreak_duration_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/outbreak_duration/outbreak_duration_seed_region_spatial_map_behav"
 plot_spatial_map_fn(outbreak_duration_results_temp_incl_no_outbreak[:,:,map_result_array_slice_to_plot_order],
                quantile_vals_spatial_maps,
                n_quantile_vals,
                n_seed_region_scenarios,
                n_behav_group_configs,
                n_counties_UAs,
                region_seed_infection_scenario_ID,
                plot_config_title_vec,
                outbreak_duration_cmap,
                outbreak_duration_clim_vals,
                outbreak_duration_colorbar_string,
                plot_fontsize_default,
                save_fig_flag,
                outbreak_duration_save_filename_prefix,
                fig_width_spatial_map,
                fig_height_spatial_map)

#========================
SPATIAL MAPS BY SEED INFECTION LOCATION
THRESHOLD COST PER UNIT INTERVENTION
(Plot title: Configuration)
========================#

# Load intervention unit threshold cost data
intervention_unit_threshold_cost_array = load("JLD2_files/cost_related_statistics.jld2", "intervention_unit_threshold_cost_array")

# Replace NaN values with zeros
intervention_unit_threshold_cost_array[:,:,1] = zeros(Float64,n_simn_replicates,n_seed_region_scenarios)

# Specify plot related values
n_behav_group_configs_intervention_unit_threshold_cost = 8 # Config 8 for all holdings using both interventions & zero infections
intervention_unit_threshold_cost_cmap = :Purples
intervention_unit_threshold_cost_clim_vals = (0.,2.)
intervention_unit_threshold_cost_colorbar_string = "  \nVaccine dose threshold cost"

# Call plot function (behavioural configuration as titles)
intervention_unit_threshold_cost_behav_config_plot_title_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/intervention_unit_threshold_cost/intervention_unit_threshold_cost_seed_region_spatial_map_behav"
plot_spatial_map_fn(intervention_unit_threshold_cost_array[:,:,map_result_array_slice_to_plot_order],
               quantile_vals_spatial_maps,
               n_quantile_vals,
               n_seed_region_scenarios,
               n_behav_group_configs_intervention_unit_threshold_cost,
               n_counties_UAs,
               region_seed_infection_scenario_ID,
               plot_config_title_vec,
               intervention_unit_threshold_cost_cmap,
               intervention_unit_threshold_cost_clim_vals,
               intervention_unit_threshold_cost_colorbar_string,
               plot_fontsize_default,
               save_fig_flag,
               intervention_unit_threshold_cost_behav_config_plot_title_save_filename_prefix,
               fig_width_spatial_map,
               fig_height_spatial_map)


#========================
OUTBREAK THRESHOLD STATISTICS
SPATIAL MAPS
(FIGURES YYY)
========================#

# Mapping from result arrays slices to this configuration order
map_result_array_slice_to_plot_order = [1;2;3;8;4;5;6;7]

# Set up the quantiles to be evaluated
quantile_vals = [0.025;0.975]
n_quantile_vals = length(quantile_vals)

# Set number of configurations in use
n_behav_group_configs_incl_no_outbreak = 8

# Assign the number of configurations in use
n_configs = size(total_infected_holdings_results_incl_no_outbreak,3)

# Set titles for panels in the figures
plot_config_title_vec = ["Uncooperative";
                        "Strong parasitism";
                        "Weak parasitism";
                        "Mutual cooperation";
                        "Coop-Parasitism-FR";
                        "Coop-Parasitism";
                        "Trust-Expectancy";
                        "Herd size dependent"]


### Outbreak size (percentage) ###

# Set outbreak threshold sizes
outbreak_threshold_size = [1;10;20]
n_outbreak_threshold_sizes = length(outbreak_threshold_size)

# Mapping from result arrays slices to above configuration order
reordered_total_infected_holdings_results_incl_no_outbreak = total_infected_holdings_results_incl_no_outbreak[:,:,map_result_array_slice_to_plot_order].*(100/n_holdings)

propn_outbreaks_size_above_threshold,
quantiles_propn_outbreaks_size_above_threshold = get_propn_simn_satisfying_threshold_summ_stats(reordered_total_infected_holdings_results_incl_no_outbreak,
                                                                                                        outbreak_threshold_size,
                                                                                                        n_outbreak_threshold_sizes,
                                                                                                        n_seed_region_scenarios,
                                                                                                        n_configs,
                                                                                                        n_simn_replicates,
                                                                                                        quantile_vals,
                                                                                                        n_quantile_vals)

# Set colourmap and limits
outbreaks_within_threshold_size_cmap = :Blues
outbreaks_within_threshold_size_clim_vals = (0.,100.)

# Set colourbar label
colourbar_string_prefix_outbreak_threshold_size = "\nInfected holdings"
colourbar_string_suffix_outbreak_threshold_size = "%"

# Generate the spatial maps
outbreaks_within_threshold_size_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/threshold_outbreak_size/threshold_outbreak_size_empirical_estimates/threshold_outbreak_size_seed_region_spatial_map_behav"
plot_outbreak_threshold_spatial_map_fn(propn_outbreaks_size_above_threshold,
                                           outbreak_threshold_size,
                                           n_outbreak_threshold_sizes,
                                           n_seed_region_scenarios,
                                           n_behav_group_configs_incl_no_outbreak,
                                           n_counties_UAs,
                                           region_seed_infection_scenario_ID,
                                           plot_config_title_vec,
                                           outbreaks_within_threshold_size_cmap,
                                           outbreaks_within_threshold_size_clim_vals,
                                           colourbar_string_prefix_outbreak_threshold_size,
                                           colourbar_string_suffix_outbreak_threshold_size,
                                           plot_fontsize_default,
                                           save_fig_flag,
                                           outbreaks_within_threshold_size_save_filename_prefix,
                                           fig_width_spatial_map,
                                           fig_height_spatial_map)

### Outbreak duration ###

# Set outbreak threshold sizes
outbreak_threshold_duration = [30;100;180]
n_outbreak_threshold_durations = length(outbreak_threshold_duration)

# Mapping from result arrays slices to above configuration order
reordered_outbreak_duration_results_temp_incl_no_outbreak = outbreak_duration_results_temp_incl_no_outbreak[:,:,map_result_array_slice_to_plot_order]

propn_outbreaks_above_duration_threshold,
quantiles_propn_outbreaks_above_duration_threshold = get_propn_simn_satisfying_threshold_summ_stats(reordered_outbreak_duration_results_temp_incl_no_outbreak,
                                                                                                        outbreak_threshold_duration,
                                                                                                        n_outbreak_threshold_durations,
                                                                                                        n_seed_region_scenarios,
                                                                                                        n_configs,
                                                                                                        n_simn_replicates,
                                                                                                        quantile_vals,
                                                                                                        n_quantile_vals)

# Set colourmap and limits
outbreaks_within_threshold_duration_cmap = :Reds
outbreaks_within_threshold_duration_clim_vals = (0.,100.)

# Set colourbar label
colourbar_string_prefix_outbreak_threshold_duration = "\nOutbreak duration"
colourbar_string_suffix_outbreak_threshold_duration = " days"

# Generate the spatial maps
outbreaks_within_threshold_duration_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/threshold_outbreak_duration/threshold_outbreak_duration_empirical_estimates/threshold_outbreak_duration_seed_region_spatial_map_behav"
plot_outbreak_threshold_spatial_map_fn(propn_outbreaks_above_duration_threshold,
                                           outbreak_threshold_duration,
                                           n_outbreak_threshold_durations,
                                           n_seed_region_scenarios,
                                           n_behav_group_configs_incl_no_outbreak,
                                           n_counties_UAs,
                                           region_seed_infection_scenario_ID,
                                           plot_config_title_vec,
                                           outbreaks_within_threshold_duration_cmap,
                                           outbreaks_within_threshold_duration_clim_vals,
                                           colourbar_string_prefix_outbreak_threshold_duration,
                                           colourbar_string_suffix_outbreak_threshold_duration,
                                           plot_fontsize_default,
                                           save_fig_flag,
                                           outbreaks_within_threshold_duration_save_filename_prefix,
                                           fig_width_spatial_map,
                                           fig_height_spatial_map)

### Maximum intervention unit cost (whilst maintaining cost effectiveness) ###

# Set outbreak threshold sizes
intervention_unit_cost_threshold_vec = [0.5;1;2]
n_outbreak_threshold_intervention_cost = length(intervention_unit_cost_threshold_vec)

# Mapping from result arrays slices to above configuration order
reordered_intervention_unit_threshold_cost_array = intervention_unit_threshold_cost_array[:,:,map_result_array_slice_to_plot_order]

propn_outbreaks_above_intervention_cost_threshold,
quantiles_propn_outbreaks_above_intervention_cost_threshold = get_propn_simn_satisfying_threshold_summ_stats(reordered_intervention_unit_threshold_cost_array,
                                                                                                            intervention_unit_cost_threshold_vec,
                                                                                                            n_outbreak_threshold_intervention_cost,
                                                                                                            n_seed_region_scenarios,
                                                                                                            n_configs,
                                                                                                            n_simn_replicates,
                                                                                                            quantile_vals,
                                                                                                            n_quantile_vals)

# Set colourmap and limits
outbreaks_within_threshold_intervention_cost_cmap = :Purples
outbreaks_within_threshold_intervention_cost_clim_vals = (0.,100.)

# Set colourbar label
colourbar_string_prefix_outbreak_threshold_intervention_cost = "\nThreshold vaccine dose cost"
colourbar_string_suffix_outbreak_threshold_intervention_cost = ""

# Generate the spatial maps
outbreaks_within_threshold_intervention_cost_save_filename_prefix = "seed_region_spatial_maps/500_replicate_runs/plot_title_behav_config/threshold_intervention_cost/threshold_intervention_cost_empirical_estimates/threshold_intervention_cost_seed_region_spatial_map_behav"
plot_outbreak_threshold_spatial_map_fn(propn_outbreaks_above_intervention_cost_threshold,
                                           intervention_unit_cost_threshold_vec,
                                           n_outbreak_threshold_intervention_cost,
                                           n_seed_region_scenarios,
                                           n_behav_group_configs_incl_no_outbreak,
                                           n_counties_UAs,
                                           region_seed_infection_scenario_ID,
                                           plot_config_title_vec,
                                           outbreaks_within_threshold_intervention_cost_cmap,
                                           outbreaks_within_threshold_intervention_cost_clim_vals,
                                           colourbar_string_prefix_outbreak_threshold_intervention_cost,
                                           colourbar_string_suffix_outbreak_threshold_intervention_cost,
                                           plot_fontsize_default,
                                           save_fig_flag,
                                           outbreaks_within_threshold_intervention_cost_save_filename_prefix,
                                           fig_width_spatial_map,
                                           fig_height_spatial_map)