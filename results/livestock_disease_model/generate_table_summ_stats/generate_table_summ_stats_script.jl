#=
Script to generate summary statistics for the tables

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
using JLD2
using StatsBase

#========================
IMPORT RESULTS DATA
========================#

# Set up filename prefixes
behav_group_config_disease_state_filename = "../GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/holdings_per_disease_state_"
behav_group_config_outbreak_duration_filename = "../GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/outbreak_duration_"

# Set up array construction variables
n_simn_replicates = 500
n_seed_region_scenarios = 89
n_behav_group_configs_simulated = 7
n_behav_group_configs_total = 8

# Initialise arrays to consolidate results from different scenarios
total_infected_holdings_results_array = zeros(Int64,n_simn_replicates,n_seed_region_scenarios,n_behav_group_configs_simulated)
outbreak_duration_results_array = zeros(Int64,n_simn_replicates,n_seed_region_scenarios,n_behav_group_configs_simulated)

# Specify values for batch_ID_offset
batch_ID_offset_vec = [5000;5100;5200;5300;5400;5500;5600]

# Iterate over behavioural group configs & specified batch IDs.
for behav_config_idx = 1:n_behav_group_configs_simulated

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
REORDER DATA TO MATCH POSITIONING OF BEHAVIOURAL CONFIGS IN THE FIGURES
========================#

# Load intervention unit threshold cost data
intervention_unit_threshold_cost_array = load("../generate_figures/JLD2_files/cost_related_statistics.jld2", "intervention_unit_threshold_cost_array")

# Replace NaN values with zeros (for the baseline scenario)
intervention_unit_threshold_cost_array[:,:,1] = zeros(Float64,n_simn_replicates,n_seed_region_scenarios)

# To outbreak size and duration result arrays, add a zero slice to correspond
# to all precautionary behaviour configuration (no outbreak scenario)
percentage_infected_holdings_results_incl_no_outbreak = cat(percentage_infected_holdings_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)
outbreak_duration_results_temp_incl_no_outbreak = cat(outbreak_duration_results_array,zeros(Int64,n_simn_replicates,n_seed_region_scenarios);dims=3)

# Assign data independent of seed location to 2d array
percentage_infected_holdings_results_all_seed_locs = reshape(percentage_infected_holdings_results_incl_no_outbreak,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs_total))
outbreak_duration_results_all_seed_locs = reshape(outbreak_duration_results_temp_incl_no_outbreak,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs_total))
intervention_unit_threshold_cost_array_all_seed_locs = reshape(intervention_unit_threshold_cost_array,(n_simn_replicates*n_seed_region_scenarios,n_behav_group_configs_total))

# Reorder the configurations

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

# Use slice mapping to reorder columns into plot x-axis order
reordered_percentage_infected_holdings_results_all_seed_locs = percentage_infected_holdings_results_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order]
reordered_outbreak_duration_results_all_seed_locs = outbreak_duration_results_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order]
reordered_intervention_unit_threshold_cost_results_all_seed_locs = intervention_unit_threshold_cost_array_all_seed_locs[:,map_result_array_slice_to_whisker_plot_order]

#========================
VIOLIN PLOTS SUMMARY STATISTICS
========================#

# Get median values
medians_percentage_infected_holdings_results_all_seed_locs = median(reordered_percentage_infected_holdings_results_all_seed_locs,dims=1)
medians_outbreak_duration_results_all_seed_locs = median(reordered_outbreak_duration_results_all_seed_locs,dims=1)
medians_intervention_unit_threshold_cost_results_all_seed_locs = median(reordered_intervention_unit_threshold_cost_results_all_seed_locs,dims=1)

# Specify the quantile values to be computed for this analysis, initialise storage arrays & compute the requested quantiles
quantile_vals = [0.025; 0.975]
n_quantile_vals = length(quantile_vals)

quantiles_percentage_infected_holdings_results_all_seed_locs = zeros(Float64,n_quantile_vals,n_behav_group_configs_total)
quantiles_outbreak_duration_results_all_seed_locs = zeros(Float64,n_quantile_vals,n_behav_group_configs_total)
quantiles_intervention_unit_threshold_cost_results_all_seed_locs = zeros(Float64,n_quantile_vals,n_behav_group_configs_total)

for behav_config_idx = 1:n_behav_group_configs_total
     
    # Extract the data relevant to the behavioural config under consideration in this loop iteration
    current_behav_config_itr_percentage_infected_holdings_results_all_seed_locs = reordered_percentage_infected_holdings_results_all_seed_locs[:,behav_config_idx]
    current_behav_config_itr_outbreak_duration_results_all_seed_locs = reordered_outbreak_duration_results_all_seed_locs[:,behav_config_idx]
    current_behav_config_itr_intervention_unit_threshold_cost_results_all_seed_locs = reordered_intervention_unit_threshold_cost_results_all_seed_locs[:,behav_config_idx]

    # Compute the requested quantile values for these outputs
    for quantile_val_idx = 1:n_quantile_vals
        quantiles_percentage_infected_holdings_results_all_seed_locs[quantile_val_idx,behav_config_idx] = quantile(current_behav_config_itr_percentage_infected_holdings_results_all_seed_locs,quantile_vals[quantile_val_idx])
        quantiles_outbreak_duration_results_all_seed_locs[quantile_val_idx,behav_config_idx] = quantile(current_behav_config_itr_outbreak_duration_results_all_seed_locs,quantile_vals[quantile_val_idx])
        quantiles_intervention_unit_threshold_cost_results_all_seed_locs[quantile_val_idx,behav_config_idx] = quantile(current_behav_config_itr_intervention_unit_threshold_cost_results_all_seed_locs,quantile_vals[quantile_val_idx])
    end
end

#========================
OUTBREAK THRESHOLD STATISTICS
========================#

### Specify total number of simulation performed per behavioural configuration ###
total_simns_per_behav_config = n_simn_replicates*n_seed_region_scenarios

### Set threshold values to assess ###
outbreak_threshold_size_vec = [1;10;20]
outbreak_threshold_durations_vec = [30;100;180]
intervention_unit_cost_threshold_vec = [0.5;1;2]

### Initialise arrays to store output statistics ###
threshold_outbreak_size_results = zeros(Float64,length(outbreak_threshold_size_vec),n_behav_group_configs_total)
threshold_outbreak_duration_results = zeros(Float64,length(outbreak_threshold_durations_vec),n_behav_group_configs_total)
threshold_intervention_cost_results = zeros(Float64,length(intervention_unit_cost_threshold_vec),n_behav_group_configs_total)

### Populate the threshold criteria output arrays ###
for behav_config_idx = 1:n_behav_group_configs_total

    # Extract the data relevant to the behavioural config under consideration in this loop iteration
    current_behav_config_itr_percentage_infected_holdings_results_all_seed_locs = reordered_percentage_infected_holdings_results_all_seed_locs[:,behav_config_idx]
    current_behav_config_itr_outbreak_duration_results_all_seed_locs = reordered_outbreak_duration_results_all_seed_locs[:,behav_config_idx]
    current_behav_config_itr_intervention_unit_threshold_cost_results_all_seed_locs = reordered_intervention_unit_threshold_cost_results_all_seed_locs[:,behav_config_idx]
    
    # Outbreak size threshold check
    for size_threshold_val_itr = 1:length(outbreak_threshold_size_vec)
        threshold_outbreak_size_results[size_threshold_val_itr,behav_config_idx] = ( sum(current_behav_config_itr_percentage_infected_holdings_results_all_seed_locs .> outbreak_threshold_size_vec[size_threshold_val_itr])*100 ) / total_simns_per_behav_config
    end

    # Outbreak duration threshold check
    for duration_threshold_val_itr = 1:length(outbreak_threshold_size_vec)
        threshold_outbreak_duration_results[duration_threshold_val_itr,behav_config_idx] = ( sum(current_behav_config_itr_outbreak_duration_results_all_seed_locs .> outbreak_threshold_durations_vec[duration_threshold_val_itr])*100 ) / total_simns_per_behav_config
    end

    # Intervention unit cost threshold check
    for cost_threshold_val_itr = 1:length(outbreak_threshold_size_vec)
        threshold_intervention_cost_results[cost_threshold_val_itr,behav_config_idx] = ( sum(current_behav_config_itr_intervention_unit_threshold_cost_results_all_seed_locs .> intervention_unit_cost_threshold_vec[cost_threshold_val_itr])*100 ) / total_simns_per_behav_config
    end
end