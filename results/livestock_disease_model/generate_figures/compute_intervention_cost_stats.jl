#===============================================================================
Purpose:
# Compute intervention cost associated statistics

# Date: August 2022
===============================================================================#

#========================
LOAD REQUIRED ENVIRONMENT
========================#

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../")

#========================
LOAD REQUIRED PACKAGES
========================#
using DelimitedFiles
using StatsBase
using DataFrames
using Shapefile
using JLD2

#========================
ALL CONTROLLED SCENARIO DATA
========================#

# Load livestock data
holding_data_import = readdlm("TBD.txt",'\t',header=true)
holding_data_array = holding_data_import[1]
holding_data_header_vals = holding_data_import[2]
df_holding_data = DataFrame(holding_data_array,vec(holding_data_header_vals))
    # Put into dataframe

holding_livestock_data = convert(Vector{Int64},df_holding_data[!,:CattlePopulation])

# Get total number of livestock
n_livestock = sum(holding_livestock_data)

#========================
IMPORT REQUIRED SCENARIO RESULT DATASETS
========================#

# Initialise simulation & scenario values
n_seed_region_scenarios = 89
n_replicates = 500
n_behav_group_configs = 8

# Specify batch ID offset value for each configuration
batch_ID_offset_vec = [2000;2400;2300;3500;3400;3600;3700]

# Set up arrays to store imported data
infected_count_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)
vacc_count_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)
livestock_housed_count_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)

# Iterate over each batchID.
# Load relevant data files
# Assign each statistic to associated array
for behav_group_scenario_ID = 1:(n_behav_group_configs-1) # Final behavioural group scenario dealt with after this loop
    for seed_region_ID = 1:n_seed_region_scenarios

        # Compute batchID for this iteration
        batchID_val = batch_ID_offset_vec[behav_group_scenario_ID] + seed_region_ID
        # batchID_val = batchID_offset_value + (behav_group_scenario_ID-1)*100 + seed_region_ID

        # Infection count for all replicates
        infected_count_array[:,seed_region_ID,behav_group_scenario_ID] = readdlm("../GB_model_epi_outputs_aggregated/cumulative_culled_batchID$(batchID_val).txt")

        # Vaccination count for all replicates
        vacc_count_array[:,seed_region_ID,behav_group_scenario_ID] = readdlm("../GB_model_epi_outputs_aggregated/cumulative_vacc_batchID$(batchID_val).txt")

        # Livestock housed count for all replicates
        livestock_housed_count_array[:,seed_region_ID,behav_group_scenario_ID] = readdlm("../GB_model_epi_outputs_aggregated/cumulative_livestock_housed_batchID$(batchID_val).txt")
    end
end

# behav_group_scenario_ID == n_behav_group_configs
# Corresponds to scenario where all holdings enact both controls, meaning there is no outbreak
infected_count_array[:,:,n_behav_group_configs] = zeros(Int64,n_replicates,n_seed_region_scenarios)
vacc_count_array[:,:,n_behav_group_configs] = n_livestock.*ones(Int64,n_replicates,n_seed_region_scenarios)
livestock_housed_count_array[:,:,n_behav_group_configs] = n_livestock.*ones(Int64,n_replicates,n_seed_region_scenarios)

#========================
INITIALISE COMPARATOR VARIABLES & ARRAYS
========================#

# Initialise cost variables
relative_vaccine_cost = 0:0.01:1
relative_livestock_housing_cost = 0:0.01:1

n_relative_vaccine_cost = length(relative_vaccine_cost)
n_relative_livestock_housing_cost = length(relative_livestock_housing_cost)

# Initialise arrays to store difference in count of livestock infected, vaccinated
# and housed between each scenario and baseline comparator scenario (IP control only)
infected_count_comparator_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)
vacc_count_comparator_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)
housed_count_comparator_array = zeros(Int64,n_replicates,n_seed_region_scenarios,n_behav_group_configs)
    # Column per seed region.
    # Row per replicate.
    # 3rd dimension slice per behavioural group config.

# Initialise relative cost associated arrays
cost_difference_array = zeros(Float64,n_replicates,n_seed_region_scenarios,n_behav_group_configs,n_relative_vaccine_cost,n_relative_livestock_housing_cost)
intervention_unit_threshold_cost_array = zeros(Float64,n_replicates,n_seed_region_scenarios,n_behav_group_configs,n_relative_vaccine_cost,n_relative_livestock_housing_cost)
    # Row per replicate.
    # Column per seed region.
    # 3rd dimension slice per behavioural group config.
    # 4th dimension slice per relative vaccination cost
    # 5th dimension slice per relative housing livestock cost

#========================
POPULATE COMPARATOR ARRAYS
========================#
for behav_group_scenario_ID = 1:n_behav_group_configs

        # Infection count for all replicates
        infected_count_comparator_array[:,:,behav_group_scenario_ID] = infected_count_array[:,:,behav_group_scenario_ID] - infected_count_array[:,:,1]

        # Vaccination count for all replicates
        vacc_count_comparator_array[:,:,behav_group_scenario_ID] = vacc_count_array[:,:,behav_group_scenario_ID] - vacc_count_array[:,:,1]

        # Livestock housed count for all replicates
        housed_count_comparator_array[:,:,behav_group_scenario_ID] = livestock_housed_count_array[:,:,behav_group_scenario_ID] - livestock_housed_count_array[:,:,1]
end

#========================
POPULATE COST DIFFERENCE ARRAYS
========================#

# Iterate over each relative vaccine cost and relative housing livestock cost
# Assign value for each cow avoiding infection of 1.
for rel_livestock_housed_cost_idx = 1:n_relative_livestock_housing_cost
    for rel_vacc_cost_idx = 1:n_relative_vaccine_cost

        # Populate cost difference array
        cost_difference_array[:,:,:,rel_vacc_cost_idx,rel_livestock_housed_cost_idx] =
            1 .*infected_count_comparator_array .+  # Relative cost contribution from infected livestock (will be negative if infections have been averted)
            (relative_vaccine_cost[rel_vacc_cost_idx].*vacc_count_comparator_array) .+ # Relative cost contribution from vaccinted livestock
            (relative_livestock_housing_cost[rel_livestock_housed_cost_idx].*housed_count_comparator_array)  # Relative cost contribution from housed livestock
    end
end

#========================
POPULATE THRESHOLD INTERVENTION UNIT COST ARRAYS
========================#

# Construct total intervention units used array
total_intervention_units_used_array = vacc_count_array + livestock_housed_count_array

# Compute array of infetions averted relative to comparator strategy
infections_averted = -1*infected_count_comparator_array

# Compute intevention unit threshold cost price
intervention_unit_threshold_cost_array = infections_averted./total_intervention_units_used_array

#========================
POPULATE percentage_simns_cost_effective_array
========================#

# Find where cost_difference_array < 0
# Corresponds to an overall cost saving.
# Sum over rows (i.e. summing over replicates for a fixed seed infetion region, behavioural group scenario & intervention costs)
propn_simns_cost_effective_array = sum(cost_difference_array .< 0,dims=1)/n_replicates
percentage_simns_cost_effective_array = dropdims(propn_simns_cost_effective_array; dims=1)*100
    # Row per seed region.
    # Column per behavioural group config.
    # 3rd dimension slice per relative vaccination cost
    # 4th dimension slice per relative housing livestock cost

#========================
SAVE ARRAYS TO FILE
========================#
# Save config details to JLD 2 file
jldsave("JLD2_files/cost_related_statistics.jld2";
        intervention_unit_threshold_cost_array)