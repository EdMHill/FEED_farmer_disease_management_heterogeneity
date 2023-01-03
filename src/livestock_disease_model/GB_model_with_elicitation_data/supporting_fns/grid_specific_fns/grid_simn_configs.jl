#=
Purpose:
File to house variable configurations to input to Great Britain spatial
livestock disease transmission model simulation that uses the conditional subsample
gridding algorithm

Configuration list:
- GB_grid_config (Great Britain, FMD-like disease with non-IP interventions)
- GB_grid_config_IP_control_only (Great Britain, FMD-like disease with IP intervention only)
=#

#-------------------------------------------------------------------------------
### Great Britain, FMD-like disease with non-IP interventions
#-------------------------------------------------------------------------------
"""
    GB_grid_config(rng::AbstractRNG)

Setup configuration to be run using grid simulation method: Great Britain, FMD-like disease.
Interventions are applied to non-IP holdings.

Inputs:
- `rng::AbstractRNG`: The random number generator.

Outputs:
- `holding_info_vec::Vector{HoldingInfo}`: Structure with fields associated with holding specific data.
- `landscape_grid_data::GridData`: Spatial landscape information.
- `time_params::TimeParams`: Timestep per iteration & timeframe simulation is run over.
- `epi_param_vals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
- `control_params::ControlParams`: Variables related to implementing control measures
- `function_params::FunctionCollection`: Composite type containing collection of functions used in model simulation.
- `save_grid_config_filenames::Array{String,1}`: Names for three files. (i) Save array defining boundary limits of each cell within obtained grid configuration
                                  Array columns correspond to [xMin, xMax, yMin, yMax];
                              (ii) Grid ID each premises resides in;
                              (iii) For each grid ID, total number of premises in the cell.
- `replicate_file_prefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `output_file_objs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: grid\\_simn\\_configs.jl
"""
function GB_grid_config(rng::AbstractRNG)

    #Set bounding box limits for entire landscape
    bounding_box_left = 0.
    bounding_box_bottom = 0.
    bounding_box_right = 250000.
    bounding_box_top = 1000000.

    bounding_box_var = [bounding_box_left,bounding_box_right,bounding_box_bottom,bounding_box_top]
    #---------------------------------------------------------------------------

    ### Get holding location and herd size information ###
    holding_data_import = readdlm("../../../Data/dummy_data/GB_cattle_synthetic_data.txt",'\t',header=true)
    holding_data_array = holding_data_import[1]
    holding_data_header_vals = holding_data_import[2]
    df_holding_data = DataFrame(holding_data_array,vec(holding_data_header_vals))
        # Put into dataframe

    # Get the number of holdings
    n_holdings = size(df_holding_data,1)
    println("n_holdings: $n_holdings")

    # Assign location and livestock count data
    holding_loc_X_vals = convert(Vector{Float64},df_holding_data[!,:Easting])
    holding_loc_Y_vals = convert(Vector{Float64},df_holding_data[!,:Northing])
    holding_loc_X_and_Y_vals = [holding_loc_X_vals holding_loc_Y_vals]

    holding_livestock_data = convert(Vector{Int64},df_holding_data[!,:CattlePopulation])

    # Initialise HoldingInfo vector
    # Entry for each holding, with data fields associated with holding specific information
    holding_info_vec = Array{HoldingInfo,1}(undef,n_holdings)
    for holding_idx = 1:n_holdings
        holding_info_vec[holding_idx] = HoldingInfo(holding_ID = holding_idx,
                                                     country_ID = "NA",
                                                     county_ID = 0,
                                                     county_name = "NA",
                                                     X_loc = holding_loc_X_vals[holding_idx],
                                                     Y_loc = holding_loc_Y_vals[holding_idx],
                                                     livestock_count = [holding_livestock_data[holding_idx]])
    end
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    ### ERROR CHECKS
    #---------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING holding LOCATION DATA & LIVESTOCK POPN DATA!
    if size(holding_loc_X_and_Y_vals,1) != size(holding_livestock_data,1)
        error("Inconsistency in number of records in location dataset($(size(holding_loc_dataCoords,1))) and number of records in livestock dataset, $(size(holding_livestock_data,1)).")
    end

    check_landscape_valid(bounding_box_var,
                                holding_loc_X_and_Y_vals[:,1],
                                holding_loc_X_and_Y_vals[:,2])
    #---------------------------------------------------------------------------

    ### Construct holding location data composite type
    landscape_grid_data = GridData(coord_type = 1,
                                    bounding_box_vals = bounding_box_var,
                                    grid_optim_method = 2,
                                    nMax = 0,
                                    n_cells_per_side = 0)
    #---------------------------------------------------------------------------

    ### Set up time_params structure - Entries: timestep_val, max_time. ###
    time_params = TimeParams(timestep_val = 1.,
                                    max_time = 10*365.)
    #---------------------------------------------------------------------------

    ### Specify epidemiological parameters ###
    epi_params = EpiParams(incubation_time = 5.,
                                    detection_time = 9.,
                                    removal_time = 13.,
                                    per_livestock_type_suscep = [1.],
                                    per_livestock_type_transmiss = [1.0e6],
                                    suscept_exponent = [0.41],
                                    transmiss_exponent = [0.42])
    #---------------------------------------------------------------------------

    ### Specify control related functions ###
    control_params = ControlParams(run_controls_fn = control_housing_and_vacc_fn!,
                                            vacc_efficacy = 1.,
                                            holding_time_to_inoculation_fn = five_day_inoculation_fn!,
                                            house_livestock_effectiveness = 0.5)
    #---------------------------------------------------------------------------

    ### Specify functions to be used ###

    function_params = FunctionCollection(calc_holding_suscep_transmiss_fn = compute_holding_suscep_transmiss_cattle_only,
                                                kernel_fn = construct_power_law_transmission_kernel,
                                                iterate_outbreak_fn = iterate_outbreak_grid_simn!,
                                                run_simn_replicate_fn = run_outbreak_grid_simn)
    #---------------------------------------------------------------------------

    ### Savefile name and objects ###

    # save_grid_config_filenames - Names for three files:
    #                              (i) Save array defining boundary limits of each cell within obtained grid configuration
    #                                   Array columns correspond to [xMin, xMax, yMin, yMax]
    #                               (ii) Grid ID each premises resides in
    #                               (iii) For each grid ID, total number of premises in the cell
    save_grid_config_filenames = ["../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_cell_boundaries_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_gridIDs_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_per_cell_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_sort_idx_batchID$(batchID).txt"]

    # replicate_file_prefix - Directory location to be used with output files storing data with individual file per replicate
    replicate_file_prefix = "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs"

    #  output_file_objs - Filename identifiers. Used for files written to by all replicates.
    outbreak_duration_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/outbreak_duration_batchID$(batchID).txt", "a")
    holding_per_disease_state_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/holdings_per_disease_state_batchID$(batchID).txt", "a")
    cumulative_culled_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_culled_batchID$(batchID).txt", "a")
    cumulative_vacc_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_vacc_batchID$(batchID).txt", "a")
    cumulative_livestock_housed_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_livestock_housed_batchID$(batchID).txt", "a")
    cumulative_cases_livestock_type_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_cases_livestock_batchID$(batchID).txt", "a")
    output_file_objs = [outbreak_duration_file,
                        holding_per_disease_state_file,
                        cumulative_culled_file,
                        cumulative_vacc_file,
                        cumulative_livestock_housed_file,
                        cumulative_cases_livestock_type_file]
    #-------------------------------------------------------------------------------

    # Specify what is returned by the function
    return holding_info_vec::Vector{HoldingInfo},
            landscape_grid_data::GridData,
            time_params::TimeParams,
            epi_params::EpiParams,
            control_params::ControlParams,
            function_params::FunctionCollection,
            save_grid_config_filenames::Array{String,1},
            replicate_file_prefix::String,
            output_file_objs::Array{IOStream,1}
end

#-------------------------------------------------------------------------------
### Great Britain, FMD-like disease with IP intervention only
#-------------------------------------------------------------------------------
"""
    GB_grid_config_IP_control_only(rng::AbstractRNG)

Setup configuration to be run using grid simulation method: Great Britain, FMD-like disease.
No additional intervention baseline scenario (i.e. no interventions applied to non-IP holdings).

Inputs:
- `rng::AbstractRNG`: The random number generator.

Outputs:
- `holding_info_vec::Vector{HoldingInfo}`: Structure with fields associated with holding specific data.
- `landscape_grid_data::GridData`: Spatial landscape information.
- `time_params::TimeParams`: Timestep per iteration & timeframe simulation is run over.
- `epi_param_vals::Array{Float64,1}`: [Incubation time, Detection time, Recovery time (natural recovery), Removal time (culled)]
- `control_params::ControlParams`: Variables related to implementing control measures
- `function_params::FunctionCollection`: Composite type containing collection of functions used in model simulation.
- `save_grid_config_filenames::Array{String,1}`: Names for three files. (i) Save array defining boundary limits of each cell within obtained grid configuration
                                  Array columns correspond to [xMin, xMax, yMin, yMax];
                              (ii) Grid ID each premises resides in;
                              (iii) For each grid ID, total number of premises in the cell.
- `replicate_file_prefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `output_file_objs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Location: grid\\_simn\\_configs.jl
"""
function GB_grid_config_IP_control_only(rng::AbstractRNG)

    ### Get GB landscape bounding box information ###
    path = "../../../Data/Shapefiles/Counties_and_Unitary_Authorities_(December_2021)_UK_BGC/CTYUA_DEC_2021_UK_BGC.shp"
    table = Shapefile.Table(path)
    geoms = Shapefile.shapes(table)

    #Aggregate information of bounding box for each region
    n_regions = length(geoms)
    region_left_vals = zeros(n_regions); region_right_vals = zeros(n_regions);
    region_bottom_vals = zeros(n_regions); region_top_vals = zeros(n_regions);
    for ii = 1:n_regions
        region_left_vals[ii] = geoms[ii].MBR.left
        region_right_vals[ii] = geoms[ii].MBR.right
        region_bottom_vals[ii] = geoms[ii].MBR.bottom
        region_top_vals[ii] = geoms[ii].MBR.top
    end

    #Get bounding box limits for entire landscape
    bounding_box_left = minimum(region_left_vals)
    bounding_box_bottom = minimum(region_bottom_vals)
    bounding_box_right = maximum(region_right_vals)
    bounding_box_top = maximum(region_top_vals)

    bounding_box_var = [bounding_box_left,bounding_box_right,bounding_box_bottom,bounding_box_top]
    #---------------------------------------------------------------------------

    ### Get holding location and herd size information ###
    holding_data_import = readdlm("../../../Data/dummy_data/GB_cattle_synthetic_data.txt",'\t',header=true)
    holding_data_array = holding_data_import[1]
    holding_data_header_vals = holding_data_import[2]
    df_holding_data = DataFrame(holding_data_array,vec(holding_data_header_vals))
        # Put into dataframe

    # Get the number of holdings
    n_holdings = size(df_holding_data,1)
    println("n_holdings: $n_holdings")

    # Assign location and livestock count data
    holding_country_IDs = df_holding_data[!,:Country]
    holding_county_IDs = df_holding_data[!,"County ID"]
    holding_county_names = df_holding_data[!,"County Name"]
    holding_loc_X_vals = convert(Vector{Float64},df_holding_data[!,:Easting])
    holding_loc_Y_vals = convert(Vector{Float64},df_holding_data[!,:Northing])
    holding_loc_X_and_Y_vals = [holding_loc_X_vals holding_loc_Y_vals]

    holding_livestock_data = convert(Vector{Int64},df_holding_data[!,:CattlePopulation])

    # Initialise HoldingInfo vector
    # Entry for each holding, with data fields associated with holding specific information
    holding_info_vec = Array{HoldingInfo,1}(undef,n_holdings)
    for holding_idx = 1:n_holdings
        holding_info_vec[holding_idx] = HoldingInfo(holding_ID = holding_idx,
                                                     country_ID = holding_country_IDs[holding_idx],
                                                     county_ID = holding_county_IDs[holding_idx],
                                                     county_name = holding_county_names[holding_idx],
                                                     X_loc = holding_loc_X_vals[holding_idx],
                                                     Y_loc = holding_loc_Y_vals[holding_idx],
                                                     livestock_count = [holding_livestock_data[holding_idx]])
    end
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    ### ERROR CHECKS
    #---------------------------------------------------------------------------
    #CHECK THAT ROW COUNTS MATCH WHEN COMPARING holding LOCATION DATA & LIVESTOCK POPN DATA!
    if size(holding_loc_X_and_Y_vals,1) != size(holding_livestock_data,1)
        error("Inconsistency in number of records in location dataset($(size(holding_loc_dataCoords,1))) and number of records in livestock dataset, $(size(holding_livestock_data,1)).")
    end

    check_landscape_valid(bounding_box_var,
                                holding_loc_X_and_Y_vals[:,1],
                                holding_loc_X_and_Y_vals[:,2])
    #---------------------------------------------------------------------------

    ### Construct holding location data composite type
    landscape_grid_data = GridData(coord_type = 1,
                                    bounding_box_vals = bounding_box_var,
                                    grid_optim_method = 2,
                                    nMax = 0,
                                    n_cells_per_side = 0)
    #---------------------------------------------------------------------------

    ### Set up time_params structure - Entries: timestep_val, max_time. ###
    time_params = TimeParams(timestep_val = 1.,
                                    max_time = 10*365.)
    #---------------------------------------------------------------------------

    ### Specify epidemiological parameters ###
    epi_params = EpiParams(incubation_time = 5.,
                                    detection_time = 9.,
                                    removal_time = 13.,
                                    per_livestock_type_suscep = [1.],
                                    per_livestock_type_transmiss = [1.0e6],
                                    suscept_exponent = [0.41],
                                    transmiss_exponent = [0.42])
    #---------------------------------------------------------------------------

    ### Specify control related functions ###
    control_params = ControlParams(run_controls_fn = remove_IPs_only_control_fn!,
                                            vacc_efficacy = 1.,
                                            holding_time_to_inoculation_fn = five_day_inoculation_fn!,
                                            house_livestock_effectiveness = 0.5)
    #---------------------------------------------------------------------------

    ### Specify functions to be used ###

    function_params = FunctionCollection(calc_holding_suscep_transmiss_fn = compute_holding_suscep_transmiss_cattle_only,
                                                kernel_fn = construct_power_law_transmission_kernel,
                                                iterate_outbreak_fn = iterate_outbreak_grid_simn!,
                                                run_simn_replicate_fn = run_outbreak_grid_simn)
    #---------------------------------------------------------------------------

    ### Savefile name and objects ###

    # save_grid_config_filenames - Names for three files:
    #                              (i) Save array defining boundary limits of each cell within obtained grid configuration
    #                                   Array columns correspond to [xMin, xMax, yMin, yMax]
    #                               (ii) Grid ID each premises resides in
    #                               (iii) For each grid ID, total number of premises in the cell
    save_grid_config_filenames = ["../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_cell_boundaries_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_gridIDs_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_per_cell_batchID$(batchID).txt",
                                "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_grid_configs/grid_config_holding_sort_idx_batchID$(batchID).txt"]

    # replicate_file_prefix - Directory location to be used with output files storing data with individual file per replicate
    replicate_file_prefix = "../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs"

    #  output_file_objs - Filename identifiers. Used for files written to by all replicates.
    outbreak_duration_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/outbreak_duration_batchID$(batchID).txt", "a")
    holding_per_disease_state_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/holdings_per_disease_state_batchID$(batchID).txt", "a")
    cumulative_culled_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_culled_batchID$(batchID).txt", "a")
    cumulative_vacc_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_vacc_batchID$(batchID).txt", "a")
    cumulative_livestock_housed_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_livestock_housed_batchID$(batchID).txt", "a")
    cumulative_cases_livestock_type_file = open("../../../results/livestock_disease_model/GB_model_with_behaviour_groups_grid_simn_outputs/GB_model_epi_outputs_aggregated/cumulative_cases_livestock_batchID$(batchID).txt", "a")
    output_file_objs = [outbreak_duration_file,
                        holding_per_disease_state_file,
                        cumulative_culled_file,
                        cumulative_vacc_file,
                        cumulative_livestock_housed_file,
                        cumulative_cases_livestock_type_file]
    #-------------------------------------------------------------------------------

    # Specify what is returned by the function
    return holding_info_vec::Vector{HoldingInfo},
            landscape_grid_data::GridData,
            time_params::TimeParams,
            epi_params::EpiParams,
            control_params::ControlParams,
            function_params::FunctionCollection,
            save_grid_config_filenames::Array{String,1},
            replicate_file_prefix::String,
            output_file_objs::Array{IOStream,1}
end
