#=
Purpose:
Run an individual-based livestock infection amongst Great Britain cattle holdings
using the conditional subsample grid-based algorithm

Julia version: 1.8
=#

#-------------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT
#-------------------------------------------------------------------------------

#Set path to directory this file resides in
cd(dirname(@__FILE__))

using Pkg
Pkg.activate("../../../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using DataFrames
using DelimitedFiles
using Distributions
using Statistics
using StatsBase
using Random
using JLD2
using Shapefile
using Parameters
using QuadGK

#-------------------------------------------------------------------------------
# IMPORT REQUIRED FUNCTION FILES
#-------------------------------------------------------------------------------
include("../../common_fns/CommonFunctions.jl")
include("supporting_fns/define_structs.jl")
include("supporting_fns/distance_fns.jl")
include("supporting_fns/initialise_behavioural_groups.jl")
include("supporting_fns/supporting_fns_transmission_dynamics.jl")
include("supporting_fns/supporting_fns_interventions.jl")
include("supporting_fns/grid_optimisation/adaptive_grid_alg_fns.jl")
include("supporting_fns/grid_optimisation/approx_node_per_cell_threshold_fns.jl")
include("supporting_fns/grid_specific_fns/grid_fns.jl")
include("supporting_fns/grid_specific_fns/grid_simn_configs.jl")
include("supporting_fns/grid_specific_fns/supporting_fns_run_model_grid_simns.jl")

#-------------------------------------------------------------------------------
# FUNCTIONS FOR RUNNING SPATIAL MODEL
#-------------------------------------------------------------------------------
"""
    spatial_livestock_grid_infection_model(RNGseed::Int64,
                                            batchID::String,
                                            n_reps::Int64,
                                            holding_info_vec::Array{HoldingInfo,1},
                                            landscape_grid_data::GridData,
                                            behaviour_group_config_fn::Function,
                                            time_params::TimeParams,
                                            epi_params::EpiParams,
                                            control_params::ControlParams,
                                            function_params::FunctionCollection,
                                            save_grid_config_filenames::Array{String,1},
                                            replicate_file_prefix::String,
                                            output_file_objs::Array{IOStream,1})

Wrapper to run outbreak model on selected configuration using the grid method.

Inputs:
- `RNGseed::Int64`: Input to seed RNG
- `batchID::String`: Used as prefix for output files.
- `n_reps::Int64`: Number of replicates to run.
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `landscape_grid_data::GridData`: Composite type with info such as co-ordinate system type & bounding box information [Min_x,Max_x,Min_y,Max_y]
- `time_params::TimeParams`: Timestep per iteration & timeframe simulation is run over.
- `epi_params::EpiParams`: Composite type containing epimidological parameters.
- `control_params::ControlParams`: Variables related to implementing control measures
- `behaviour_group_config_fn::Function`: Load desired behavioural group configuration from a premade function
- `function_params::FunctionCollection`: Composite type containing collection of functions used in model simulation.
- `save_grid_config_filenames::Array{String,1}`: Names for three files. (i) Save array defining boundary limits of each cell within obtained grid configuration
                                  Array columns correspond to [xMin, xMax, yMin, yMax];
                              (ii) Grid ID each holding resides in;
                              (iii) For each grid ID, total number of holding in the cell.
- `replicate_file_prefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `output_file_objs::Array{IOStream,1}`: Filename identifiers written to by all replicates.

Outputs: None \n
Location: run\\_generic\\_livestock\\_disease\\_spatial\\_model.jl
"""
function spatial_livestock_grid_infection_model(RNGseed::Int64,
                                            batchID::String,
                                            n_reps::Int64,
                                            holding_info_vec::Array{HoldingInfo,1},
                                            landscape_grid_data::GridData,
                                            behaviour_group_config_fn::Function,
                                            time_params::TimeParams,
                                            epi_params::EpiParams,
                                            control_params::ControlParams,
                                            function_params::FunctionCollection,
                                            save_grid_config_filenames::Array{String,1},
                                            replicate_file_prefix::String,
                                            output_file_objs::Array{IOStream,1})

    #---------------------------------------------------------------------------
    # CHECK CO-ORD SYSTEM AND GET NODE LOCATION ATTRIBUTES
    #---------------------------------------------------------------------------

    # Co-ordinate type variable
    if landscape_grid_data.coord_type != 1 && landscape_grid_data.coord_type != 2 && landscape_grid_data.coord_type != 3     #Check value, throw error if invalid
        error("coord_type has value $coord_type, but coord_type must take value 1 ('Cartesian (metres)'), 2 ('Cartesian (km)') or 3 ('LatLong').")
    end

    # Holding location attributes
    n_holdings = length(holding_info_vec)::Int64
    holding_loc_X_vals::Array{Float64,1} = getfield.(holding_info_vec,:X_loc)
    holding_loc_Y_vals::Array{Float64,1} = getfield.(holding_info_vec,:Y_loc)
    holding_loc_X_and_Y_vals = hcat(holding_loc_X_vals,holding_loc_Y_vals)

    # Holding livestock counts
    n_livestock_types::Int64 = length(holding_info_vec[1].livestock_count)
    if n_livestock_types == 1
        holding_livestock_data = zeros(Int64,n_holdings) #Get number of species/livestock types in use
        for holding_idx = 1:n_holdings
            holding_livestock_data[holding_idx] = holding_info_vec[holding_idx].livestock_count[1]
        end
    elseif n_livestock_types > 1
        holding_livestock_data = zeros(n_holdings,n_livestock_types) #Get number of species/livestock types in use
        for holding_idx = 1:n_holdings
            holding_livestock_data[holding_idx,:] = holding_info_vec[holding_idx].livestock_count
        end
    else
        # If n_livestock_types not positive, throw an error
        error("n_livestock_types has value $n_livestock_types. Invalid.")
    end

    #---------------------------------------------------------------------------
    # GET NODE TRANSMISSIBILITY & SUSCEPTIBILITY
    # ALSO GET SPECIES LEVEL TRANSMISSIBILITY & SUSCEPTIBILITY (PER NODE)
    #---------------------------------------------------------------------------

    # Call function to calculate holding susceptibilities and transmissability
    # that will be used to reinitialise variables in each replicate
    holding_suscept_start_val::Vector{Float64},
    holding_transmiss_start_val::Vector{Float64},
    holding_suscept_by_livestock_type_start_val::Vector{Float64},
    holding_transmiss_by_livestock_type_start_val::Vector{Float64} = function_params.calc_holding_suscep_transmiss_fn(holding_livestock_data,
                                                                                                    epi_params)

    #---------------------------------------------------------------------------
    # SET UP GRID
    #---------------------------------------------------------------------------

    #Call function to get landscape dimensions
    bounding_box_width::Float64,
    bounding_box_height::Float64,
    longest_landscape_edge::Float64,
    max_dist_within_landscape::Int64 = get_landscape_size(landscape_grid_data,
                                                          holding_loc_X_vals,
                                                          holding_loc_Y_vals)

    #---------------------------------------------------------------------------
    # CONSTRUCT KERNEL LOOKUP ARRAY
    #---------------------------------------------------------------------------
    kernel_lookup_vec::Vector{Float64} = function_params.kernel_fn(max_dist_within_landscape)

    #---------------------------------------------------------------------------
    # ASSIGN GRIDID FOR EACH PREMISES. ORDER PREMISES BY GRIDID
    #---------------------------------------------------------------------------

    #Get flag variable corresponding to grid optimisation method in use
    if landscape_grid_data.grid_optim_method != 1 &&  landscape_grid_data.grid_optim_method != 2 #Error check
        error("landscape_grid_data.grid_optim_method has value $(landscape_grid_data.grid_optim_method)). Must take value 1 or 2.")
    end

    #Grid ID assignment
    if landscape_grid_data.grid_optim_method == 1 #fixed

        holding_cell_IDs::Vector{Int64},
        cell_xMin::Vector{Float64},
        cell_xMax::Vector{Float64},
        cell_yMin::Vector{Float64},
        cell_yMax::Vector{Float64} = static_grid_construction(landscape_grid_data.n_cell_per_side::Int64,
                                                                                            landscape_grid_data.bounding_box_vals::Array{Float64,1},
                                                                                            holding_loc_X_vals::Array{Float64,1},
                                                                                            holding_loc_Y_vals::Array{Float64,1})
    elseif landscape_grid_data.grid_optim_method == 2 #dynamic

        #First, estimate number of nodes per grid cell that will optimise efficiency of
        #spatial simulation

        #Requires construction of landscape type variable
        get_trans = maximum(holding_transmiss_start_val)  #Maximum transmissibility
        get_susc = median(holding_suscept_start_val)  #Median susceptibility

        L::LandscapeData = LandscapeData(n_holdings,
                                        longest_landscape_edge,
                                        get_trans,
                                        get_susc,
                                        landscape_grid_data.coord_type)


        #Call node threshold function
        approx_node_per_cell_grid_size_threshold = 50 #Set max num of cells per grid side to test
        approx_node_per_cell_threshold::Float64 = approx_node_per_cell_threshold_fn(L,
                                                                    kernel_lookup_vec,
                                                                    time_params.timestep_val,
                                                                    approx_node_per_cell_grid_size_threshold)
            #Note, could replace with nMax as function input!
        println("approx_node_per_cell_threshold calculated.")

        #Run adaptive grid construction approach.
        holding_cell_IDs,
         cell_xMin,
         cell_xMax,
         cell_yMin,
         cell_yMax = adaptive_grid_construct_fn(approx_node_per_cell_threshold,
                                                                                landscape_grid_data.bounding_box_vals,
                                                                                longest_landscape_edge,
                                                                                holding_loc_X_vals,
                                                                                holding_loc_Y_vals)
        println("Adaptive grid construction completed.")
    end

    #---------------------------------------------------------------------------
    # GET GRID-LEVEL MAXMIUM SUSCEPTIBILITY AND TRANSMISSIBILITY
    # ALSO ASSIGN NUMBER OF LOCATIONS PER GRID TO VARIABLE
    #---------------------------------------------------------------------------
    n_holdings_per_cell::Vector{Int64},
    holding_in_cell_IDs::Vector{Vector{Int64}},
    Max_Sus_grid::Vector{Float64},
    Max_Trans_grid::Vector{Float64} = get_cell_attributes(holding_cell_IDs,
                                            holding_loc_X_and_Y_vals,
                                            holding_suscept_start_val,
                                            holding_transmiss_start_val)

    #---------------------------------------------------------------------------
    # OUTPUT GRID CONFIGURATION TO FILE
    #---------------------------------------------------------------------------
    #println("holding_cell_IDs: $holding_cell_IDs")
    grid_limits_array = [cell_xMin cell_xMax cell_yMin cell_yMax]
    writedlm(save_grid_config_filenames[1],grid_limits_array, ',')
    writedlm(save_grid_config_filenames[2],holding_cell_IDs)
    writedlm(save_grid_config_filenames[3],n_holdings_per_cell)

    #---------------------------------------------------------------------------
    # SET UP BINOMIAL RNG VARIABLES
    #---------------------------------------------------------------------------

    #Pre-set probabilities. Includes 0 as first entry so opening interval is defined
    P_CS = [0,5.0e-9,1.0e-8,1.0e-7,1.0e-6,1.0e-5,1.0e-4,1.0e-3,1.0e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

    #Call function to construct binomial RNG for each value of P_CS and with total
    #number of trials between 1 and maximum holding number in a cell
    binomial_RNG_array::Array{Binomial{Float64},2} = construct_binomial_RNG(n_holdings_per_cell::Array{Int64,1},
                                                                            P_CS::Vector{Float64})

    #---------------------------------------------------------------------------
    # PRECALCULATE GRID TO GRID TRANSMISSION PROBABILITIES
    #---------------------------------------------------------------------------
    max_cell_prob::Array{Float64,2} = max_cell_prob_calc(grid_limits_array,
                                                        n_holdings_per_cell,
                                                        landscape_grid_data.coord_type,
                                                        kernel_lookup_vec,
                                                        Max_Sus_grid,
                                                        Max_Trans_grid,
                                                        time_params.timestep_val)

    #---------------------------------------------------------------------------
    # INITIALISE STORAGE ARRAYS - COUNTS FOR ENTITIES THAT WERE BOTH VACCINATED AND INFECTED
    #---------------------------------------------------------------------------
    n_holdings_both_vacc_and_infected = zeros(Int64,n_reps)
    n_livestock_both_vacc_and_infected = zeros(Int64,n_reps,n_livestock_types)

    #---------------------------------------------------------------------------
    # THE MAIN ITERATION. RUN OUTBREAK. OUTPUT DESIRED STATISTICS PER TIMESTEP
    # BASED ON KEELING&ROHANI BOOK, PROGRAM_7_6.M.
    #---------------------------------------------------------------------------
    for itr_idx = 1:n_reps     #Run n_reps total replicates

        #-----------------------------------------------------------------------
        # SET UP RANDOM NUMBER GENERATOR AT START OF REPLICATE
        #-----------------------------------------------------------------------
        rng = MersenneTwister(RNGseed + itr_idx)
        println(rng)

        #-----------------------------------------------------------------------
        # REINITIALISE STATUS VARIABLES
        #-----------------------------------------------------------------------

        #Initialise holding_status, holding_has_had_infection_flag and holding_vacc_status
        holding_status = zeros(n_holdings)
        holding_has_had_infection_flag = zeros(Bool, n_holdings)
        holding_vacc_status = zeros(Bool, n_holdings)
        livestock_type_vacc_status_by_holding = zeros(Bool, n_holdings)

        # Reset holding-level and livestock-level susceptiblity and transmissibility
        # to unmodified values (when no controls applied)
        holding_suscept = copy(holding_suscept_start_val)
        holding_transmiss = copy(holding_transmiss_start_val)
        holding_suscept_by_livestock_type = copy(holding_suscept_by_livestock_type_start_val)
        holding_transmiss_by_livestock_type = copy(holding_transmiss_by_livestock_type_start_val)

        #-----------------------------------------------------------------------
        # SEED INFECTION
        #-----------------------------------------------------------------------

        #Initialise disease status indicator array and event time array.
        #Events depends on control function in use
        if (control_params.run_controls_fn == control_housing_and_vacc_fn!)  # Vaccination and housing controls in use
            states_to_track = 7
                # Movements between each disease state from latently infected onward, +4
                # & implementation of vaccination control (+2, vaccination administered time
                #                                              &  vaccination becomes effective time)
                # & implementation of housing (+1)
        else #Culling/infectious period ending only
            states_to_track = 4
            # Movements between each disease state from latently infected onward (+4)
        end
        event_indicator_array = [1:1:n_holdings zeros(Int64,n_holdings,states_to_track)] #First column is holding ID.
        event_time_array = [1:1:n_holdings  -1 .*ones(Float64,n_holdings,states_to_track)]

        #Call function to seed infection
        seed_holdings_IDs::Vector{Int64} = seed_infection(initial_infection_info,
                                            holding_info_vec,
                                            holding_loc_X_and_Y_vals,
                                            landscape_grid_data.coord_type,
                                            rng)
        println("seed_holdings_IDs: $seed_holdings_IDs")

        #Use seed_holdings_IDs to update holding infection status
        if ndims(seed_holdings_IDs) == 0 #Only a single holding status needs updating
            holding_status[seed_holdings_IDs] = time_params.timestep_val  #Update current disease status value
            holding_has_had_infection_flag[seed_holdings_IDs] = 1 #Update record of infection ever being present
        elseif ndims(seed_holdings_IDs) == 1 #Multiple holding statuses to be updated
            holding_status[seed_holdings_IDs] .= time_params.timestep_val
            holding_has_had_infection_flag[seed_holdings_IDs] .= 1
            println("holding_status[seed_holdings_IDs]: $(holding_status[seed_holdings_IDs]))")
        else
            error("seed_holdings_IDs should have dimension 0 or 1. Returns dimension $(ndims(seed_holdings_IDs)).")
        end

        #-----------------------------------------------------------------------
        # REINITIALISE BEHAVIOURAL GROUP VARIABLES
        #-----------------------------------------------------------------------
        if control_params.run_controls_fn == control_housing_and_vacc_fn!
            for holding_idx = 1:n_holdings
                # Reset intervention usage statuses
                holding_info_vec[holding_idx].house_livestock_flag = false
                holding_info_vec[holding_idx].house_livestock_status = false
                holding_info_vec[holding_idx].vacc_status = false
            end

            # Regenerate behavioural group assignment
            behaviour_group_config_fn(rng,holding_info_vec)
        end

        #-----------------------------------------------------------------------
        # If applicable, assign vaccination stage to holding
        # & initialise specified proportion of holding as vaccinated if appropriate
        #-----------------------------------------------------------------------
        if (control_params.run_controls_fn == control_housing_and_vacc_fn!)

            # Sample the inoculation times for each holding
            # Time between vaccination being administered and becoming effective (if not infected)
            control_params.holding_time_to_inoculation_fn(rng,
                                                            holding_info_vec)

            # For those assigned as applying intervention prior to simulation start
            # i.e. intervention_threshold_distance = 1e10
            # update intervention related variables (if not assigned as a seed infected)
            for holding_idx = 1:n_holdings
                if (0<=holding_status[holding_idx]<epi_params.detection_time) &&
                        (holding_info_vec[holding_idx].intervention_threshold_distance == 1e10) &&
                        (holding_idx ∉ seed_holdings_IDs)
                    # Holding to be vaccinated from start of simn
                    # AND holding is not a seed infected (if seed infected, assume did not have vaccine)

                    # Vaccinated holding, update holding_vacc_status
                    holding_vacc_status[holding_idx] = true
                    holding_info_vec[holding_idx].vacc_status = true

                    # Amend holding_remaining_time_to_vacc_becoming_effective to be zero
                    holding_info_vec[holding_idx].holding_remaining_time_to_vacc_becoming_effective = 0

                    # Update livestock_type_vacc_status_by_holding
                    livestock_type_vacc_status_by_holding[holding_idx] = true

                    # Check holding is susceptible & not a seed infected
                    # If holding is a seed infected, interventions are not applied
                    if (holding_status[holding_idx] == 0)
                        # Adjustments to livestock type suscepbiility and transmissibility due to vaccination
                        holding_suscept_by_livestock_type[holding_idx,:] = holding_suscept_by_livestock_type[holding_idx,:].*(1-control_params.vacc_efficacy)
                        holding_transmiss_by_livestock_type[holding_idx,:] = holding_transmiss_by_livestock_type[holding_idx,:].*(1-control_params.vacc_efficacy)

                        # If applicable, adjustments due to housing livestock
                        if holding_info_vec[holding_idx].house_livestock_flag == true
                            # Error check. House livestock status should be false
                            if holding_info_vec[holding_idx].house_livestock_status == true
                                error("holding_ID $(holding_idx) already has house_livestock_status of true. Invalid, all holdings should have house livestock status of false.")
                            end

                            # Adjustments by livestock type to suscepbiility and transmissibility
                            holding_suscept_by_livestock_type[holding_idx,:] = holding_suscept_by_livestock_type[holding_idx,:].*(1-control_params.house_livestock_effectiveness)
                            holding_transmiss_by_livestock_type[holding_idx,:] = holding_transmiss_by_livestock_type[holding_idx,:].*(1-control_params.house_livestock_effectiveness)

                            # Update housing livestock status
                            holding_info_vec[holding_idx].house_livestock_status = true
                        end

                        # Update holding-level susceptibility and transmissibility
                        holding_suscept[holding_idx] = sum(holding_suscept_by_livestock_type[holding_idx,:])
                        holding_transmiss[holding_idx] = sum(holding_transmiss_by_livestock_type[holding_idx,:])
                    end

                end
            end
        end

        #-----------------------------------------------------------------------
        # RUN REPLICATE
        #-----------------------------------------------------------------------
        #Profile.clear() #Clear profiler

        # Produce HoldingVectorsAndArrays composite type
        holding_vectors_and_arrays_params::HoldingVectorsAndArrays = HoldingVectorsAndArrays(holding_cell_IDs = holding_cell_IDs,
                                                                             holding_loc_X_vals = holding_loc_X_vals,
                                                                             holding_loc_Y_vals = holding_loc_Y_vals,
                                                                             holding_loc_X_and_Y_vals = holding_loc_X_and_Y_vals,
                                                                             holding_livestock_data = holding_livestock_data,
                                                                             holding_suscept = holding_suscept,
                                                                             holding_transmiss = holding_transmiss,
                                                                             holding_suscept_by_livestock_type = holding_suscept_by_livestock_type,
                                                                             holding_transmiss_by_livestock_type = holding_transmiss_by_livestock_type,
                                                                             holding_status = holding_status,
                                                                             holding_has_had_infection_flag = holding_has_had_infection_flag,
                                                                             holding_vacc_status = holding_vacc_status,
                                                                             livestock_type_vacc_status_by_holding = livestock_type_vacc_status_by_holding,
                                                                             holding_infectious_flag_vec = zeros(Bool, n_holdings),
                                                                             cull_holding_during_current_timestep_vec = zeros(Bool, n_holdings),
                                                                             vacc_holding_during_current_timestep_vec = zeros(Bool, n_holdings),
                                                                             vacc_becomes_effective_during_current_timestep_vec = zeros(Bool, n_holdings),
                                                                             livestock_housed_during_current_timestep_vec = zeros(Bool, n_holdings))

        # Run simulation replicate
        @time function_params.run_simn_replicate_fn(rng,
                                                    batchID,
                                                    itr_idx,
                                                    max_cell_prob,
                                                    epi_params,
                                                    control_params,
                                                    landscape_grid_data.coord_type,
                                                    holding_in_cell_IDs,
                                                    holding_info_vec,
                                                    holding_vectors_and_arrays_params,
                                                    event_indicator_array,
                                                    event_time_array,
                                                    time_params,
                                                    kernel_lookup_vec,
                                                    replicate_file_prefix,
                                                    output_file_objs,
                                                    binomial_RNG_array,
                                                    P_CS,
                                                    function_params)
    end

    #Close files (for files storing data across all n_reps)
    for file_idx = 1:length(output_file_objs)
        close(output_file_objs[file_idx])
    end

    return nothing
end


#-------------------------------------------------------------------------------
# LOAD & RUN DESIRED CONFIGURATION
#-------------------------------------------------------------------------------
"""
    run_spatial_grid_simn(grid_config_fn::Function,
                        behaviour_group_config_fn::Function,
                        initial_infection_info::Vector{Any},
                        batchID::String,
                        n_reps::Int64,
                        RNGseed::Int64,
                        config_txt_file_name::String,
                        config_JLD2_file_name::String)

Load and run spatial simulation on selected configuration using the grid method.

Inputs:
- `grid_config_fn::Function`: Load desired variable configuration from a premade function
- `behaviour_group_config_fn::Function`: Load desired behavioural group configuration from a premade function
- `initial_infection_info::Vector{Any}`: (tuple) [seed_method,NumOfNodes/NodeIDs]
                             Seed method.
                              - 1 = random,
                              - 2/3 = single/group specific node id(s),
                              - 4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
                              - 5 = seed a random site and it's N nearest neighbours.
                              - 6 = seed a random site and 2 nearest neighbours
                             Number of nodes to seed each replicate if seed_method = 1 or 5; or if seed_method = 2/3, seed this specific node every replicate.
                                seed_method = 4 from file, node ids to seed given by one id/row, number of lines must be == number of replicates.
                                seed method = 6, list of IDs corresponds to counties/unitary authorities where index infected node may be selected
- `batchID::String`: String to be used in filenames as an identifier
- `n_reps::Int64`: Number of replicates to run using selected configuration
- `RNGseed::Int64`: Input to seed RNG
- `config_txt_file_name::String`: Text file to keep record of variable values used in each batchID
- `config_JLD2_file_name::String`: JLD2 file to keep record of variable values used in each batchID

Outputs: None \n
Location: run\\_generic\\_livestock\\_disease\\_spatial\\_model.jl
"""
function run_spatial_grid_simn(grid_config_fn::Function,
                            behaviour_group_config_fn::Function,
                            initial_infection_info::Vector{Any},
                            batchID::String,
                            n_reps::Int64,
                            RNGseed::Int64,
                            config_txt_file_name::String,
                            config_JLD2_file_name::String)

    # Initialise RNG seed
    Random.seed!(RNGseed)
    rng = MersenneTwister(RNGseed)

    # Load variable configuration
    holding_info_vec::Vector{HoldingInfo},
    landscape_grid_data::GridData,
    time_params::TimeParams,
    epi_params::EpiParams,
    control_params::ControlParams,
    function_params::FunctionCollection,
    save_grid_config_filenames::Array{String,1},
    replicate_file_prefix::String,
    output_file_objs::Array{IOStream,1} = grid_config_fn(rng)

    # Assign behavioural groups and intervention response types
    behaviour_group_config_fn(rng,holding_info_vec)

    println("Loaded variable configuration.")

    # Save config details to txt file
    io = open(config_txt_file_name, "w")
    writedlm(io, ["grid_config_fn: $grid_config_fn"])
    writedlm(io, ["batchID: $batchID"])
    writedlm(io, ["n_reps: $n_reps"])
    writedlm(io, ["RNGseed: $RNGseed"])
    writedlm(io, ["time_params: $time_params"])
    writedlm(io, ["landscape_grid_data: $landscape_grid_data"])
    writedlm(io, ["epi_params: $epi_params"])
    writedlm(io, ["control_params: $control_params"])
    writedlm(io, ["function_params: $function_params"])
    writedlm(io, ["behaviour_group_config_fn: $behaviour_group_config_fn"])
    writedlm(io, ["replicate_file_prefix: $replicate_file_prefix"])
    close(io)

    # Save config details to JLD 2 file
    jldsave(config_JLD2_file_name;
            grid_config_fn,
            behaviour_group_config_fn,
            batchID,
            n_reps,
            RNGseed,
            holding_info_vec,
            landscape_grid_data,
            time_params,
            epi_params,
            control_params,
            function_params,
            save_grid_config_filenames,
            replicate_file_prefix,
            output_file_objs)

    #Call model simulation function
    spatial_livestock_grid_infection_model(RNGseed,
                                            batchID,
                                            n_reps,
                                            holding_info_vec,
                                            landscape_grid_data,
                                            behaviour_group_config_fn,
                                            time_params,
                                            epi_params,
                                            control_params,
                                            function_params,
                                            save_grid_config_filenames,
                                            replicate_file_prefix,
                                            output_file_objs)
    return nothing
end


#-------------------------------------------------------------------------------
# SET VARIABLES FROM ARGS
#-------------------------------------------------------------------------------
args = ARGS

# If running locally from REPL, will not have any input ARGS
# Can set values for input parameters to main run function here
# args/ARGS list
# args[1] job_ID
# args[2] batchID_offset: Value that batchID is offset by
# args[3] RNGseed: To be used to initialise the random number generator
# args[4] grid_config_fn: Location and pathogen configuration
# args[5] behaviour_group_config_fn: Behavioural response configuration
# args[6] seed_infection_params: Seed infection method (see "seed_infection" in "supporting_fns_transmission_dynamics.jl" for options)
# args[7] n_reps: Number of replicates requested
if length(ARGS)==0
    args = [ "1", "0", "1234", "GB_grid_config", "assign_three_behavioural_groups_empirical_data!", "[6,[100]]", "3"]
end

# Set identifier for job
const job_ID = parse(Int64, args[1])

# Set value that batchID will be offset by
const batchID_offset = parse(Int64, args[2])

# Set RNG seed
const RNGseed = parse(Int64, args[3])

# # Specify configuration function
# # Relevant files in "../GridSimnFns/SimnVarConfigs.jl"
s_config = Symbol(args[4])
const grid_config_fn = getfield(Main, s_config) #Make Symbol a callable function

# Specify behavioural group configuration to apply to population
s_behaviour_groups = Symbol(args[5])
const behaviour_group_config_fn = getfield(Main, s_behaviour_groups) #Make Symbol a callable function

# Assign seed infection parameters to variable
initial_infection_info = eval(Meta.parse(args[6]))

# Set number of replicates to be run
const n_reps = parse(Int64, args[7])

# Assign batchID from list of possible batchIDs
const batchID_vec = string.(batchID_offset .+ collect(1:100))
const batchID = batchID_vec[job_ID]

# If seed infection method based on region, set region IDs for seed infection to 0
# With dummy data, all premises can be selected for initial infection
if initial_infection_info[1] == 6
    # Load regions to seed infection based on job_ID
    initial_infection_info[2] = [0]
        # seed_region_scenarios pads entries with "" to give an array.
        # Therefore, for selected row remove the "" elements to leave the integers corresponding to region IDs

    println("Region IDs with seed infections: $(initial_infection_info[2])")
end

# File to store variable values
if (grid_config_fn == GB_grid_config) || (grid_config_fn == GB_grid_config_IP_control_only)
    config_txt_file_name = string("config_log_files/grid_simn/GB_model_grid_simn_config_batch_ID",batchID,".txt")
    config_JLD2_file_name = string("config_JLD2_files/grid_simn/GB_model_grid_simn_config_batch_ID",batchID,".jld2")
else
    error("Invalid grid_config_fn provided.")
end

# Call the main function
run_spatial_grid_simn(grid_config_fn,
                    behaviour_group_config_fn,
                    initial_infection_info,
                    batchID,
                    n_reps,
                    RNGseed,
                    config_txt_file_name,
                    config_JLD2_file_name)
