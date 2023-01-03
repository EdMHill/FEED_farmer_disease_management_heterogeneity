#=
Purpose:
This file houses functions to perform single replicate of spatial simulation
of livestock infectious disease.
Version that uses the conditional subsample grid algorithm
=#

#-------------------------------------------------------------------------------
### Update holding status following infection updates & before control carried out
#-------------------------------------------------------------------------------
"""
    update_holding_status!(holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                            delta_t::Float64)

Update the disease status of each holding

Inputs:
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: supporting\\_fns\\_run\\_model\\_grid\\_simns.jl
"""
function update_holding_status!(holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                delta_t::Float64)

    for holding_idx_itr = 1:length(holding_vectors_and_arrays_params.holding_status)
        if holding_vectors_and_arrays_params.holding_status[holding_idx_itr]>0 #Premises previously infected. Update status
            holding_vectors_and_arrays_params.holding_status[holding_idx_itr] = holding_vectors_and_arrays_params.holding_status[holding_idx_itr] + delta_t
        elseif holding_vectors_and_arrays_params.holding_status[holding_idx_itr]==0 #Premises susceptible. Add infection flag value. Those infected on current timestep will update
            holding_vectors_and_arrays_params.holding_status[holding_idx_itr] = holding_vectors_and_arrays_params.holding_status[holding_idx_itr] + (holding_vectors_and_arrays_params.holding_has_had_infection_flag[holding_idx_itr]*delta_t)
        end
    end
    return nothing
end

#-------------------------------------------------------------------------------
### Perform single timestep of outbreak
#-------------------------------------------------------------------------------
"""
    iterate_outbreak_grid_simn!(rng::AbstractRNG,
                                binomial_RNG_array::Array{Binomial{Float64},2},
                                P_CS::Array{Float64,1},
                                n_cells::Int64,
                                max_grid_prob::Array{Float64,2},
                                epi_params::EpiParams,
                                control_params::ControlParams,
                                n_holdings::Int64,
                                holding_info_vec::Array{HoldingInfo,1},
                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                suscept_holding_by_cell_IDs::Array{Array{Int64,1},1},
                                infectious_holding_by_cell_IDs::Array{Array{Int64,1},1},
                                notified_holding_timeseries::Array{Int64,1},
                                coord_type::Int64,
                                kernel_lookup_vec::Array{Float64,1},
                                delta_t::Float64)

Perform a single timestep of the spatial outbreak simulation using the grid method.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `binomial_RNG_array::Array{Binomial{Float64},2}`: Array of Binomial RNGs. Row per susceptible number, column per probability threshold.
- `P_CS::Array{Float64,1})`: Preset probabilities to initialise RNG with.
- `n_cells::Int64`: Number of cells in the entire grid.
- `max_grid_prob::Array{Float64,2}`: Array of grid to grid transmission probabilities.
- `epi_params::EpiParams`: Composite type containing epimidological parameters.
- `control_params::ControlParams`: Variables related to implementing control measures
- `n_holdings::Int64`: Number of holding in landscape.
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `suscept_holding_by_cell_IDs::Array{Array{Int64,1},1}`: Vector of vectors. Vector per cell, with a vector of susceptible holding in that cell.
- `infectious_holding_by_cell_IDs::Array{Array{Int64,1},1}`: Vector of vectors. Vector per cell, with a vector of infected holding in that cell.
- `notified_holding_timeseries::Array{Int64,1}`: Per timestep, the number of holding that reported infection
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep between each iteration.

Outputs:
- `cull_holding_during_current_timestep_vec::Array{Bool,1}`:  Indicator vector, 1 if culling occurs on holding during current timestep, 0 otherwise.
- `livestock_housed_during_current_timestep_vec::Array{Bool,1}`: Indicator vector. 1 if livestock at holding housed during current timestep. 0 otherwise.
- `vacc_holding_during_current_timestep_vec::Array{Bool,1}`: Indicator vector. 1 if holding culled during current timestep. 0 otherwise.
- `vacc_becomes_effective_during_current_timestep_vec::Array{Bool,1}`: Indicator vector, 1 if vaccination now due to become effective during current timestep, 0 otherwise.

Location: supporting\\_fns\\_run\\_model\\_grid\\_simns.jl
"""
function iterate_outbreak_grid_simn!(rng::AbstractRNG,
                                        binomial_RNG_array::Array{Binomial{Float64},2},
                                        P_CS::Array{Float64,1},
                                        n_cells::Int64,
                                        max_grid_prob::Array{Float64,2},
                                        epi_params::EpiParams,
                                        control_params::ControlParams,
                                        n_holdings::Int64,
                                        holding_info_vec::Array{HoldingInfo,1},
                                        holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                        suscept_holding_by_cell_IDs::Array{Array{Int64,1},1},
                                        infectious_holding_by_cell_IDs::Array{Array{Int64,1},1},
                                        notified_holding_timeseries::Array{Int64,1},
                                        coord_type::Int64,
                                        kernel_lookup_vec::Array{Float64,1},
                                        delta_t::Float64)

    #---------------------------------------------------------------------------
    # RUN TRANSMISSION PROCESS
    #---------------------------------------------------------------------------

    # Find infectious grids. Remove any repeat entries
    cell_IDs_with_infectious_holdings = unique(holding_vectors_and_arrays_params.holding_cell_IDs[holding_vectors_and_arrays_params.holding_infectious_flag_vec.==1])
    n_infectious_cells = length(cell_IDs_with_infectious_holdings)

    if n_infectious_cells == 0 # If no infectious cells, update status for those exposed

        # Update status for those previously infected
        for prem_idx_itr = 1:n_holdings
            if holding_vectors_and_arrays_params.holding_status[prem_idx_itr]>0 # Premises previously infected. Update status
                holding_vectors_and_arrays_params.holding_status[prem_idx_itr] = holding_vectors_and_arrays_params.holding_status[prem_idx_itr] + delta_t
            end
        end
    else #Infectious holding present.

        # Iterate over each cell with infectious holding and check for infection events
        for infectious_cell_itr_idx = 1:n_infectious_cells
            infectious_cell_ID = cell_IDs_with_infectious_holdings[infectious_cell_itr_idx]
                #Get ID of cell containing infected units

            #Have infectious nodes in cell A. Iterate over each cell B
            #When A and B are not the same cell, use conditional subsample algorithm
            #When A and B are the same cell, use pairwise algorithm
            for suscept_cell_itr_idx = 1:n_cells
                n_sus::Int64 = length(suscept_holding_by_cell_IDs[suscept_cell_itr_idx])
                if n_sus > 0 #If no susceptibles can move onto next cell
                    if infectious_cell_ID == suscept_cell_itr_idx #use pairwise algorithm
                        local_pairwise_check_sucep_infectious_paris_in_cell!(rng,
                                                                            n_sus,
                                                                            holding_vectors_and_arrays_params,
                                                                            suscept_holding_by_cell_IDs[suscept_cell_itr_idx],
                                                                            infectious_holding_by_cell_IDs[infectious_cell_ID],
                                                                            coord_type,
                                                                            kernel_lookup_vec,
                                                                            delta_t)

                    else     #use conditional subsample algorithm
                        max_grid_prob_val = max_grid_prob[infectious_cell_ID,suscept_cell_itr_idx] #Cell-to-cell (over-estimated) rate of transmission

                        #If max_grid_prob_val = 0, no infection possible. Skip remainder of function
                        #Otherwise, perform infection event check.
                        if max_grid_prob_val > 0
                            conditional_subsample_alg!(rng,
                                                        n_sus,
                                                        max_grid_prob_val,
                                                        holding_vectors_and_arrays_params,
                                                        suscept_holding_by_cell_IDs[suscept_cell_itr_idx],
                                                        infectious_holding_by_cell_IDs[infectious_cell_ID],
                                                        binomial_RNG_array,
                                                        P_CS,
                                                        coord_type,
                                                        kernel_lookup_vec,
                                                        delta_t)
                        end
                    end
                end
            end
        end

        # End of timestep, pre-control: Update status for holding
        update_holding_status!(holding_vectors_and_arrays_params,
                                delta_t)
    end

    #---------------------------------------------------------------------------
    # RUN CONTROLS
    #---------------------------------------------------------------------------

    #Initiate controls at end of timestep
    holding_vectors_and_arrays_params.cull_holding_during_current_timestep_vec .= zeros(Bool,n_holdings)
    if (control_params.run_controls_fn == control_housing_and_vacc_fn!)
        # Housing and vaccination based strategy.

        # Initialise flag vectors for housing and vaccination being administered
        # during current timestep
        holding_vectors_and_arrays_params.livestock_housed_during_current_timestep_vec .= zeros(Bool,n_holdings)
        holding_vectors_and_arrays_params.vacc_holding_during_current_timestep_vec .= zeros(Bool,n_holdings)

        #Initialise flag vector for vaccination becoming effective current timestep
        holding_vectors_and_arrays_params.vacc_becomes_effective_during_current_timestep_vec .= zeros(Bool,n_holdings)

        #Run control function
        control_params.run_controls_fn(holding_info_vec,
                        n_holdings,
                        holding_vectors_and_arrays_params,
                        epi_params,
                        control_params,
                        coord_type,
                        delta_t)
    else # End of IP infectious period only
        control_params.run_controls_fn(n_holdings,
                        holding_vectors_and_arrays_params,
                        epi_params)
    end

    return nothing
end

#-------------------------------------------------------------------------------
### Perform complete outbreak replicate
#-------------------------------------------------------------------------------
"""
    run_outbreak_grid_simn(rng::AbstractRNG,
                                    batchID::String,
                                    itr_idx::Int64,
                                    max_grid_prob::Array{Float64,2},
                                    epi_params::EpiParams,
                                    control_params::ControlParams,
                                    coord_type::Int64,
                                    holding_in_cell_IDs::Array{Array{Int64,1},1},
                                    holding_info_vec::Array{HoldingInfo,1},
                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                    event_indicator_array::Array{Int64,2},
                                    event_time_array::Array{Float64,2},
                                    time_params::TimeParams,
                                    kernel_lookup_vec::Array{Float64,1},
                                    replicate_file_prefix::String,
                                    output_file_objs::Array{IOStream,1},
                                    binomial_RNG_array::Array{Binomial{Float64},2},
                                    P_CS::Array{Float64,1},
                                    function_params::FunctionCollection)

Perform a complete replicate of the spatial outbreak simulation using the grid method.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `batchID::String`: Used as prefix for output files.
- `itr_idx::Int64`: ID for replicate
- `max_grid_prob::Array{Float64,2}`: Array of grid to grid transmission probabilities.
- `epi_params::EpiParams`: Composite type containing epimidological parameters.
- `control_params::ControlParams`: Variables related to implementing control measures
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `holding_in_cell_IDs::Array{Array{Int64,1},1}`: For each cell, vector of node IDs that reside in that cell.
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `holding_vectors_and_arrays_params::Array{Float64,1}`: Premises-level disease status. Entry per holding. Susceptible: 0. infected: >0. Increases by 1 each day infected until culled.
- `event_indicator_array::Array{Int64,2}`: Flag recording if event occurs. Row per holding. First column gives holding ID.
- `event_time_array::Array{Float64,2}`: Array recording time the event occurred at. Row per holding. First column gives holding ID.
- `time_params::TimeParams`: Timestep per iteration & timeframe simulation is run over.
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `replicate_file_prefix::String`: Directory location to be used with output files storing data with individual file per replicate
- `output_file_objs::Array{IOStream,1}`: Filename identifiers written to by all replicates.
- `binomial_RNG_array::Array{Binomial{Float64},2}`: Array of Binomial RNGs. Row per susceptible number, column per probability threshold.
- `P_CS::Array{Float64,1})`: Preset probabilities to initialise RNG with.
- `function_params::FunctionCollection`: Composite type containing collection of functions used in model simulation.

Outputs: None \n
Location: supporting\\_fns\\_run\\_model\\_grid\\_simns.jl
"""
function run_outbreak_grid_simn(rng::AbstractRNG,
                                batchID::String,
                                itr_idx::Int64,
                                max_grid_prob::Array{Float64,2},
                                epi_params::EpiParams,
                                control_params::ControlParams,
                                coord_type::Int64,
                                holding_in_cell_IDs::Array{Array{Int64,1},1},
                                holding_info_vec::Array{HoldingInfo,1},
                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                event_indicator_array::Array{Int64,2},
                                event_time_array::Array{Float64,2},
                                time_params::TimeParams,
                                kernel_lookup_vec::Array{Float64,1},
                                replicate_file_prefix::String,
                                output_file_objs::Array{IOStream,1},
                                binomial_RNG_array::Array{Binomial{Float64},2},
                                P_CS::Array{Float64,1},
                                function_params::FunctionCollection)

    #---------------------------------------------------------------------------
    #Intialise variables to be used throughout script
    #---------------------------------------------------------------------------

    #Compute number of occupied grids
    #Equal to number of rows/columns of cell-to-cell transmission array
    n_cells = size(max_grid_prob,1)

    #Get number of holding in use
    n_holdings = size(holding_vectors_and_arrays_params.holding_status,1)

    #Get number of livestock species in use
    n_livestock_types = size(holding_vectors_and_arrays_params.holding_livestock_data,2)

    #---------------------------------------------------------------------------
    #Open files to store outputs, for specific replicate
    #---------------------------------------------------------------------------
    event_indicator_array_file = string(replicate_file_prefix,"_event_arrays/event_indicator_array_batchID$(batchID)_replicate$(itr_idx).txt")
    event_time_array_file = string(replicate_file_prefix,"_event_arrays/event_time_array_batchID$(batchID)_replicate$(itr_idx).txt")

    holding_per_disease_state_per_timestep_file = string(replicate_file_prefix,"_per_timestep/holding_per_disease_state_per_timestep_batchID$(batchID)_replicate$(itr_idx).txt")
    cumulative_culled_per_timestep_file = string(replicate_file_prefix,"_per_timestep/cumulative_culled_per_timestep_batchID$(batchID)_replicate$(itr_idx).txt")
    cumulative_vacc_per_timestep_file = string(replicate_file_prefix,"_per_timestep/cumulative_vacc_per_timestep_batchID$(batchID)_replicate$(itr_idx).txt")
    cumulative_livestock_housed_per_timestep_file = string(replicate_file_prefix,"_per_timestep/cumulative_livestock_housed_per_timestep_batchID$(batchID)_replicate$(itr_idx).txt")
    cumulative_cases_livestock_per_timestep_file = string(replicate_file_prefix,"_per_timestep/cumulative_cases_livestock_per_timestep_batchID$(batchID)_replicate$(itr_idx).txt")

    #---------------------------------------------------------------------------
    #Initialise simulation tracking variables
    #---------------------------------------------------------------------------
    t = 0.
    simn_itr = 1

    #---------------------------------------------------------------------------v
    #Initialise compartment status variables
    #---------------------------------------------------------------------------
    req_itr = convert(Int64,ceil(time_params.max_time/time_params.timestep_val) + 1)  #Add 1 as also storing value at intial time, then every timestep_val until max_time exceeded
    println("req_itr: $req_itr")
    S = zeros(Int64, req_itr); E = zeros(Int64, req_itr); I = zeros(Int64, req_itr);
    R = zeros(Int64, req_itr); R2 = zeros(Int64, req_itr); Culled = zeros(Int64, req_itr);

    # Premises level storage vectors
    vacc_holding = zeros(Int64, req_itr)
    n_holding_housed_livestock = zeros(Int64, req_itr)
    cumulative_cases_holding = zeros(Int64, req_itr)
    cumulative_cases_without_control_holding = zeros(Int64, req_itr)
    cumulative_cases_with_control_holding = zeros(Int64, req_itr)

    # Livestock type level storage vectors
    culled_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)
    vacc_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)
    housed_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)
    cumulative_cases_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)
    cumulative_cases_without_control_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)
    cumulative_cases_with_control_by_livestock_type = zeros(Int64,req_itr,n_livestock_types)

    #---------------------------------------------------------------------------
    #Assign first entry to compartment status variables
    #---------------------------------------------------------------------------
    #Equivalent to:
        # S[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.==0)
        # E[simn_itr] = sum(0 .< holding_vectors_and_arrays_params.holding_status .<=incubation_time)
        # I[simn_itr] = sum(incubation_time .<holding_vectors_and_arrays_params.holding_status .<=detection_time)
        # R[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.==(detection_time+1))
        # R2[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.>detection_time)
        # Culled[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.<0)
        # culled_by_livestock_type[simn_itr,:] = zeros(Int64,size(holding_livestock_data,2))

    # Initialise temporary storage vectors
    logic_temp_bit_vec = BitArray{1}(undef,n_holdings)
    holding_suscept_flag_vec = BitArray{1}(undef,n_holdings)

    # Update node type variables & suscept flag value
    for holding_idx = 1:n_holdings

        # Reinitialise holding_suscept_flag_vec entry & holding_infectious_flag_vec
        holding_suscept_flag_vec[holding_idx] = 0
        holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 0

        # Check holding status and update appropriate variables
        if holding_vectors_and_arrays_params.holding_status[holding_idx] < 0
            # Removed
            Culled[simn_itr] += 1

            # Increment culled by livestock type array
            for livestock_type_idx = 1:n_livestock_types
                culled_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
            end
        elseif holding_vectors_and_arrays_params.holding_status[holding_idx] == 0
            # Susceptible
            S[simn_itr] += 1
            holding_suscept_flag_vec[holding_idx] = 1
        elseif holding_vectors_and_arrays_params.holding_status[holding_idx] <= epi_params.incubation_time
            # Latent infected
            E[simn_itr] += 1
        elseif holding_vectors_and_arrays_params.holding_status[holding_idx] <= epi_params.detection_time
            # Infectious, not reported
            I[simn_itr] += 1
            holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 1
        elseif holding_vectors_and_arrays_params.holding_status[holding_idx] > epi_params.detection_time
            # Infectious, reported
            R2[simn_itr] += 1
            holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 1

            # Also record if report occurred on current timestep
            if holding_vectors_and_arrays_params.holding_status[holding_idx]==(epi_params.detection_time+time_params.timestep_val)
                R[simn_itr] += 1
            end
        else
            error("Invalid holding_status value for holding ID $holding_idx. holding_status[$holding_idx] = $(holding_status[holding_idx]).")
        end

        # Check for new infection events
        # Seed sites, latently infected at initial timestep
        if holding_vectors_and_arrays_params.holding_status[holding_idx] == time_params.timestep_val
            event_indicator_array[holding_idx,2] = 1
            event_time_array[holding_idx,2] = t
        end

        # Check livestock housed status
        if holding_info_vec[holding_idx].house_livestock_status == true
            n_holding_housed_livestock[simn_itr] += 1
            for livestock_type_idx = 1:n_livestock_types
                housed_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
            end

            # Update event arrays
            event_indicator_array[holding_idx,8] = 1
            event_time_array[holding_idx,8] = t
        end

        # Check holding vaccination status
        if holding_info_vec[holding_idx].vacc_status == true
            vacc_holding[simn_itr] += 1
            for livestock_type_idx = 1:n_livestock_types
                vacc_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
            end

            # Update event arrays
            # Record that vaccine has been given
            event_indicator_array[holding_idx,6] = 1
            event_time_array[holding_idx,6] = t

            # Update event arrays
            # Record that vaccine is effective
            event_indicator_array[holding_idx,7] = 1
            event_time_array[holding_idx,7] = t
        end

        # Update cumulative cases statistics
        if holding_vectors_and_arrays_params.holding_has_had_infection_flag[holding_idx] == true
            cumulative_cases_holding[simn_itr] += 1
            for livestock_type_idx = 1:n_livestock_types
                cumulative_cases_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                if holding_info_vec[holding_idx].vacc_status == false
                    cumulative_cases_without_control_holding[simn_itr] += 1
                    cumulative_cases_without_control_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                else
                    cumulative_cases_with_control_holding[simn_itr] += 1
                    cumulative_cases_with_control_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                end
            end
        end
    end

    #---------------------------------------------------------------------------
    #Initialise suscept_holding_by_cell_IDs, infectious_holding_by_cell_IDs, SusPremByCell_Locs, InfPremByCell_Locs
    # Tuples, with each entry giving ID/location of holding that is susceptible/infected
    #---------------------------------------------------------------------------
    suscept_holding_by_cell_IDs::Vector{Vector{Int64}},
    infectious_holding_by_cell_IDs::Vector{Vector{Int64}} =  construct_holding_by_cell_tuples(n_cells,
                                                                        holding_in_cell_IDs,
                                                                        holding_vectors_and_arrays_params.holding_status,
                                                                        holding_vectors_and_arrays_params.holding_infectious_flag_vec,
                                                                        holding_suscept_flag_vec)

    #---------------------------------------------------------------------------
    #Enter main iteration
    #---------------------------------------------------------------------------
    simn_itr = simn_itr + 1
    iterate_flag = 1 # Flag variable that when altered will lead exit from loop
    while (t < time_params.max_time && iterate_flag == 1)
        # Call event update function
        function_params.iterate_outbreak_fn(rng,
                                            binomial_RNG_array,
                                            P_CS,
                                            n_cells,
                                            max_grid_prob,
                                            epi_params,
                                            control_params,
                                            n_holdings,
                                            holding_info_vec,
                                            holding_vectors_and_arrays_params,
                                            suscept_holding_by_cell_IDs,
                                            infectious_holding_by_cell_IDs,
                                            R,
                                            coord_type,
                                            kernel_lookup_vec,
                                            time_params.timestep_val)

        #-----------------------------------------------------------------------
        ### Amend compartment based on status
        #-----------------------------------------------------------------------
        # Equivalent to:
            # S[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.==0)
            # E[simn_itr] = sum(0 .< holding_vectors_and_arrays_params.holding_status .<=incubation_time)
            # I[simn_itr] = sum(incubation_time .<holding_vectors_and_arrays_params.holding_status .<=detection_time)
            # R[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.==(detection_time+1))
            # R2[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.>detection_time)
            # Culled[simn_itr] = sum(holding_vectors_and_arrays_params.holding_status.<0)

        # Update simulation time
        t = t + time_params.timestep_val

        # Update node type variables & suscept flag value
        for holding_idx = 1:n_holdings

            # Reinitialise holding_suscept_flag_vec entry & holding_infectious_flag_vec
            holding_suscept_flag_vec[holding_idx] = 0
            holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 0

            # Check holding status and update appropriate variables
            if holding_vectors_and_arrays_params.holding_status[holding_idx] < 0
                # Removed
                Culled[simn_itr] += 1

                # Increment culled by livestock type array
                for livestock_type_idx = 1:n_livestock_types
                    culled_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                end

                # If required, update event arrays
                if holding_vectors_and_arrays_params.cull_holding_during_current_timestep_vec[holding_idx] == 1
                    event_indicator_array[holding_idx,5] = 1
                    event_time_array[holding_idx,5] = t
                end
            elseif holding_vectors_and_arrays_params.holding_status[holding_idx] == 0
                # Susceptible
                S[simn_itr] += 1
                holding_suscept_flag_vec[holding_idx] = 1
            elseif holding_vectors_and_arrays_params.holding_status[holding_idx] <= epi_params.incubation_time
                # Latent infected
                E[simn_itr] += 1

                # Check for new infection events
                # Seed sites, latently infected at initial timestep
                if holding_vectors_and_arrays_params.holding_status[holding_idx] == time_params.timestep_val
                    event_indicator_array[holding_idx,2] = 1
                    event_time_array[holding_idx,2] = t
                end
            elseif holding_vectors_and_arrays_params.holding_status[holding_idx] <= epi_params.detection_time
                # Infectious, not reported
                I[simn_itr] += 1
                holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 1

                #Infectious during timestep
                if holding_vectors_and_arrays_params.holding_status[holding_idx] == (epi_params.incubation_time+time_params.timestep_val)
                    event_indicator_array[holding_idx,3] = 1
                    event_time_array[holding_idx,3] = t
                end
            elseif holding_vectors_and_arrays_params.holding_status[holding_idx] > epi_params.detection_time
                # Infectious, reported
                R2[simn_itr] += 1
                holding_vectors_and_arrays_params.holding_infectious_flag_vec[holding_idx] = 1

                # Also record if report occurred on current timestep
                if holding_vectors_and_arrays_params.holding_status[holding_idx]==(epi_params.detection_time+time_params.timestep_val)
                    R[simn_itr] += 1
                    event_indicator_array[holding_idx,4] = 1
                    event_time_array[holding_idx,4] = t
                end
            else
                error("Invalid holding_status value for holding ID $holding_idx. holding_status[$holding_idx] = $(holding_status[holding_idx]).")
            end

            # Check livestock housed status
            if holding_info_vec[holding_idx].house_livestock_status == true
                n_holding_housed_livestock[simn_itr] += 1
                for livestock_type_idx = 1:n_livestock_types
                    housed_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                end

                # Updata event arrays, if required
                if holding_vectors_and_arrays_params.livestock_housed_during_current_timestep_vec[holding_idx] == true
                    event_indicator_array[holding_idx,8] = 1
                    event_time_array[holding_idx,8] = t
                end
            end

            # Check holding vaccination status
            if holding_info_vec[holding_idx].vacc_status == true
                vacc_holding[simn_itr] += 1

                for livestock_type_idx = 1:n_livestock_types
                    vacc_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                end

                # Updata event arrays
                # Record that vaccine has been given
                if holding_vectors_and_arrays_params.vacc_holding_during_current_timestep_vec[holding_idx] == 1
                    event_indicator_array[holding_idx,6] = 1
                    event_time_array[holding_idx,6] = t
                end

                # Updata event arrays
                # Record that vaccine is effective
                if holding_vectors_and_arrays_params.vacc_becomes_effective_during_current_timestep_vec[holding_idx] == 1
                    event_indicator_array[holding_idx,7] = 1
                    event_time_array[holding_idx,7] = t
                end
            end

            # Update cumulative cases statistics
            if holding_vectors_and_arrays_params.holding_has_had_infection_flag[holding_idx] == true
                cumulative_cases_holding[simn_itr] += 1
                for livestock_type_idx = 1:n_livestock_types
                    cumulative_cases_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                    if holding_info_vec[holding_idx].vacc_status == false
                        cumulative_cases_without_control_holding[simn_itr] += 1
                        cumulative_cases_without_control_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                    else
                        cumulative_cases_with_control_holding[simn_itr] += 1
                        cumulative_cases_with_control_by_livestock_type[simn_itr,livestock_type_idx] += holding_vectors_and_arrays_params.holding_livestock_data[holding_idx,livestock_type_idx]
                    end
                end
            end
        end

        #-----------------------------------------------------------------------
        #Update increment variables
        #-----------------------------------------------------------------------
        simn_itr = simn_itr + 1

        #-----------------------------------------------------------------------
        #Exit condition
        #-----------------------------------------------------------------------
        if t > 5
            if (E[simn_itr-3] + I[simn_itr-3] + R2[simn_itr-3]) == 0
                println("E[simn_itr-3]: $(E[simn_itr-3])")
                println("I[simn_itr-3]: $(I[simn_itr-3])")
                println("R2[simn_itr-3]: $(R2[simn_itr-3])")
                iterate_flag = 0
                println("Infection no longer present. Simulation terminated.")
            end
        end

        #-----------------------------------------------------------------------
        # Update suscept_holding_by_cell_IDs, infectious_holding_by_cell_IDs, SusPremByCell_Locs, InfPremByCell_Locs
        #-----------------------------------------------------------------------
        suscept_holding_by_cell_IDs,
        infectious_holding_by_cell_IDs = construct_holding_by_cell_tuples(n_cells,
                                                                        holding_in_cell_IDs,
                                                                        holding_vectors_and_arrays_params.holding_status,
                                                                        holding_vectors_and_arrays_params.holding_infectious_flag_vec,
                                                                        holding_suscept_flag_vec)


        if mod(simn_itr,20) .== 0
            println("simn_itr: $simn_itr")
        end
    end

    #-------------------------------------------------------------------
    #Write data to files, replicate specific
    #-------------------------------------------------------------------
    writedlm(holding_per_disease_state_per_timestep_file, [S E I R R2 Culled vacc_holding n_holding_housed_livestock cumulative_cases_holding cumulative_cases_without_control_holding cumulative_cases_with_control_holding])
    writedlm(cumulative_culled_per_timestep_file, culled_by_livestock_type)
    writedlm(cumulative_vacc_per_timestep_file, vacc_by_livestock_type)
    writedlm(cumulative_livestock_housed_per_timestep_file, housed_by_livestock_type)
    writedlm(cumulative_cases_livestock_per_timestep_file, [cumulative_cases_by_livestock_type cumulative_cases_without_control_by_livestock_type cumulative_cases_with_control_by_livestock_type])

    #-------------------------------------------------------------------
    #Write event indicator arrays to file
    #-------------------------------------------------------------------

    # Find those holding (rows) with non-zero entries for events
    # Column 1 is holding ID, so take sum of each row (dims=2), over cols 2:end
    n_cols_in_event_arrays = size(event_indicator_array,2)
    for holding_idx = 1:n_holdings

        # Get number of events to have occurred at the holding
        n_events = 0
        for event_indicator_col_idx = 2:n_cols_in_event_arrays
            n_events += event_indicator_array[holding_idx,event_indicator_col_idx]
        end

        # Update logic vector, based on whether events have occurred at that holding
        if n_events > 0
            logic_temp_bit_vec[holding_idx] = 1
        else
            logic_temp_bit_vec[holding_idx] = 0
        end
    end

    #Write only those rows, corresponding to holding where events occurred, to file
    n_holdings_with_events = sum(logic_temp_bit_vec)
    logic_temp_indicator_array = zeros(Int64,n_holdings_with_events,n_cols_in_event_arrays)
        # Set up temporary array to be allocated to

        # Event indicator arrays
    holding_with_event_idx = 1
    for holding_idx = 1:n_holdings
        if logic_temp_bit_vec[holding_idx] == 1
            for event_idx = 1:size(event_indicator_array,2)
                logic_temp_indicator_array[holding_with_event_idx,event_idx] = event_indicator_array[holding_idx,event_idx]
            end

            # Increement holding_with_event_idx
            holding_with_event_idx += 1
        end
    end
    # logic_temp_indicator_array .= event_indicator_array[logic_temp_bit_vec,:]
    writedlm(event_indicator_array_file, logic_temp_indicator_array)

        # Event time array
    holding_with_event_idx = 1
    for holding_idx = 1:n_holdings
        if logic_temp_bit_vec[holding_idx] == 1
            for event_idx = 1:size(event_indicator_array,2)
                logic_temp_indicator_array[holding_with_event_idx,event_idx] = event_time_array[holding_idx,event_idx]
            end

            # Increement holding_with_event_idx
            holding_with_event_idx += 1
        end
    end
    # logic_temp_indicator_array .= event_time_array[logic_temp_bit_vec,:]
    writedlm(event_time_array_file, logic_temp_indicator_array)

    #-------------------------------------------------------------------
    #Write data to files storing data across all nReps
    #-------------------------------------------------------------------
    writedlm(output_file_objs[1], t)
    writedlm(output_file_objs[2], [S[simn_itr-1] E[simn_itr-1] I[simn_itr-1] R[simn_itr-1] R2[simn_itr-1] Culled[simn_itr-1] vacc_holding[simn_itr-1] n_holding_housed_livestock[simn_itr-1] cumulative_cases_holding[simn_itr-1] cumulative_cases_without_control_holding[simn_itr-1] cumulative_cases_with_control_holding[simn_itr-1]])
    writedlm(output_file_objs[3], culled_by_livestock_type[simn_itr-1])
    writedlm(output_file_objs[4], vacc_by_livestock_type[simn_itr-1])
    writedlm(output_file_objs[5], housed_by_livestock_type[simn_itr-1])
    writedlm(output_file_objs[6], [cumulative_cases_by_livestock_type[simn_itr-1] cumulative_cases_without_control_by_livestock_type[simn_itr-1] cumulative_cases_with_control_by_livestock_type[simn_itr-1]])

    #-------------------------------------------------------------------
    #Return outputs (if any)
    #-------------------------------------------------------------------
    return nothing
end
