#=
Purpose:
File to house supporting functions for running GB model,
associated with interventions

Function list:
- five_day_inoculation_fn!: Five day delay in vaccine becoming effective
- remove_IPs_only_control_fn!
- control_housing_and_vacc_fn!
- apply_livestock_vacc_and_housing_fn!
- cull_holding_and_update_time_to_all_livestock_vacc_effective_fn!(
- activate_livestock_vacc_fn!
=#


#= FUNCTION GENERIC OUTLINE

    inoculation_fn(rng,n_holdings)

Randomly draw sample from specified discrete probability distribution.
Used to randomly set infection, adherence and testing waiting times.

Inputs: `rng` - random number generator,
        `n_holdings` - Number of holdings in the simulation
Outputs: `holding_time_to_inoculation` - For each holding, delay in vaccine becoming effective once administered
=#

#===============================================================================
FIVE_DAY_INOCULATION_FN - FIVE DAY DELAY IN VACCINE BECOMING EFFECTIVE
===============================================================================#
"""
    five_day_inoculation_fn!(rng::AbstractRNG,
                             holding_info_vec::Array{HoldingInfo,1})

Inoculation function: Single day delay in vaccine becoming effective.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function five_day_inoculation_fn!(rng::AbstractRNG,
                                    holding_info_vec::Array{HoldingInfo,1})

    # Update field attributes for inoculation time per holding and initialise
    # time tracker for remaining time until vaccine becomes effective
    setfield!.(holding_info_vec,:holding_vacc_inoclulation_time,5.)
    setfield!.(holding_info_vec,:holding_remaining_time_to_vacc_becoming_effective,5.)

end

#===============================================================================
FNS FOR REMOVING IPs ONLY. NO ADDTIONAL INTERVENTIONS APPLIED
===============================================================================#
"""
    remove_IPs_only_control_fn!(n_holdings::Int64,
                                        holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                        epi_params::EpiParams)

Apply control: Infected holding no longer infectious after specified infectious period duration elapses.

Inputs:
- `n_holdings::Int64`: Number of holding in landscape.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `epi_params::EpiParams`: Composite type containing epimidological parameters..

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function remove_IPs_only_control_fn!(n_holdings::Int64,
                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                    epi_params::EpiParams)

    #Find and update holding to be culled
    holding_infectious_period_ending_IDs = findall(holding_vectors_and_arrays_params.holding_status.>epi_params.removal_time)
    for ii = 1:length(holding_infectious_period_ending_IDs)
        #Get ID of holding to be culled
        holding_to_remove_ID = holding_infectious_period_ending_IDs[ii]

        #Cull IP
        holding_vectors_and_arrays_params.holding_status[holding_to_remove_ID] = -1.

        #Update holding culled during current timestep vector
        holding_vectors_and_arrays_params.cull_holding_during_current_timestep_vec[holding_to_remove_ID] = 1
    end
    return nothing
end

#===============================================================================
FNS FOR REMOVING IPs, USING VACCINATION AND HOUSING OF LIVESTOCK
===============================================================================#
"""
    control_housing_and_vacc_fn!(holding_info_vec::Array{HoldingInfo,1},
                                n_holdings::Int64,
                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                epi_params::EpiParams,
                                control_params::ControlParams,
                                coord_type::Int64,
                                delta_t::Float64)

Apply control: Remove infectious holding + vaccination (& possibly housing livestock) in response to notification of infection within a specified distance.

Inputs:
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `n_holdings::Int64`: Number of holding in landscape.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `epi_params::EpiParams`: Composite type containing epimidological parameters..
- `control_params::ControlParams`: Variables related to implementing control measures
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function control_housing_and_vacc_fn!(holding_info_vec::Array{HoldingInfo,1},
                                        n_holdings::Int64,
                                        holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                        epi_params::EpiParams,
                                        control_params::ControlParams,
                                        coord_type::Int64,
                                        delta_t::Float64)

    #---------------------------------------------------------------------------
    # FIND AND UPDATE HOLDING TO BE CULLED
    # FIND AND UPDATE TIME TO VACCINE BECOMING EFFECTIVE
    #---------------------------------------------------------------------------
    cull_holding_and_update_time_to_all_livestock_vacc_effective_fn!(holding_info_vec,
                                                                    n_holdings,
                                                                    holding_vectors_and_arrays_params,
                                                                    epi_params.removal_time,
                                                                    control_params.vacc_efficacy,
                                                                    delta_t)

    #---------------------------------------------------------------------------
    # EVALUATE IF VACCINATION OCCURS
    # FOR HOLDINGS NOT PREVIOUSLY VACCINATED, HAS A NOTIFIED INFECTION OCCURRED ON
    # CURRENT TIMESTEP WITHIN THE THRESHOLD DISTANCE?
    #---------------------------------------------------------------------------

    # NEED HOLDING NOTIFYING INFECTION ON CURRENT TIMESTEP
    newly_notified_infection_holding_IDs = findall(holding_vectors_and_arrays_params.holding_status.==(epi_params.detection_time+delta_t))
    n_newly_notified_infection_holding_IDs = length(newly_notified_infection_holding_IDs)

    # If there are newly notified infections, enter loop to check whether
    # any holdings vaccinate
    if n_newly_notified_infection_holding_IDs > 0

        # Iterate over each holding
        #   - Check if unvaccinated & vaccinates based on risk
        #   - Check distance to holding now reporting infection
        #   - If within the threshold distance, apply vaccination
        for holding_idx = 1:n_holdings

            # Check holding is unvaccinated,
            # part of a group that vaccinates in response to reported infections within given distance
            # (i.e. holding_info_vec[holding_idx].intervention_threshold_distance > 0)
            # and not had notified infection
            if ( (holding_vectors_and_arrays_params.holding_vacc_status[holding_idx] == 0) &&
                 (holding_info_vec[holding_idx].intervention_threshold_distance > 0) &&
                 (0<=holding_vectors_and_arrays_params.holding_status[holding_idx]<epi_params.detection_time) )

                # Iterate over each holding with newly notified infection
                for newly_notified_holding_idx = 1:n_newly_notified_infection_holding_IDs

                    # Get index of the newly notified holding in this iteration of for loop
                    newly_notified_holding_ID = newly_notified_infection_holding_IDs[newly_notified_holding_idx]

                    # Check distance to holding now reporting infection
                    if coord_type == 1 #Cartesian co-ords (metres)
                        dist_between_holdings =  eucl_distance(holding_vectors_and_arrays_params.holding_loc_X_vals[newly_notified_holding_ID],
                                            holding_vectors_and_arrays_params.holding_loc_Y_vals[newly_notified_holding_ID],
                                            holding_vectors_and_arrays_params.holding_loc_X_vals[holding_idx],
                                            holding_vectors_and_arrays_params.holding_loc_Y_vals[holding_idx])
                    elseif coord_type == 2 #Cartesian co-ords (metres)
                        dist_between_holdings = eucl_distance_convert_to_metres(holding_vectors_and_arrays_params.holding_loc_X_vals[newly_notified_holding_ID],
                                                            holding_vectors_and_arrays_params.holding_loc_Y_vals[newly_notified_holding_ID],
                                                            holding_vectors_and_arrays_params.holding_loc_X_vals[holding_idx],
                                                            holding_vectors_and_arrays_params.holding_loc_Y_vals[holding_idx])
                    elseif coord_type == 3 #Lat/Long co-ords
                        dist_between_holdings = great_circle_distance(holding_vectors_and_arrays_params.holding_loc_Y_vals[newly_notified_holding_ID],
                                                                        holding_vectors_and_arrays_params.holding_loc_X_vals[newly_notified_holding_ID],  #lat1, lon1
                                                                        holding_vectors_and_arrays_params.holding_loc_Y_vals[holding_idx],
                                                                        holding_vectors_and_arrays_params.holding_loc_X_vals[holding_idx]) #lat2, lon2
                    end

                    # If within the threshold distance, apply vaccination
                    if dist_between_holdings <= holding_info_vec[holding_idx].intervention_threshold_distance
                        apply_livestock_vacc_and_housing_fn!(holding_info_vec[holding_idx],
                                                            holding_idx,
                                                            n_holdings,
                                                            holding_vectors_and_arrays_params,
                                                            epi_params.detection_time,
                                                            control_params)

                        # Break out the inner for loop.
                        # No need to check other newly notified holding
                        break
                    end
                end
            end
        end
    end

    return nothing
end

"""
    apply_livestock_vacc_and_housing_fn!(holding_data_fields::HoldingInfo,
                                            holding_idx::Int64,
                                            n_holdings::Int64,
                                            holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                            detection_time::Float64,
                                            control_params::ControlParams)

For a single holding, check if intervention usage occurs.
If yes, all livestock are vaccinated.
Then also checks if livestock are housed.

Inputs:
- `holding_data_fields::HoldingInfo`: Fields associated with holding specific data.
- `holding_idx::Int64`: ID of holding to be checked.
- `n_holdings::Int64`: Number of holdings in landscape.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `detection_time::Float64`: Time from infection until holding notifies that infection is present.
- `control_params::ControlParams`: Variables related to implementing control measures

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function apply_livestock_vacc_and_housing_fn!(holding_data_fields::HoldingInfo,
                                                holding_idx::Int64,
                                                n_holdings::Int64,
                                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                detection_time::Float64,
                                                control_params::ControlParams)

    if ( (holding_vectors_and_arrays_params.holding_vacc_status[holding_idx]) == 0 &&
         (0<=holding_vectors_and_arrays_params.holding_status[holding_idx]<detection_time) &&
         (holding_data_fields.intervention_threshold_distance > 0) )

        # Vaccinate holding, update holding_vacc_status
        holding_vectors_and_arrays_params.holding_vacc_status[holding_idx] = 1
        holding_data_fields.vacc_status = true

        # All livestock types vaccinated
        # Update livestock_type_vacc_status_by_holding
        n_livestock_types = size(holding_vectors_and_arrays_params.holding_suscept_by_livestock_type,2)
        for livestock_type_idx = 1:n_livestock_types
            holding_vectors_and_arrays_params.livestock_type_vacc_status_by_holding[holding_idx,livestock_type_idx] = 1
        end

        # Update holding vaccinated during current timestep vector
        holding_vectors_and_arrays_params.vacc_holding_during_current_timestep_vec[holding_idx] = true

        # If holding houses livestock, amend susceptibility and transmissibility
        if holding_data_fields.house_livestock_flag == true

            # By livestock type adjustments
            for col_idx = 1:n_livestock_types
                holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx].*(1-control_params.house_livestock_effectiveness)
                holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx].*(1-control_params.house_livestock_effectiveness)
            end
            # Holding-level adjustments
            holding_vectors_and_arrays_params.holding_suscept[holding_idx] = holding_vectors_and_arrays_params.holding_suscept[holding_idx].*(1-control_params.house_livestock_effectiveness)
            holding_vectors_and_arrays_params.holding_transmiss[holding_idx] = holding_vectors_and_arrays_params.holding_transmiss[holding_idx].*(1-control_params.house_livestock_effectiveness)

            # Update housing livestock status variables
            holding_data_fields.house_livestock_status = true
            holding_vectors_and_arrays_params.livestock_housed_during_current_timestep_vec[holding_idx] = true
        end

        # Vaccine only successful in modifying susceptiblity & transmissiblity immediately IF
        # no delay in vaccine becoming effective post being administered AND
        # holding is susceptible
        if (holding_data_fields.holding_remaining_time_to_vacc_becoming_effective == 0)

            # Check holding is still susceptible.
            # If susceptible, apply vaccination effect
            if (holding_vectors_and_arrays_params.holding_status[holding_idx] == 0)
                # Modify holding-level values.
                # Update vaccinated holding susceptibility and transmissibility for all livestock types
                for col_idx = 1:n_livestock_types
                    holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx].*(1-vacc_efficacy)
                    holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx].*(1-vacc_efficacy)
                end

                holding_vectors_and_arrays_params.holding_suscept[holding_idx] = holding_vectors_and_arrays_params.holding_suscept[holding_idx].*(1-control_params.vacc_efficacy)
                holding_vectors_and_arrays_params.holding_transmiss[holding_idx] = holding_vectors_and_arrays_params.holding_transmiss[holding_idx].*(1-control_params.vacc_efficacy)

                # Update vaccination becomes effective during current timestep vector
                holding_vectors_and_arrays_params.vacc_becomes_effective_during_current_timestep_vec[holding_idx] = 1
            end
        end
    end
    return nothing
end

"""
    cull_holding_and_update_time_to_all_livestock_vacc_effective_fn!(holding_info_vec::Array{HoldingInfo,1},
                                                                    n_holdings::Int64,
                                                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                                    removal_time::Float64,
                                                                    vacc_efficacy::Float64,
                                                                    delta_t::Float64)

Find and update: \n
(i) holding to be culled/ending infectious period; \n
(ii) time to vaccine becoming effective (vaccine administered to all livestock types).

Inputs:
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `n_holdings::Int64`: Number of holdings in landscape.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `removal_time::Float64`: Time from infection until livestock on holding are culled.
- `vacc_efficacy::Float64`: Efficacy of vaccine.
- `delta_t::Float64`: Timestep.

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function cull_holding_and_update_time_to_all_livestock_vacc_effective_fn!(holding_info_vec::Array{HoldingInfo,1},
                                                                        n_holdings::Int64,
                                                                        holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                                        removal_time::Float64,
                                                                        vacc_efficacy::Float64,
                                                                        delta_t::Float64)
    for holding_idx = 1:n_holdings
        # If applicable, find and update holding to be culled
        if holding_vectors_and_arrays_params.holding_status[holding_idx]>removal_time
            #Cull IP/end of IP infectious period
            holding_vectors_and_arrays_params.holding_status[holding_idx] = -1.

            #Update holding culled during current timestep vector
            holding_vectors_and_arrays_params.cull_holding_during_current_timestep_vec[holding_idx] = 1
        end

        # If applicable, find and update time to vaccine becoming effective
        if holding_vectors_and_arrays_params.holding_vacc_status[holding_idx] == 1 # Find holding that are vaccinated

            # Check if time to inoculation is above zero.
            if holding_info_vec[holding_idx].holding_remaining_time_to_vacc_becoming_effective > 0
                # If satisfied, decrease by timestep
                holding_info_vec[holding_idx].holding_remaining_time_to_vacc_becoming_effective -= delta_t

                # Error check
                # holding_time_to_inoculation[holding_idx] should not decrease below zero
                if holding_info_vec[holding_idx].holding_remaining_time_to_vacc_becoming_effective < 0
                    error("holding_time_to_inoculation[$holding_idx]: $(holding_time_to_inoculation[holding_idx]). Invalid.")
                end

                # If time to inoculation reaches zero, apply effect of vaccination (if valid)
                if holding_info_vec[holding_idx].holding_remaining_time_to_vacc_becoming_effective == 0

                    # Modify holding-level values for susceptibility and transmissibility for each livestock type.
                    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
                    activate_livestock_vacc_fn!(holding_idx,
                                                holding_vectors_and_arrays_params,
                                                vacc_efficacy)
                end
            end
        end
    end
end

"""
    activate_livestock_vacc_fn!(holding_idx::Int64,
                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                    vacc_efficacy::Float64)

For a single holding, apply effect of livestock vaccination strategy.

Inputs:
- `holding_idx::Int64`: ID of holding to be checked.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `vacc_efficacy::Float64`: Efficacy of vaccine.

Outputs: None \n
Location: supporting\\_fns\\_interventions.jl
"""
function activate_livestock_vacc_fn!(holding_idx::Int64,
                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                    vacc_efficacy::Float64)

    # Vaccine only successful in modifying susceptiblity & transmissiblity if susceptible
    if holding_vectors_and_arrays_params.holding_status[holding_idx] == 0

        # Update vaccinated holding susceptibility and transmissibility for each livestock type.
        n_livestock_types = size(holding_vectors_and_arrays_params.holding_suscept_by_livestock_type,2)
        for col_idx = 1:n_livestock_types
            holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_suscept_by_livestock_type[holding_idx,col_idx].*(1-vacc_efficacy)
            holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx] = holding_vectors_and_arrays_params.holding_transmiss_by_livestock_type[holding_idx,col_idx].*(1-vacc_efficacy)
        end

        # Updating holding-level susceptibility and transmissibility values
        holding_vectors_and_arrays_params.holding_suscept[holding_idx] = holding_vectors_and_arrays_params.holding_suscept[holding_idx].*(1-vacc_efficacy)
        holding_vectors_and_arrays_params.holding_transmiss[holding_idx] = holding_vectors_and_arrays_params.holding_transmiss[holding_idx].*(1-vacc_efficacy)

        # Update vaccination due to become effective during current timestep vector
        holding_vectors_and_arrays_params.vacc_becomes_effective_during_current_timestep_vec[holding_idx] = 1
    end
    return nothing
end
