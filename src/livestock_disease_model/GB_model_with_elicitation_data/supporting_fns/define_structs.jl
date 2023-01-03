#=
Purpose:
Define composite types (structured data types)
=#

"""
    Structure: LandscapeData

Contains spatial landscape information.

Location: define\\_structs.jl
"""
@with_kw struct LandscapeData
    get_n_nodes::Int64  #Number of nodes in the landscaspe
    longest_side::Float64  #Length of longest dimension of box containing all premises
    get_trans::Float64
    get_susc::Float64
    coord_type::Int64 #Type of coordinate system in use
end

"""
    Structure: HoldingInfo

Fields associated with holding-level data.

Location: define\\_structs.jl
"""
@with_kw mutable struct HoldingInfo
    const holding_ID::Int64 # Identifier for the holding
    const country_ID::String # Identifier for country
    const county_ID::Int64  # Numeric Identifier for county
    const county_name::String  # Name of county the holding resides in
    const X_loc::Float64 # Easting location info
    const Y_loc::Float64 # Northing location info
    const livestock_count::Vector{Int64} # Livestock at holding. Entry per livestock type.
    behavioural_group_ID::Int64 = 0 # Determines distribution of intevention trigger timings to sample from
    intervention_timing_group_ID::Int64 = 0 # Determines the intervention_threshold_distance
    intervention_threshold_distance::Float64 = -1. # Once there is a reported holding within this distance, implement interventions at this holding
    linked_holdings_for_control::Vector{Int64} = Vector{Int64}(undef,0) # IDs of premises to undergo control if index holding is infected
    house_livestock_flag::Bool = false  # Specifies if holding will house livestock when intervention use is triggered
    house_livestock_status::Bool = false # True if livestock are housed, false if not
    vacc_status::Bool = false # True if vaccinated, false if not
    holding_vacc_inoclulation_time::Float64 = -1.
    holding_remaining_time_to_vacc_becoming_effective::Float64 = -1.
end

"""
    Structure: TimeParams

Contains simulation time associated fields.

Location: define\\_structs.jl
"""
@with_kw struct TimeParams
    timestep_val::Float64 = 1.  # Time increment per iteration of the model
    max_time::Float64 = 10*365. # Length of the time horizon for the simulation
end

"""
    Structure: EpiParams

Composite type containing epimidological parameters.

Location: define\\_structs.jl
"""
@with_kw struct EpiParams
    incubation_time::Float64 = 5.
    detection_time::Float64 = 9.
    removal_time::Float64 = 13.
    per_livestock_type_suscep::Vector{Float64} = [1.] # Susceptiblity scale parameters for each livestock type
    per_livestock_type_transmiss::Vector{Float64} = [1.] # Transmissibility scale parameters for each livestock type
    suscept_exponent::Vector{Float64} = [1.] # Susceptibility exponents for each livestock type
    transmiss_exponent::Vector{Float64} = [1.] # Transmissibility exponents for each livestock type
end

"""
    Structure: ControlParams

Composite type containing control/intervention associated parameters.

Location: define\\_structs.jl
"""
@with_kw struct ControlParams
    run_controls_fn::Function = remove_IPs_only_control_fn! # Specify function to be called to enact specified controls
    vacc_efficacy::Float64 = 1.  # As a proportion: 0 = no effect, 1 = maximum protection.
    holding_time_to_inoculation_fn::Function = five_day_inoculation_fn! # Set up delay in vaccine becoming effective.
    house_livestock_effectiveness::Float64 = 0.5 # Housing livestock effectiveness. (0 = no effect, 1 = maximum protection)
end


"""
    Structure: GridData

Fields associated with grid overlay on the landscape.

Location: define\\_structs.jl
"""
@with_kw struct GridData
    coord_type::Int64  #Type of coordinate system in use
    bounding_box_vals::Vector{Float64} # Limits to the grid the landscape is contained within
    grid_optim_method::Int64 = 2 # Specify whether grid optimisation method should be used. (1 = fixed, 2 = dynamic)
    nMax::Int64 = 0 # For adaptive grid construction method, number of nodes in one cell when using dynamic cell size.
    n_cells_per_side::Int64 = 0 # For static grid construction method, number of cells along one dimension of the grid (total number will be this^2).
end

"""
    Structure: FunctionCollection

Composite type containing collection of functions used in model simulation.

Location: define\\_structs.jl
"""
@with_kw struct FunctionCollection
    calc_holding_suscep_transmiss_fn::Function # Function to calculate holding-level susceptibility and transmissibility
    kernel_fn::Function # Defines risk of transmission w.r.t. distance
    iterate_outbreak_fn::Function # Function to perform disease transitions and control implementation per timestep
    run_simn_replicate_fn::Function # Function to perform single outbreak replicate
end

"""
    Structure: HoldingVectorsAndArrays

Set of vectors and arrays associated with holding attributes.

In each vector, each entry corresponds to a distinct holding.

Location: define\\_structs.jl
"""
@with_kw struct HoldingVectorsAndArrays
    holding_cell_IDs::Vector{Int64} = Vector{Int64}(undef,0) # Entry per holding. Cell ID the holding resides in.
    holding_loc_X_vals::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. East-West axis location.
    holding_loc_Y_vals::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. North-South axis location.
    holding_loc_X_and_Y_vals::Array{Float64,2} = Array{Float64,2}(undef,0,0) # Row per holding. Columns: [East-West axis location, North-South axis location]
    holding_livestock_data::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Number of livestock.
    holding_suscept::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Overall susceptibility at the holding.
    holding_transmiss::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Overall transmissibility at the holding.
    holding_suscept_by_livestock_type::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Susceptibility at the holding per livestock type.
    holding_transmiss_by_livestock_type::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Transmissibility at the holding per livestock type.
    holding_status::Vector{Float64} = Vector{Float64}(undef,0) # Entry per holding. Disease status.
    holding_has_had_infection_flag::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding has had infection, false otherwise.
    holding_vacc_status::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding has had livestock vaccinated, false otherwise.
    livestock_type_vacc_status_by_holding::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding has had livestock vaccinated, false otherwise
    holding_infectious_flag_vec::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding is currently infectious, false otherwise.
    cull_holding_during_current_timestep_vec::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding no longer infectious on current timestep, false otherwise.
    vacc_holding_during_current_timestep_vec::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if holding vaccinated during current timestep, false otherwise.
    vacc_becomes_effective_during_current_timestep_vec::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if vaccine becomes effective at holding during current timestep, false otherwise.
    livestock_housed_during_current_timestep_vec::Vector{Bool} = Vector{Bool}(undef,0) # Indicator vector, entry per holding. True if livestock housed at holding during current timestep, false otherwise.
end
