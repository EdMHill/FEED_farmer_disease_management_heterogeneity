#=
Purpose:
File to house supporting functions for running GB model,
associated with transmission dynamics

Contains:
 - compute_holding_suscep_transmiss_cattle_only
 - power_law_transmission_kernel (& construct_power_law_transmission_kernel)
 - seed_infection (& supporting funcction "order_holdings_by_distance_to_index_infected_holding")
 - compute_local_pairwise_prob (infection event probability between two holdings)
 - conditional_subsample_alg! (run grid simulation using conditional subsample algorithm)
 - local_pairwise_check_sucep_infectious_paris_in_cell! (used with conditional subsample algorithm, run infection event check for holdings within same cell)
=#

"""
    compute_holding_suscep_transmiss_cattle_only(holding_livestock_data::Array{Int64,1},
                                                 epi_params::EpiParams)

Calculate premises-level susceptibility and transmissibility for pathogen.

Inputs:
- ``holding_livestock_data::Array{Int64,1}`: Number of cattle per premises. Row per holding.
- `epi_params::EpiParams`: Composite type containing epimidological parameters.

Outputs:
- `holding_suscept::Vector{Float64}`, `holding_transmiss::Vector{Float64}`: Premises-level susceptibility and transmissibility.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function compute_holding_suscep_transmiss_cattle_only(holding_livestock_data::Array{Int64,1},
                                                        epi_params::EpiParams)

    # Use holding_livestock_data input data,(Number of each livestock type per premises. Row per premises, column per animal)
    # with per_cow_suscept, per_cow_suscept, suscept_exponent, transmiss_exponent

    # For each holding, get contribution to susceptibility and transmissibility from each livestock type
    holding_suscept = epi_params.per_livestock_type_suscep[1].*(holding_livestock_data.^epi_params.suscept_exponent[1])
    holding_transmiss = epi_params.per_livestock_type_transmiss[1].*(holding_livestock_data.^epi_params.transmiss_exponent[1])
        # Calculation outline:
        # Exponent i applied to col_i of holding_livestock_data
        # Multiply col_i of per_cow_suscept/per_cow_transmiss by col_i of holding_livestock_data
        # Sum across columns to get overall premises-level value

    # Populate outputs for livestock type variables
    holding_suscept_by_livestock_type = 1. *holding_suscept
    holding_transmiss_by_livestock_type = 1. *holding_transmiss

    return  holding_suscept::Vector{Float64},
            holding_transmiss::Vector{Float64},
            holding_suscept_by_livestock_type::Array{Float64,1},
            holding_transmiss_by_livestock_type::Array{Float64,1}
end


#===============================================================================
TRANSMISSION KERNEL FNS
===============================================================================#
"""
    construct_power_law_transmission_kernel(MaxDist::Int64)

Create lookup vector for power law transmission kernel.

Inputs:
- `MaxDist::Int64`: Maximum distance to compute value of transmission kernel for.

Outputs:
- `kernel_lookup_vec::Array{Float64,1}`: Kernel term, with entry for each one metre increment.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function construct_power_law_transmission_kernel(MaxDist::Int64)

    #Initialise lookup array
    kernel_lookup_vec = zeros(Float64,MaxDist)

    #Iterate over desired distances. Assign kernel value to array
    for DistItr = 1:MaxDist
        kernel_lookup_vec[DistItr] = power_law_transmission_kernel(convert(Float64,DistItr))
    end

    return kernel_lookup_vec::Array{Float64,1}
end

"""
    power_law_transmission_kernel(d::Float64)

Calculate value of power law transmission kernel.

Inputs:
- `d::Float64`: Distance to compute value of transmission kernel for.

Outputs:
- `kernel_val::Float64`: Value of transmission kernel at distance d.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function power_law_transmission_kernel(d::Float64)

    #List of required parameter values
    # k1::Float64 = 0.08912676
    # k2::Float64 = 1600.0
    # k3::Float64 = 4.6
    k2::Float64 = 2000.
    k3::Float64 = 2.
        # non_normalised_integral: 8.091926327895512e7

    # Set limit on transmission distance to 50km (50,000m)
    max_transmission_distance = 50000

    # Set kernel value based on distance d
    # If above max transmission distance, set to zero
    if d > max_transmission_distance
        kernel_val = 0.
    else
        # Get normalisation constant, based on max_transmission_distance, k2 & k3 values
        non_normalised_integral, err = quadgk(d -> ((2*pi*d) / (1 + ((d/k2)^k3))), 0, max_transmission_distance, rtol=1e-10)

        # Set normalisation constant
        k1 = 1/non_normalised_integral

        # Express transmission kernel formulation
        kernel_val = k1 / (1 + ((d/k2)^k3))
    end

    # Return kernal value for distance d from function
    return kernel_val::Float64
end

#===============================================================================
SEED INFECTION FNS
===============================================================================#

"""
    seed_infection(initial_infection_info::Array{Any,1},
                    holding_info_vec::Array{HoldingInfo,1},
                    holding_loc_vals::Array{Float64,2},
                    coord_type::Int64,
                    rng::AbstractRNG)

Seed initial infection cases according to a specified method.

Inputs:
- `initial_infection_info`: (tuple) [seed_method,NumOfNodes/NodeIDs]
                             Seed method.
                              - 1 = random,
                              - 2/3 = single/group specific node id(s),
                              - 4 = from file (one id per row, will seed all the node ids given in the file each replicate.).
                              - 5 = seed a random site and it's N nearest neighbours.
                              - 6 = seed a random site and 2 nearest neighbours
                             Number of nodes to seed each replicate if seed_method = 1 or 5; or if seed_method = 2/3, seed this specific node every replicate.
                                seed_method = 4 from file, node ids to seed given by one id/row, number of lines must be == number of replicates.
                                seed method = 6, list of IDs corresponds to counties/unitary authorities where index infected node may be selected
- `holding_info_vec::Array{HoldingInfo,1}`: Structure with fields associated with holding specific data.
- `holding_loc_vals::Array{Float64,2}`: Coordinates for each holding.
- `coord_type::Int64`: Value 1 ("Cartesian", metres) or 2 ("Cartesian", km).
- `rng::AbstractRNG`: The random number generator.

Outputs:
- `seed_holding_IDs::Array{Int64,1}`: Vector of 1's (infected) and 0's (not infected). Entry per node that has infection.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function seed_infection(initial_infection_info::Array{Any,1},
                        holding_info_vec::Array{HoldingInfo,1},
                        holding_loc_vals::Array{Float64,2},
                        coord_type::Int64,
                        rng::AbstractRNG)

    #Disaggregate initial_infection_info
    seed_method = initial_infection_info[1]::Int64

    # Get number of holdings
    n_holdings::Int64 = length(holding_info_vec)

    if seed_method == 1 #Random
        #Get number of premises to be seeded with infection each replicate
        n_holdings_seed_infections = initial_infection_info[2][1]::Int64

        # Draw IDs based on number of premises in use.
        # Use rng to get consistency in selected premises IDs across control runs
        seed_holding_IDs = rand(rng,1:n_holdings,n_holdings_seed_infections)
    elseif seed_method == 2 #Single specific node ID
        seed_holding_IDs = [initial_infection_info[2]]::Array{Int64,1}
    elseif seed_method == 3 #Set of specific node ID
        seed_holding_IDs = initial_infection_info[2]::Array{Int64,1}
    elseif seed_method == 4 #From file (filename stored in initial_infection_info[2])
        #One id per row, will seed all the node ids given in the file each replicate.
        seed_holding_IDs = readdlm(initial_infection_info[2], Int64)
    elseif seed_method == 5
        # Single random node and it's (n_holdings_seed_infections-1) nearest neighbours

        # Disaggregate location data
        holding_loc_X_val::Array{Float64,1} = holding_loc_vals[:,1]
        holding_loc_Y_val::Array{Float64,1} = holding_loc_vals[:,2]

        #Get number of premises to be seeded with infection each replicate
        n_holdings_seed_infections = initial_infection_info[2][1]::Int64

        # Draw the index case from which cluster will be generated
        # Use rng to get consistency in selected premises IDs across control runs
        index_seed_holding_ID = rand(rng,1:n_holdings,1)[1]
        println("index_seed_holding_ID: $index_seed_holding_ID")

        # Get holding IDs in ascending order (by distance) to index infected holding
        ordered_distance_holding_IDs = order_holdings_by_distance_to_index_infected_holding(n_holdings,
                                                                                            index_seed_holding_ID,
                                                                                            coord_type,
                                                                                            holding_loc_X_val,
                                                                                            holding_loc_Y_val)

        # Retain IDs of index case & (n_holdings_seed_infections-1) nearest neighbours
        seed_holding_IDs = [index_seed_holding_ID;ordered_distance_holding_IDs[1:(n_holdings_seed_infections-1)]]

        # Error check
        # Should be no duplicates in seed_holding_IDs
        if length(unique(seed_holding_IDs)) < n_holdings_seed_infections
            error("length(unique(seed_holding_IDs)): $(length(unique(seed_holding_IDs)))); Less than n_holdings_seed_infections ($n_holdings_seed_infections). Invalid.")
        end
    elseif seed_method == 6
        # Single random node within specified regions
        # & N nearest neighbours

        # Set the total number of seed infections
        n_holdings_seed_infections = 3

        # Disaggregate location data
        holding_loc_X_val = holding_loc_vals[:,1]
        holding_loc_Y_val = holding_loc_vals[:,2]

        # Get IDs of regions where index infection may be selected from
        seed_infections_region_IDs = initial_infection_info[2]::Vector{Int64}

        # Extract holding ID and county ID from holding paramter structure
        # Then select index infection from eligible regions
        holding_IDs_vec = getfield.(holding_info_vec,:holding_ID)
        holding_county_IDs_vec = getfield.(holding_info_vec,:county_ID)
        eligible_holding_IDs = holding_IDs_vec[holding_county_IDs_vec .∈ [seed_infections_region_IDs]]

        # Draw the index case from which cluster will be generated
        # Use rng to get consistency in selected premises IDs across control runs
        index_seed_holding_ID = rand(rng,eligible_holding_IDs,1)[1]
        println("index_seed_holding_ID: $index_seed_holding_ID")

        # Get holding IDs in ascending order (by distance) to index infected holding
        ordered_distance_holding_IDs = order_holdings_by_distance_to_index_infected_holding(n_holdings,
                                                                                            index_seed_holding_ID,
                                                                                            coord_type,
                                                                                            holding_loc_X_val,
                                                                                            holding_loc_Y_val)

        # Retain IDs of index case & (n_holdings_seed_infections-1) nearest neighbours
        seed_holding_IDs = [index_seed_holding_ID;ordered_distance_holding_IDs[1:(n_holdings_seed_infections-1)]]
    else #Invalid value, throw error
        error("seed_method has value $seed_method. seed_method should take value 1, 2, 3, 4 or 5.")
    end

    return seed_holding_IDs::Array{Int64,1}
end

"""
        order_holdings_by_distance_to_index_infected_holding(n_holdings::Int64,
                                                            index_seed_holding_ID::Int64,
                                                            coord_type::Int64,
                                                            holding_loc_X_val::Vector{Float64},
                                                            holding_loc_Y_val::Vector{Float64})

Get holding IDs in ascending order (by distance) to index infected holding.

Inputs:
- `n_holdings::Int64`: Number of cattle per premises. Row per holding.
- `index_seed_holding_ID::Int64`: ID of holding selected as initial infected location.
- `coord_type::Int64`: Value 1 ("Cartesian", metres) or 2 ("Cartesian", km).
- `holding_loc_X_val::Vector{Float64}`: East-West plane location data for each holding.
- `holding_loc_Y_val::Vector{Float64}`: North-South plane location data for each holding.

Outputs:
- `ordered_distance_holding_IDs::Vector{Int64}`: IDs of holdings in sorted ascending distance order (relative to index infected holding).

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function order_holdings_by_distance_to_index_infected_holding(n_holdings::Int64,
                                                                index_seed_holding_ID::Int64,
                                                                coord_type::Int64,
                                                                holding_loc_X_val::Vector{Float64},
                                                                holding_loc_Y_val::Vector{Float64})

    # Get distances from index case to all other premises
    # Check distance to premises now reporting infection
    dist_to_index_prem = zeros(Float64,n_holdings) # Initialise vector to store each premises to index premises distance
    for holding_ID = 1:n_holdings
        if coord_type == 1 #Cartesian co-ords (metres)
            dist_to_index_prem[holding_ID] = eucl_distance(holding_loc_X_val[index_seed_holding_ID],
                                holding_loc_Y_val[index_seed_holding_ID],
                                holding_loc_X_val[holding_ID],
                                holding_loc_Y_val[holding_ID])
        elseif coord_type == 2 #Cartesian co-ords (km)
            dist_to_index_prem[holding_ID] = eucl_distance_convert_to_metres(holding_loc_X_val[index_seed_holding_ID],
                                                holding_loc_Y_val[index_seed_holding_ID],
                                                holding_loc_X_val[holding_ID],
                                                holding_loc_Y_val[holding_ID])
        end
    end

    # Get indexes of holdings in sorted ascending distance order
    ordered_distance_holding_IDs = sortperm(dist_to_index_prem)

    # Remove the index case ID
    filter!(x->x≠index_seed_holding_ID,ordered_distance_holding_IDs)

    return ordered_distance_holding_IDs::Vector{Int64}

end

#===============================================================================
INFECTION EVENT CHECK FNS (used with Sellke construction simn)
===============================================================================#

"""
    compute_local_pairwise_prob(suscep_holding_X_loc::Float64,
                                suscep_holding_Y_loc::Float64,
                                InfNode_xLoc::Float64,
                                InfNode_yLoc::Float64,
                                holding_transmiss::Vector{Float64},
                                holding_suscept::Vector{Float64},
                                coord_type::Int64,
                                kernel_lookup_vec::Array{Float64,1},
                                delta_t::Float64,
                                selected_infected_holding_ID::Int64,
                                selected_suscep_holding_ID::Int64)

Calculates directly for each infectious-susceptible pair in the population the
probability of infection occurring, which is evaluated as a Benoulli process.

Inputs:
- `suscep_holding_X_loc::Array{Float64,1}`: East-west plane co-ordinate of susceptible premises.
- `suscep_holding_Y_loc::Array{Float64,1}`: North-south plane co-ordinate of susceptible premises.
- `InfNode_xLoc::Array{Float64,1}`: East-west plane co-ordinate of infectious premises.
- `InfNode_yLoc::Array{Float64,1}`: North-south plane co-ordinate of infectious premises.
- `holding_transmiss::Vector{Float64}`: Premises-level transmissibility.
- `holding_suscept::Vector{Float64}`: Premises-level susceptibility.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep increment.
- `selected_infected_holding_ID::Int64`: ID of the infectious node.
- `selected_suscep_holding_ID::Int64`: ID of the susceptible node.

Outputs:
- `PairwiseInfProb::Float64`: Probability of infection.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function compute_local_pairwise_prob(suscep_holding_X_loc::Float64,
                                    suscep_holding_Y_loc::Float64,
                                    InfNode_xLoc::Float64,
                                    InfNode_yLoc::Float64,
                                    holding_transmiss::Vector{Float64},
                                    holding_suscept::Vector{Float64},
                                    coord_type::Int64,
                                    kernel_lookup_vec::Array{Float64,1},
                                    delta_t::Float64,
                                    selected_infected_holding_ID::Int64,
                                    selected_suscep_holding_ID::Int64)

    #Calculate distance between the two points
    if coord_type == 1 #Cartesian co-ords (metres)
        d =  eucl_distance(InfNode_xLoc,
                            InfNode_yLoc,
                            suscep_holding_X_loc,
                            suscep_holding_Y_loc)
    elseif coord_type == 2 #Cartesian co-ords (metres)
        d = eucl_distance_convert_to_metres(InfNode_xLoc,
                                            InfNode_yLoc,
                                            suscep_holding_X_loc,
                                            suscep_holding_Y_loc)
    end

    #Calculate rate of infection
    dist_idx = return_dist_idx_for_kernel(d)
    FOI_rate = holding_transmiss[selected_infected_holding_ID]*holding_suscept[selected_suscep_holding_ID]*kernel_lookup_vec[dist_idx]

    #Compute probability of infection (Use 1-exp function!)
    PairwiseInfProb = oneMinusExp(-FOI_rate*delta_t)

    return PairwiseInfProb::Float64
end


#===============================================================================
ALGORITHM: Cumulative pairwise infection probability, against single suscetpible unit
===============================================================================#
"""
    compute_pairwise_cumulative_infection_prob(holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                        coord_type::Int64,
                                        kernel_lookup_vec::Array{Float64,1},
                                        delta_t::Float64,
                                        infectious_holding_in_transmission_cell_IDs::Array{Int64,1},
                                        selected_suscep_holding_ID::Int64)

Compute cumulative pairwise infection probability against a single susceptible unit.

Inputs:
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep increment.
- `infectious_holding_in_transmission_cell_IDs::Int64`: IDs of units that are infectious within current cell.
- `selected_suscep_holding_ID::Int64`: ID of the susceptible node.

Outputs:
- `cumulative_p::Float64`: Probability that SusNode is infected.

Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function compute_pairwise_cumulative_infection_prob(holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                    coord_type::Int64,
                                                    kernel_lookup_vec::Array{Float64,1},
                                                    delta_t::Float64,
                                                    infectious_holding_in_transmission_cell_IDs::Array{Int64,1},
                                                    selected_suscep_holding_ID::Int64)

    #Get location of susceptible node
    suscep_holding_X_loc = holding_vectors_and_arrays_params.holding_loc_X_vals[selected_suscep_holding_ID] #Susceptible node x location
    suscep_holding_Y_loc = holding_vectors_and_arrays_params.holding_loc_Y_vals[selected_suscep_holding_ID]

    #println("suscep_holding_X_loc, suscep_holding_Y_loc: $suscep_holding_X_loc, $suscep_holding_Y_loc")

    #Initialise no infection probability variable
    no_infection_prob = 1.0

    #Iterate over each infectious node. Get probability of susceptible not being infected by ANY infectious node
    for infectious_holding_idx = 1:length(infectious_holding_in_transmission_cell_IDs)

        #Current infectious node of interest
        selected_infectious_holding_ID = infectious_holding_in_transmission_cell_IDs[infectious_holding_idx]

        #Location of infectious unit
        InfNode_xLoc = holding_vectors_and_arrays_params.holding_loc_X_vals[selected_infectious_holding_ID]
        InfNode_yLoc = holding_vectors_and_arrays_params.holding_loc_Y_vals[selected_infectious_holding_ID]

        #Calculate distance between the two points
        if coord_type == 1 #Cartesian co-ords (metres)
            d =  eucl_distance(InfNode_xLoc,
                                InfNode_yLoc,
                                suscep_holding_X_loc,
                                suscep_holding_Y_loc)
        elseif coord_type == 2 #Cartesian co-ords (metres)
            d = eucl_distance_convert_to_metres(InfNode_xLoc,
                                            InfNode_yLoc,
                                            suscep_holding_X_loc,
                                            suscep_holding_Y_loc)
        elseif coord_type == 3 #Lat/Long co-ords
            d = great_circle_distance(InfNode_yLoc, InfNode_xLoc,  #lat1, lon1
                                    suscep_holding_Y_loc, suscep_holding_X_loc) #lat2, lon2
        end


        #Calculate rate of infection
        dist_idx = return_dist_idx_for_kernel(d)
        FOI_rate = holding_vectors_and_arrays_params.holding_transmiss[selected_infectious_holding_ID]*holding_vectors_and_arrays_params.holding_suscept[selected_suscep_holding_ID]*kernel_lookup_vec[dist_idx]

        #Compute probability of infection (Use 1-exp function!)
        local_prob_val = oneMinusExp(-FOI_rate*delta_t)

        #println("holding_transmiss, holding_suscept, KernalVal: $(holding_transmiss[selected_infectious_holding_ID]), $(holding_suscept[selected_suscep_holding_ID]), $(kernel_lookup_vec[dist_idx])")
        #println("d, FOI_rate, local_prob_val: $d, $FOI_rate, $local_prob_val")

        #Revise no infection probability
        no_infection_prob = no_infection_prob * (1.0 -  local_prob_val)  #p_aj = p_aj * (1 - p_ij)
    end

    # One susceptible node now analysed.
    # Update cumulative probability of at least one node in
    # infectious cell infecting this susceptible holding
    # Equivalent to (1 - prob. no transmission)
    cumulative_p = 1.0 - no_infection_prob

    return cumulative_p::Float64
end

#===============================================================================
ALGORITHM: Conditional subsample algorithm
===============================================================================#
#Utilise gridding approach to partition the spatial landscape into grids.

#-------------------------------------------------------------------------------
### Parameter definitions
#-------------------------------------------------------------------------------

# p_ij: Exact probability that node i infects node j
# p_ij = 1 - exp(-S_j T_i K(d_ij) deltat)

# v_ab: Upper boundary for the probability of one infectious node i in cell a infecting one susceptible
#       node j in cell b based on over-estimated transmission parameters.
# v_ab = 1 - exp(-HatT_a * HatS_b * K(d_ab))

# HatT_a: Maximum transmissibility of any node in cell a
# HatS_b: Maximum susceptibility of any node in cell b

# w_ab: The probability at least one of the infectious nodes in cell a infects one node in b given max suscep
# and transmissibility
# w_ab = 1 - (1 - v_ab)^|I_a|

#-------------------------------------------------------------------------------
### Steps in words
#-------------------------------------------------------------------------------

# Step 1: For each cell a containing at least one infectious node, iterate over all other
#         cells b and generate a random number nb from a binomial distribution with
#         number of trials |Jb| and probability ωab.

# Step 2: Randomly select a subset of n_b nodes from b, each with equal probability of being selected

# Step 3: For each of the nodes j in the selected subset of nodes, calculate the
#         cumulative true probability of the event that j is infected by one or more of the infectious nodes in a.

# Step 4:  Evaluate if infection actually takes place using the
#         true probability conditional on the upper boundary of the probability used in step 1

#-------------------------------------------------------------------------------
### Psuedocode (also includes case where Pairwise Algorithm is needed.)
#-------------------------------------------------------------------------------

#for each cell with infectious nodes a
    #for each cell with susceptible nodes b
        #if a is not b
            #N_b = set of susceptibles in b
            #|J_b| = size of N_b
            #n_b = binomial random variate(|J_b|, omega_ab)
            #if (n_b > 0)
                #N_a = set of infectious in a
                #N_sample = random subsample of size n_b from N_b
                #for each j in N_sample
                    #p_aj = 1.0
                    #for each i in N_a
                        #p_aj = p_aj * (1 - p_ij)
                    #end
                    #R ~ U(0,1)
                    #if R < ()(1-p_aj)/omega_ab)
                        #N_sample[j] is infected
                    #end
                #end
            #end
        #else
            #Pairwise algorithm
        #end
    #end
#end


#-------------------------------------------------------------------------------
### C++ code
#-------------------------------------------------------------------------------

# //Conditional subsample method used in main paper.
# unsigned long long int Local_spread::localBinomial_CtoC(Cell* inf_cell, Cell* recipient_cell)
# {
#     unsigned long long int n_kernel_calls = 1;
#     size_t n_suscept = recipient_cell->get_n_suscept();
#     if(n_suscept == 0)
#     {
#         return 0;
#     }
#     //Probability of the most infectious node in c1 infecting the most susc. node in c2 at the shortest poss. distance.
#     double p_over = inf_cell->get_p_over(recipient_cell);
#     //Probability of one farm becomes infected by any of the infectious nodes in c1 using p_over
#     size_t n_inf = inf_cell->get_n_inf();
#     double cumulative_p_over = 1.0 - std::pow(1.0-p_over, n_inf);
#     //Draw the number of nodes that would get infected using the cumulative p_over.
#     size_t n_over = 0;
#     cumulative_p_over = this->generate_bin_rv_lookup(n_suscept, cumulative_p_over, n_over);
# //    size_t n_over = common::draw_binom(n_suscept, cumulative_p_over);
#     if(n_over > 0)
#     {
#         std::vector<Node*> sus_nodes = recipient_cell->get_nodes(0);
#         std::vector<Node*> random_sample;
#         common::random_unique(sus_nodes, n_over, random_sample);
#         std::vector<Node*> inf_nodes = inf_cell->get_nodes(2);
#         for(Node* sus_node : random_sample)
#         {
#             double cumulative_p = 1.0;
#             for(Node* inf_node : inf_nodes)
#             {
#                 cumulative_p *= 1.0 - (this->*p_function)(inf_node, sus_node);
#                 n_kernel_calls++;
#             }
#             cumulative_p = 1.0 - cumulative_p; //Cumulative probability of at least one node in inf_cell infecting this sus_node.
#             if(common::unif_rand() <= cumulative_p / cumulative_p_over)
#             {
#                 sus_node->become_exposed();
#             }
#         }
#     }
#     return n_kernel_calls;
# }

#-------------------------------------------------------------------------------
### Our version
#-------------------------------------------------------------------------------
"""
    conditional_subsample_alg!(rng::AbstractRNG,
                                n_suscept::Int64,
                                p_over::Float64,
                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                suscept_holding_in_recip_cell_IDs::Array{Int64,1},
                                infectious_holding_in_transmitting_cell_IDs::Array{Int64,1},
                                binomial_RNG_array::Array{Binomial{Float64},2},
                                P_CS::Array{Float64,1},
                                coord_type::Int64,
                                kernel_lookup_vec::Array{Float64,1},
                                delta_t::Float64)

Run the conditional subsample algorithm for testing for infection events between epidemiological units in different grid cells.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `n_suscept::Int64`: Number of susceptibles in recipient cell.
- `p_over::Float64`: Probability of the most infectious node in InfCell infecting the most susceptible node in RecipCell at the shortest possible distance.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `suscept_holding_in_recip_cell_IDs::Array{Int64,1}`: For recipient cell, list of nodes that are susceptible.
- `infectious_holding_in_recip_cell_IDs::Array{Int64,1}`: For recipient cell, list of nodes that are infectious.
- `binomial_RNG_array::Array{Binomial{Float64},2}`: Array of Binomial RNGs. Row per susceptible number, column per probability threshold.
- `P_CS::Array{Float64,1})`: Preset probabilities to initialise RNG with.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: supporting\\_fns\\_transmission\\_dynamics.jl
"""
function conditional_subsample_alg!(rng::AbstractRNG,
                                    n_suscept::Int64,
                                    p_over::Float64,
                                    holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                    suscept_holding_in_recip_cell_IDs::Array{Int64,1},
                                    infectious_holding_in_transmitting_cell_IDs::Array{Int64,1},
                                    binomial_RNG_array::Array{Binomial{Float64},2},
                                    P_CS::Array{Float64,1},
                                    coord_type::Int64,
                                    kernel_lookup_vec::Array{Float64,1},
                                    delta_t::Float64)

    #Probability one farm becomes infected by any of the infectious nodes in c1 using p_over
    n_inf::Int64 = length(infectious_holding_in_transmitting_cell_IDs) #Get number of infected in cell with infectious nodes
    cumulative_p_over::Float64 = 1.0 - ((1.0 - p_over)^n_inf)

    #Work around the computationally costly construction of binomial random generators each time
    #Instead, use a previously set up number of generators for all combinations of range of fixed probabilities
    #P_CS = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 5.0e-9}
    #Note: SET UP P_CS IN ASCENDING ORDER!
    binomial_RNG_array_col_p_over_idx = searchsortedlast(P_CS,cumulative_p_over)

    # If probability is 1.0, binomial_RNG_array_col_p_over_idx value is one larger than number of columns in binomial_RNG_array
    # Decrement binomial_RNG_array_col_p_over_idx value so will access the final column of binomial_RNG_array (that holds binomial RNG with probability of success per trial of 1)
    if binomial_RNG_array_col_p_over_idx == length(P_CS)
        binomial_RNG_array_col_p_over_idx = binomial_RNG_array_col_p_over_idx - 1
    end

    # Get index to relevant overestimated probability value from P_CS vector
    # Required index corresponds to first '0' index in searchsortedlast(P_CS,cumulative_p_over)
    #      - Can be obtained by taking next entry after binomial_RNG_array_col_p_over_idx (i.e. binomial_RNG_array_col_p_over_idx + 1)
    rounded_P_CS_idx = binomial_RNG_array_col_p_over_idx + 1

    # Assign overestimated probability that one farm becomes infected by any of the infectious nodes in index infected cell
    overestimate_cumulative_p_over = P_CS[rounded_P_CS_idx]


    #Draw the number of nodes that would get infected using the cumulative p_over.
    n_over::Int64 = rand(rng,binomial_RNG_array[n_suscept,binomial_RNG_array_col_p_over_idx])

    #If nodes would have been infected using p_over, check if infected using actual transmission rate
    if (n_over > 0)

        #Sample n_over nodes from susceptible nodes
        sample_suscep_nodes_r::Array{Int64,1} = StatsBase.sample(1:n_suscept,n_over; replace=false)

        #Iterate over each sample susceptible node.
        #Test if infected by collection of infectious node in InfCell
        for sus_node_idx = 1:n_over
            #Get index of node selected from pool of susceptibles
            sampled_sus_node_idx = sample_suscep_nodes_r[sus_node_idx]

            #Get population level ID of the sampled susceptible node
            selected_sus_node_ID = suscept_holding_in_recip_cell_IDs[sampled_sus_node_idx]

            #Calculate cumulative_p (Cumulative probability of at least one node in inf_cell infecting this sus_node.)
            cumulative_p = compute_pairwise_cumulative_infection_prob(holding_vectors_and_arrays_params,
                                                                        coord_type,
                                                                        kernel_lookup_vec,
                                                                        delta_t,
                                                                        infectious_holding_in_transmitting_cell_IDs,
                                                                        selected_sus_node_ID)

            #Event-driven stochastic event. If infection occurs, update node status
            if  (rand(rng) <= (cumulative_p / overestimate_cumulative_p_over))

                #SelectedSusNode is infected. Signify infection event on become_exposed_flag indivator vector
                holding_vectors_and_arrays_params.holding_has_had_infection_flag[selected_sus_node_ID] = 1 #Index retreived from nodeID field
            end
        end
    end

    return nothing
end


#===============================================================================
ALGORITHM: Pairwise infection, within same cell (used with conditional subsample algorithm)
===============================================================================#
"""
    local_pairwise_check_sucep_infectious_paris_in_cell!(rng::AbstractRNG,
                                                            n_suscept::Int64,
                                                            holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                            suscept_holding_in_recip_cell_IDs::Array{Int64,1},
                                                            infectious_holding_in_recip_cell_IDs::Array{Int64,1},
                                                            coord_type::Int64,
                                                            kernel_lookup_vec::Array{Float64,1},
                                                            delta_t::Float64)

Check for infection across infectious-susceptible pair in the same cell.

Inputs:
- `rng::AbstractRNG`: The random number generator.
- `n_suscept::Int64`: Number of susceptibles in recipient cell.
- `holding_vectors_and_arrays_params::HoldingVectorsAndArrays`: Set of vectors and arrays associated with holding attributes.
- `suscept_holding_in_recip_cell_IDs::Array{Int64,1}`: For recipient cell, list of nodes that are susceptible.
- `infectious_holding_in_recip_cell_IDs::Array{Int64,1}`: For recipient cell, list of nodes that are infectious.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `delta_t::Float64`: Timestep increment.

Outputs: None \n
Location: LocalSpreadFns.jl
"""
function local_pairwise_check_sucep_infectious_paris_in_cell!(rng::AbstractRNG,
                                                                n_suscept::Int64,
                                                                holding_vectors_and_arrays_params::HoldingVectorsAndArrays,
                                                                suscept_holding_in_recip_cell_IDs::Array{Int64,1},
                                                                infectious_holding_in_recip_cell_IDs::Array{Int64,1},
                                                                coord_type::Int64,
                                                                kernel_lookup_vec::Array{Float64,1},
                                                                delta_t::Float64)

    #Iterate over susceptible holdings
    for suscept_holding_idx = 1:n_suscept

        #Assign current susceptible node of interest to variable
        selected_suscept_holding_ID = suscept_holding_in_recip_cell_IDs[suscept_holding_idx]

        #Get probability susceptible unit infected by any infectious node in cell
        cumul_prob_against_suscept_holding = compute_pairwise_cumulative_infection_prob(holding_vectors_and_arrays_params,
                                                                                coord_type,
                                                                                kernel_lookup_vec,
                                                                                delta_t,
                                                                                infectious_holding_in_recip_cell_IDs,
                                                                                selected_suscept_holding_ID)

        #Draw random number. If less than cumul_prob_against_suscept_holding, infection event was succesful
        if rand(rng) < cumul_prob_against_suscept_holding
            #Update status
            holding_vectors_and_arrays_params.holding_has_had_infection_flag[selected_suscept_holding_ID] = 1
        end
    end

    return nothing
end
