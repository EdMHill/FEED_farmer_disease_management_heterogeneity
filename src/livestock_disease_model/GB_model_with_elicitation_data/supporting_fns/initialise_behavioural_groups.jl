#=
Purpose:
File to house supporting functions for assigning each holding to a behavioural group
and specify the distance threshold the intervention will be enacted.

Function list:
- assign_three_behavioural_groups_empirical_data!
- assign_four_behavioural_groups_empirical_data!
- all_early_reaction_behavioural_group! (holding implements interventions when reported case within 200 miles/320km)
- all_later_reaction_behavioural_group! (holding implements interventions when reported case within 30 miles/50km)

    These two functions also use "assign_holding_intervention_timing_group_and_livestock_house_status!")
- uniform_three_intervention_timing_groups! (third of holdings assigned to each of precautionary, early reactors, late reactors)
- uniform_four_intervention_stance_groups! (quarter of holdings assigned to each of precautionary, early reactors, late reactors, non-users)
=#

"""
    assign_three_behavioural_groups_empirical_data!(rng::AbstractRNG,
                                                    holding_info_vec::Array{HoldingInfo,1},
                                                    sample_from_distribution_flag::Bool)

Assign three behavioural groups and intervention response based on empirical data
from model using the five most stable variables.
Applied to all of Great Britain (homogoneous across nations).

Behavioural group assignment has a dependence on herd size.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.
- `sample_from_distribution_flag::Bool`: Specify whether values are taken from empirical estimates (false) or sampled from distribution (true)

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function assign_three_behavioural_groups_empirical_data!(rng::AbstractRNG,
                                                         holding_info_vec::Array{HoldingInfo,1})

    #=================
    SET UP BEHAVIOURAL GROUP & INTERVENTION TIMING DISTRIBUTION EMPIRICAL ESTIMATES
    AND UNCERTAINTY ASSOCIATED PARAMETERS
    =================#

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Specify proportions and counts in each behavioural group
    # Entries: [Group 1; Group 2; Group 3]
    behavioural_group_propns_three_group_empirical_estimates = [0.43;0.42;0.15]
    n_behavioural_groups = length(behavioural_group_propns_three_group_empirical_estimates)

    # Per behavioural group, the split in intervention timing
    # Column per behavioural group, row per intervention time group
    # Four distinct trigger times
    # 1 - From first emergence in GB (i.e. from outset of simulation)
    # 2 - Within 200 miles (we use 320km)
    # 3 - Within 30 miles (we use 50km)
    # 4 - Never (non-users)
    intervention_timing_propns_group_1_empirical_estimates = [0.46;0.39;0.15;0.00]
    intervention_timing_propns_group_2_empirical_estimates = [0.40;0.32;0.20;0.08]
    intervention_timing_propns_group_3_empirical_estimates = [0.22;0.22;0.56;0.00]
    n_intervention_timing_groups = length(intervention_timing_propns_group_1_empirical_estimates)

    # Per behavioural group, the proportion who house livestock conditional on
    # also vaccinating
    # Entries: [Group 1; Group 2; Group 3]
    also_house_livestock_propn_per_behavioural_group_empirical_estimates = [0.54;0.30;0.78]

    # Specify herd size interval threshold values
    #herd_size_bin_upper_bounds = [100;200;400;600;800;10000]
    herd_size_bin_upper_bounds = [100;250;500;1000;10000]
    n_herd_size_bins = length(herd_size_bin_upper_bounds)

    # Load the data for herd size and behavioural cluster assignment
    herd_size_behavioural_group_data = readdlm("../../../data/elicitation_interviews/cluster_assignment_five_variable_model_with_herd_size.csv", ',', Int64)
    n_data_records = size(herd_size_behavioural_group_data,1)
    behavioural_group_data_col = herd_size_behavioural_group_data[:,1]
    herd_size_data_col = herd_size_behavioural_group_data[:,2]
        # Have disaggregated the behavioural group data and herd size data columns
        # into their own variables

    # Error check. If herd size in loaded data exceeds the final entry in herd_size_bin_upper_bounds,
    # exit the simulation
    if maximum(herd_size_data_col) > herd_size_bin_upper_bounds[end]
        error("herd_size_data_col maximum value, $(maximum(herd_size_data_col)), exceeds final entry in herd_size_bin_upper_bounds ($(herd_size_bin_upper_bounds[end])).")
    end

    # From the data, construct the array behavioural_group_counts_three_group_by_herd_size_empirical_data
    # Row per behavioural group, column per herd size group
    behavioural_group_counts_three_group_by_herd_size_empirical_data = zeros(Int64,n_behavioural_groups,n_herd_size_bins) # Initialise the array
    for data_record_itr = 1:n_data_records

        # Get the herd size group the data record resides in
        data_record_herd_size_bin_idx = 1
        while herd_size_data_col[data_record_itr] > herd_size_bin_upper_bounds[data_record_herd_size_bin_idx]
            # Until the herd size no longer exceeds herd_size_bin_upper_bounds[data_record_herd_size_bin_idx],
            # increment the herd size group index,
            data_record_herd_size_bin_idx += 1
        end

        # Get the behavioural group for this data record
        # Increment the relevant entry in behavioural_group_counts_three_group_by_herd_size_empirical_data
        data_record_behavioural_group = behavioural_group_data_col[data_record_itr]
        behavioural_group_counts_three_group_by_herd_size_empirical_data[data_record_behavioural_group,data_record_herd_size_bin_idx] += 1
    end

    # Use the counts array to compute behavioural_group_propns_three_group_by_herd_size_empirical_estimates
    # Gives the proportion in each behavioural group within a given herd size grouping
    # Done by taking each row, dividing each entry in the row by the sum of that row
    # Row per behavioural group, column per herd size group
    behavioural_group_propns_three_group_by_herd_size_empirical_estimates =
            behavioural_group_counts_three_group_by_herd_size_empirical_data./sum(behavioural_group_counts_three_group_by_herd_size_empirical_data,dims=1)

    # Have those in final group all be in group 3
    behavioural_group_propns_three_group_by_herd_size_empirical_estimates[:,end] = [0;0;1.]

    # Behavioural group assignment has dependence on herd size group the holding
    # resides within
    # behavioural_group_propns_three_group_by_herd_size_empirical_estimates =
    #         [0.52381  0.666667   0.2       0.0  0.0;
    #          0.47619  0.277778   0.466667  0.5  0.0;
    #          0.0      0.0555556  0.333333  0.5  1.0]

     # Error check
     # Number of rows of behavioural group propn array must match number of herd size bins/intervals
     if size(behavioural_group_propns_three_group_by_herd_size_empirical_estimates,2) != n_herd_size_bins
         error("Number of rows in behavioural_group_propns_three_group_by_herd_size_empirical_estimates array: $(size(behavioural_group_propns_three_group_by_herd_size_empirical_estimates,2)). Does not equal number of herd size bins: $(n_herd_size_bins)")
     end

     # Error check
     # Sum of each row of behavioural_group_propns_three_group_by_herd_size_empirical_estimates
     # should be 1. Sum of rows meeting that condition should equal n_herd_size_bins
     if sum(sum(behavioural_group_propns_three_group_by_herd_size_empirical_estimates,dims=1) .== 1) != n_herd_size_bins
         error("Invalid entries in behavioural_group_propns_three_group_by_herd_size_empirical_estimates")
     end

    #=================
    ASSIGN EACH HOLDING TO A BEHAVIOURAL GROUP, AN INTERVENTION TIMING GROUP & HOUSING LIVESTOCK STATUS
    =================#
    for holding_ID = 1:n_holdings

        # Get number of cattle at the holding
        holding_herd_size = holding_info_vec[holding_ID].livestock_count[1]

        # Error check. If herd size exceeds the final entry in herd_size_bin_upper_bounds,
        # exit the simulation
        if holding_herd_size > herd_size_bin_upper_bounds[end]
            error("holding_herd_size for holding_ID $holding_ID is $holding_herd_size. Exceeds final entry in herd_size_bin_upper_bounds ($(herd_size_bin_upper_bounds[end])).")
        end

        # From the herd size, allocate the behavioural group.
        holding_herd_size_bin_idx = 1
        while holding_herd_size > herd_size_bin_upper_bounds[holding_herd_size_bin_idx]
            # Until the herd size no longer exceeds herd_size_bin_upper_bounds[holding_herd_size_bin_idx],
            # increment the herd size group index,
            holding_herd_size_bin_idx += 1
        end

        # Extract relevant behavioural group proportions for the given herd size
        behavioural_group_propns_to_sample_from = behavioural_group_propns_three_group_by_herd_size_empirical_estimates[:,holding_herd_size_bin_idx]

        # Assign behavioural group ID
        holding_info_vec[holding_ID].behavioural_group_ID = sample(rng,[1,2,3], Weights(behavioural_group_propns_to_sample_from), 1)[1]

        # For given behavioural group, sample and assign intervention timing
        if holding_info_vec[holding_ID].behavioural_group_ID == 1
            holding_intervention_timing_group_ID =
                sample(rng,[1,2,3,4], Weights(intervention_timing_propns_group_1_empirical_estimates), 1)[1]
        elseif holding_info_vec[holding_ID].behavioural_group_ID == 2
            holding_intervention_timing_group_ID =
                sample(rng,[1,2,3,4], Weights(intervention_timing_propns_group_2_empirical_estimates), 1)[1]
        elseif holding_info_vec[holding_ID].behavioural_group_ID == 3
            holding_intervention_timing_group_ID =
                sample(rng,[1,2,3,4], Weights(intervention_timing_propns_group_3_empirical_estimates), 1)[1]
        else
            error("holding_info_vec[$holding_ID].behavioural_group_ID has value $(holding_info_vec[holding_ID].behavioural_group_ID). Invalid value.")
        end
        holding_info_vec[holding_ID].intervention_timing_group_ID = holding_intervention_timing_group_ID
            # Assignment to HoldingInfo structure

        # Assign intervention threshold distance (in metres)
        if holding_intervention_timing_group_ID == 1
            holding_info_vec[holding_ID].intervention_threshold_distance = 1e10
        elseif holding_intervention_timing_group_ID == 2
            # Threshold: 320km (approx. 200 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 320000
        elseif holding_intervention_timing_group_ID == 3
            # Threshold: 50km (approx. 30 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 50000
        elseif holding_intervention_timing_group_ID == 4
            # Interventions are never used
            holding_info_vec[holding_ID].intervention_threshold_distance = -1
        end

        # For given intervention timing group, assign housing livestock flag value
        if holding_intervention_timing_group_ID < 4
            if rand(rng) < also_house_livestock_propn_per_behavioural_group_empirical_estimates[holding_info_vec[holding_ID].behavioural_group_ID]
                holding_info_vec[holding_ID].house_livestock_flag = true
            else
                holding_info_vec[holding_ID].house_livestock_flag = false
            end
        else # In a non-vaccintion and non-housing intervention group.
            holding_info_vec[holding_ID].house_livestock_flag = false
        end
    end

    return nothing
end

"""
    assign_four_behavioural_groups_empirical_data!(rng::AbstractRNG,
                                                    holding_info_vec::Array{HoldingInfo,1},
                                                    sample_from_distribution_flag::Bool)

Assign four behavioural groups and intervention response based on empirical data
from model using the two most stable variables.
Applied to all of Great Britain (homogoneous across nations).

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.
- `sample_from_distribution_flag::Bool`: Specify whether values are taken from empirical estimates (false) or sampled from distribution (true)

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function assign_four_behavioural_groups_empirical_data!(rng::AbstractRNG,
                                                         holding_info_vec::Array{HoldingInfo,1})

    #=================
    SET UP BEHAVIOURAL GROUP & INTERVENTION TIMING DISTRIBUTION EMPIRICAL ESTIMATES
    AND UNCERTAINTY ASSOCIATED PARAMETERS
    =================#

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Construct holding vector access reference vector
    holding_idx_reference_vec = collect(1:1:n_holdings)

    # Specify proportions and counts in each behavioural group
    # Entries: [Group 1; Group 2; Group 3; Group 4]
    behavioural_group_propns_four_group_empirical_estimates = [0.43;0.25;0.17;0.15]
    behavioural_group_propns_four_group_sd = [0.05;0.06;0.05;0.03]
    n_behavioural_groups = length(behavioural_group_propns_four_group_empirical_estimates)

    # Per behavioural group, the split in intervention timing
    # Column per behavioural group, row per intervention time group
    # Four distinct trigger times
    # 1 - From first emergence in GB (i.e. from outset of simulation)
    # 2 - Within 200 miles (we use 320km)
    # 3 - Within 30 miles (we use 50km)
    # 4 - Never (non-users)
    intervention_timing_propns_group_1_empirical_estimates = [0.46;0.39;0.15;0.00]
    intervention_timing_propns_group_2_empirical_estimates = [0.40;0.40;0.20;0.00]
    intervention_timing_propns_group_3_empirical_estimates = [0.40;0.10;0.50;0.00]
    intervention_timing_propns_group_4_empirical_estimates = [0.22;0.34;0.22;0.22]
    n_intervention_timing_groups = length(intervention_timing_propns_group_1_empirical_estimates)

    # intervention_timing_propns_group_1_sd = [0.09;0.09;0.10;0.00]
    # intervention_timing_propns_group_2_sd = [0.1;0.1;0.05;0.03]
    # intervention_timing_propns_group_3_sd = [0.1;0.1;0.15;0.00]
    # intervention_timing_propns_group_4_sd = [0.1;0.1;0.15;0.00]

    # Per behavioural group, the proportion who house livestock conditional on
    # also vaccinating
    # Entries: [Group 1; Group 2; Group 3: Group 4]
    also_house_livestock_propn_per_behavioural_group_empirical_estimates = [0.54;0.20;0.60;0.71]
    # also_house_livestock_propn_per_behavioural_group_sd = [0.1;0.075;0.1;0.00]

    # For this simulation run, sample values for behavioural groups, intervention timing
    # and proportion housing livestock (conditional on also vaccinating)
    behavioural_group_propns_GB = 1. *behavioural_group_propns_four_group_empirical_estimates
    intervention_timing_propns_group_1 = 1. *intervention_timing_propns_group_1_empirical_estimates
    intervention_timing_propns_group_2 = 1. *intervention_timing_propns_group_2_empirical_estimates
    intervention_timing_propns_group_3 = 1. *intervention_timing_propns_group_3_empirical_estimates
    intervention_timing_propns_group_4 = 1. *intervention_timing_propns_group_4_empirical_estimates
    also_house_livestock_propn_per_behavioural_group = 1. *also_house_livestock_propn_per_behavioural_group_empirical_estimates

    # Get count per behavioural group
    n_per_behavioural_group_GB = round.(Int64,n_holdings.*behavioural_group_propns_GB)
    n_per_behavioural_group_GB[2] -= 1 # Through rounding, increases total count by 1. Manually adjust count in second group to remedy this.

    # Consolidate intervention timings into a single array
    intervention_timing_group_propns = hcat(intervention_timing_propns_group_1,
                                            intervention_timing_propns_group_2,
                                            intervention_timing_propns_group_3,
                                            intervention_timing_propns_group_4)

    # Get proportion of total population per behavioural group and intervention timing
    propn_per_behavioural_and_intervention_timing_group = behavioural_group_propns_GB'.*intervention_timing_group_propns

    # Get population size per behavioural group (overall and per intervention timing group)
    count_per_behavioural_and_intervention_timing_group = ceil.(Int64,n_holdings.*propn_per_behavioural_and_intervention_timing_group)

    # Error check
    if sum(n_per_behavioural_group_GB) != n_holdings
        error("sum(n_per_behavioural_group_GB): $(sum(n_per_behavioural_group_GB)); n_holdings: $(n_holdings). Should be equal.")
    end

    # Compute the proportion who use vaccination in each behavioural group
    # Corresponds to first three rows in intervention_timing_group_propns
    # Then compute the number who house livestock.
    propn_who_vaccinate_per_behavioural_group = vec(sum(intervention_timing_group_propns[1:3,:],dims=1))
    propn_who_house_livestock_per_behavioural_and_intervention_timing_group = also_house_livestock_propn_per_behavioural_group'.*propn_per_behavioural_and_intervention_timing_group[1:3,:]
    count_who_house_livestock_per_behavioural_and_intervention_timing_group = round.(Int64,n_holdings.*propn_who_house_livestock_per_behavioural_and_intervention_timing_group)

    #=================
    ASSIGN EACH HOLDING TO A BEHAVIOURAL GROUP, INTERVENTION TIMING GROUP & HOUSING LIVESTOCK STATUS
    =================#

    # Initialise holding attribute array
    # Row per holding
    # Columns: Behavioural group; intervention timing; if also houses livestock.
    holding_attribute_array = zeros(Int64,n_holdings,3)

    # Behavioural group assignment for each holding in country
    # First construct vector with count of each behavioural group ID reflecting
    # the proportional split in the population
    behavioural_group_IDs_to_sample = Int64[]
    for behavioural_group_idx = 1:n_behavioural_groups
        append!(behavioural_group_IDs_to_sample,behavioural_group_idx*ones(Int64,n_per_behavioural_group_GB[behavioural_group_idx]))
    end
        # behavioural_group_IDs_to_sample = [ones(Int64,n_per_behavioural_group_GB[1]);
        #                                    2*ones(Int64,n_per_behavioural_group_GB[2]);
        #                                    3*ones(Int64,n_per_behavioural_group_GB[3]);
        #                                    4*ones(Int64,n_per_behavioural_group_GB[4])]

    # Then sample and assign to first column of holding_attribute_array
    holding_behavioural_group_IDs_vec = sample(rng,
                                                 behavioural_group_IDs_to_sample,
                                                 n_holdings,
                                                 replace=false)
    holding_attribute_array[:,1] = holding_behavioural_group_IDs_vec

    # To each behavioural group, assign to each holding an intervention timing
    # and whether the holding will also house livestock
    for behavioural_group_idx = 1:n_behavioural_groups

       # Construct vector with count of each intervetion timing group reflecting
       # the proportional split in the population
       intervention_timing_group_IDs_to_sample = Int64[]
       for intervention_timing_idx = 1:n_intervention_timing_groups
           append!(intervention_timing_group_IDs_to_sample,intervention_timing_idx*ones(Int64,count_per_behavioural_and_intervention_timing_group[intervention_timing_idx,behavioural_group_idx]))
       end

       # Get IDs of holdings in behavioural_group_idx
       # Sample and assign to second column of holding_attribute_array
       holding_IDs_in_this_behavioural_group = holding_idx_reference_vec[holding_behavioural_group_IDs_vec.== behavioural_group_idx]
       holding_intervention_group_IDs_for_this_behavioural_group_vec = sample(rng,
                                                                              intervention_timing_group_IDs_to_sample,
                                                                              length(holding_IDs_in_this_behavioural_group),
                                                                              replace=false)
       holding_attribute_array[holding_IDs_in_this_behavioural_group,2] = holding_intervention_group_IDs_for_this_behavioural_group_vec

        # Assign those holdings that also house livestock
        for intervention_timing_group_idx = 1:(n_intervention_timing_groups-1)
                # Ignore final intervention timing group, as they do not use interventions at any time

                # Get IDs of holdings that are in the applicable country,
                # behavioural group and intervention timing group
                holding_IDs_in_this_behavioural_group_and_intervention_timing_group = holding_IDs_in_this_behavioural_group[holding_intervention_group_IDs_for_this_behavioural_group_vec.== intervention_timing_group_idx]

                # Construct vector of values to be sampled from
                house_livestock_bool_vec = [ones(Bool,count_who_house_livestock_per_behavioural_and_intervention_timing_group[intervention_timing_group_idx,behavioural_group_idx]);
                                            zeros(Bool,count_per_behavioural_and_intervention_timing_group[intervention_timing_group_idx,behavioural_group_idx] -
                                                        count_who_house_livestock_per_behavioural_and_intervention_timing_group[intervention_timing_group_idx,behavioural_group_idx])]

                holding_attribute_array[holding_IDs_in_this_behavioural_group_and_intervention_timing_group,3] = sample(rng,
                                                                                                                          house_livestock_bool_vec,
                                                                                                                          length(holding_IDs_in_this_behavioural_group_and_intervention_timing_group),
                                                                                                                          replace=false)
        end
    end

    #=================
    ITERATE OVER EACH HOLDING AND ASSIGN INFORMATION TO HOLDING_INFO_VEC
    =================#
    for holding_ID = 1:n_holdings

        # Assign behavioural group ID
        holding_info_vec[holding_ID].behavioural_group_ID = holding_attribute_array[holding_ID,1]

        # Assign intervention timing group distance
        holding_intervention_timing_group_ID = holding_attribute_array[holding_ID,2]
        holding_info_vec[holding_ID].intervention_timing_group_ID = holding_intervention_timing_group_ID

        # Assign intervention threshold distance (in metres)
        if holding_intervention_timing_group_ID == 1
            holding_info_vec[holding_ID].intervention_threshold_distance = 1e10
        elseif holding_intervention_timing_group_ID == 2
            # Threshold: 320km (approx. 200 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 320000
        elseif holding_intervention_timing_group_ID == 3
            # Threshold: 50km (approx. 30 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 50000
        elseif holding_intervention_timing_group_ID == 4
            # Interventions are never used
            holding_info_vec[holding_ID].intervention_threshold_distance = -1
        end

        # Assign housing livestock flag value
        holding_info_vec[holding_ID].house_livestock_flag = holding_attribute_array[holding_ID,3]
    end

    return nothing
end

"""
    all_early_reaction_behavioural_group!(rng::AbstractRNG,holding_info_vec::Array{HoldingInfo,1})

Assign a single intervention response, when reported infection is within 320km (approx 200 miles).

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function all_early_reaction_behavioural_group!(rng::AbstractRNG,
                                                holding_info_vec::Array{HoldingInfo,1})

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Iterate over each holding & assign threshold distance
    for holding_ID = 1:n_holdings

         # Assign intervention timing group distance
         holding_info_vec[holding_ID].intervention_timing_group_ID = 2
         holding_info_vec[holding_ID].intervention_threshold_distance = 320000

         # Assign eacch holding as also using housing
         holding_info_vec[holding_ID].house_livestock_flag = true
    end

     return nothing
end

"""
    all_later_reaction_behavioural_group!(rng::AbstractRNG,holding_info_vec::Array{HoldingInfo,1})

Assign a single intervention response, when reported infection is within 50km (approx 30 miles).

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function all_later_reaction_behavioural_group!(rng::AbstractRNG,
                                                holding_info_vec::Array{HoldingInfo,1})

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Iterate over each holding & assign threshold distance
    for holding_ID = 1:n_holdings

         # Assign intervention timing group distance
         holding_info_vec[holding_ID].intervention_timing_group_ID = 3
         holding_info_vec[holding_ID].intervention_threshold_distance = 50000

         # Assign eacch holding as also using housing
         holding_info_vec[holding_ID].house_livestock_flag = true
    end

     return nothing
end

"""
    uniform_three_intervention_timing_groups!(rng::AbstractRNG,holding_info_vec::Array{HoldingInfo,1})

Split holdings uniformly across three intervention timing groups.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function uniform_three_intervention_timing_groups!(rng::AbstractRNG,
                                                   holding_info_vec::Array{HoldingInfo,1})

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Split in intervention timing
    # Four distinct trigger times
    # 1 - From first emergence in GB (i.e. from outset of simulation)
    # 2 - Within 200 miles (we use 320km)
    # 3 - Within 30 miles (we use 50km)
    # 4 - Never (non-users)
    intervention_timing_propns_uniform_three_groups = [1/3;1/3;1/3;0.00]


    # Specify the proportion of this behavioural group who will also house livestock
    also_house_livestock_propn_uniform_three_groups = 1.

    # Call function to update holding_info_vec with intervention timing group
    # and housing livestock flag value
    behavioural_group_ID_input_uniform_three_groups = 1
    assign_holding_intervention_timing_group_and_livestock_house_status!(rng,
                                                                         holding_info_vec,
                                                                         n_holdings,
                                                                         behavioural_group_ID_input_uniform_three_groups,
                                                                         intervention_timing_propns_uniform_three_groups,
                                                                         also_house_livestock_propn_uniform_three_groups)


     return nothing
end

"""
    uniform_four_intervention_stance_groups!(rng::AbstractRNG,holding_info_vec::Array{HoldingInfo,1})

Split holdings uniformly across four intervention stance groups, including non-users of vaccination and housing.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function uniform_four_intervention_stance_groups!(rng::AbstractRNG,
                                                  holding_info_vec::Array{HoldingInfo,1})

    # Get number of holdings
    n_holdings = length(holding_info_vec)

    # Split in intervention timing
    # Four distinct trigger times
    # 1 - From first emergence in GB (i.e. from outset of simulation)
    # 2 - Within 200 miles (we use 320km)
    # 3 - Within 30 miles (we use 50km)
    # 4 - Never (non-users)
    intervention_timing_propns_uniform_four_groups = [0.25;0.25;0.25;0.25]


    # Specify the proportion of this behavioural group who will also house livestock
    # (from those groups that are not non-vaccinators)
    also_house_livestock_propn_uniform_four_groups = 1.

    # Call function to update holding_info_vec with intervention timing group
    # and housing livestock flag value
    behavioural_group_ID_input_uniform_four_groups = 1
    assign_holding_intervention_timing_group_and_livestock_house_status!(rng,
                                                                         holding_info_vec,
                                                                         n_holdings,
                                                                         behavioural_group_ID_input_uniform_four_groups,
                                                                         intervention_timing_propns_uniform_four_groups,
                                                                         also_house_livestock_propn_uniform_four_groups)


     return nothing
end

"""
    assign_holding_intervention_timing_group_and_livestock_house_status!(rng::AbstractRNG,
                                                                         holding_info_vec::Vector{HoldingInfo},
                                                                         n_holdings::Int64,
                                                                         behavioural_group_ID_input::Int64,
                                                                         intervention_timing_group_propns::Vector{Float64},
                                                                         also_house_livestock_propn::Float64)

Assign to holdings the intervention timing group and flag value fo housing livestock.

Inputs:
- `rng::AbstractRNG`: Random number generator.
- `holding_info_vec::Array{HoldingInfo,1}`: Holding-level attributes.
- `n_holdings::Int64`: Number of holdings in the population.
- `behavioural_group_ID_input::Int64`: Behavioural group ID to be assigned to all holdings.
- `intervention_timing_group_propns::Vector{Float64}`: Proportion in each intervention timing group.
- `also_house_livestock_propn::Float64`: Proportion of behavioural group that uses housing of livestock intervention.

Outputs: None \n
Location: initialise\\_behavioural\\_groups.jl
"""
function assign_holding_intervention_timing_group_and_livestock_house_status!(rng::AbstractRNG,
                                                                             holding_info_vec::Vector{HoldingInfo},
                                                                             n_holdings::Int64,
                                                                             behavioural_group_ID_input::Int64,
                                                                             intervention_timing_group_propns::Vector{Float64},
                                                                             also_house_livestock_propn::Float64)
    # Get count & cumulative sum of total population per intervention timing
    count_per_intervention_timing_group = ceil.(Int64,n_holdings.*intervention_timing_group_propns)
    csum_count_per_intervention_time_group = cumsum(count_per_intervention_timing_group)

    # Compute the proportion who use vaccination in the behavioural group
    # Corresponds to first three rows in intervention_timing_group_propns
    # Then compute the number who house livestock.
    propn_who_vaccinate_in_the_behavioural_group = sum(intervention_timing_group_propns[1:3])
    count_who_also_house_livestock_in_this_behavioural_group = round(Int64,n_holdings.*propn_who_vaccinate_in_the_behavioural_group.*also_house_livestock_propn)

    #=================
    ASSIGN EACH HOLDING TO INTERVENTION TIMING GROUP & HOUSING LIVESTOCK STATUS
    =================#

    # Initialise holding attribute array
    # Row per holding
    # Columns: Intervention timing; if also houses livestock.
    holding_attribute_array = zeros(Int64,n_holdings,2)

    # Get intervention timing group IDs & assign to holding_attribute_array
    intervention_timing_ID_vec_to_sample_from = [ones(Int64,count_per_intervention_timing_group[1]);
                                                         2*ones(Int64,count_per_intervention_timing_group[2]);
                                                         3*ones(Int64,count_per_intervention_timing_group[3]);
                                                         4*ones(Int64,count_per_intervention_timing_group[4])]
    holding_attribute_array[:,1] = sample(rng,intervention_timing_ID_vec_to_sample_from,n_holdings,replace=false)

    # For determining holdings house livestock status, construct vector of values to be sampled from
    idx_vec = collect(1:1:n_holdings)
    eligible_holding_ID_to_house_livestock = idx_vec[holding_attribute_array[:,1].<4]
    n_eligible_holdings_to_house_livestock = length(eligible_holding_ID_to_house_livestock)
    house_livestock_bool_vec = [ones(Bool,count_who_also_house_livestock_in_this_behavioural_group);
                                zeros(Bool,n_eligible_holdings_to_house_livestock - count_who_also_house_livestock_in_this_behavioural_group)]

    # Assign to holding_attribute_array the holdings that also house livestock
    holding_attribute_array[eligible_holding_ID_to_house_livestock,2] = sample(rng,
                                                                              house_livestock_bool_vec,
                                                                              n_eligible_holdings_to_house_livestock,
                                                                              replace=false)

    #=================
    ITERATE OVER EACH HOLDING AND ASSIGN INFORMATION TO HOLDING_INFO_VEC
    =================#
    for holding_ID = 1:n_holdings

        # Assign behavioural group ID
        holding_info_vec[holding_ID].behavioural_group_ID = behavioural_group_ID_input

        # Assign intervention timing group distance
        holding_intervention_timing_group_ID = holding_attribute_array[holding_ID,1]
        holding_info_vec[holding_ID].intervention_timing_group_ID = holding_intervention_timing_group_ID

        # Assign intervention threshold distance (in metres)
        if holding_intervention_timing_group_ID == 1
            holding_info_vec[holding_ID].intervention_threshold_distance = 1e10
        elseif holding_intervention_timing_group_ID == 2
            # Threshold: 320km (approx. 200 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 320000
        elseif holding_intervention_timing_group_ID == 3
            # Threshold: 50km (approx. 30 miles)
            holding_info_vec[holding_ID].intervention_threshold_distance = 50000
        elseif holding_intervention_timing_group_ID == 4
            # Interventions are never used
            holding_info_vec[holding_ID].intervention_threshold_distance = -1
        end

        # Assign housing livestock flag value
        holding_info_vec[holding_ID].house_livestock_flag = holding_attribute_array[holding_ID,2]
    end

    return nothing
end