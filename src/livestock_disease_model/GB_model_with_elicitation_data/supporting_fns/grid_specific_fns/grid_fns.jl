#=
Purpose:
File to house functions used in performing the simulation of infection using
conditional subsample grid algorithm

Alogorithm list:
- Check landscape & holding data are compatible
- Return landscaspe size attributes (& check values are valid)
- Return cell level attributes (susceptibility, transmissibility, number of holding)
- Construct binomial RNG
- Precalculate cell-to-cell maximum transmission rate
- Precalculate cell-to-cell maximum transmission probabilities
- Construct susceptible & infectious prem ID & location tuples
- Static grid construction
=#

#-------------------------------------------------------------------------------
### CHECK LANDSCAPE & holding DATA ARE COMPATIBLE
#-------------------------------------------------------------------------------
"""
    check_landscape_valid(bounding_box_var::Array{Float64,1},
                            holding_loc_X_vals::Array{Float64,1},
                            holding_loc_Y_vals::Array{Float64,1})

Check landscape and holding data are compatible.

Inputs:
- `bounding_box_var::Array{Float64,1}`: Vector with limits for landscape bounding box. Entries [Min_x,Max_x,Min_y,Max_y].
- `holding_loc_X_vals::Array{Float64,1}`: East-west plane co-ordinate (per holding).
- `holding_loc_Y_vals::Array{Float64,1}`: North-south plane co-ordinate (per holding).

Outputs: None \n
Location: grid\\_fns.jl
"""
function check_landscape_valid(bounding_box_var::Array{Float64,1},
                            holding_loc_X_vals::Array{Float64,1},
                            holding_loc_Y_vals::Array{Float64,1})

    #Function designed to throw error and exit programme if holding lie on boundary/outside
    #landscape bounding box

    #CHECK BOUNDING BOX WIDTH AND HEIGHT WILL BE POSITIVE!
    xMin =  bounding_box_var[1]; xMax = bounding_box_var[2];
    yMin =  bounding_box_var[3]; yMax = bounding_box_var[4];
    bounding_box_width = xMax - xMin
    bounding_box_height = yMax - yMin

    println("holding_loc_X_vals min: $(minimum(holding_loc_X_vals))")
    println("holding_loc_X_vals max: $(maximum(holding_loc_X_vals))")
    println("holding_loc_Y_vals min: $(minimum(holding_loc_Y_vals))")
    println("holding_loc_Y_vals max: $(maximum(holding_loc_Y_vals))")

    if bounding_box_width <= 0 && bounding_box_height <= 0
        error("Nonpositive bounding_box_width and bounding_box_height found.")
    elseif bounding_box_width <= 0
        error("Nonpositive bounding_box_width found.")
    elseif bounding_box_height <= 0
        error("Nonpositive bounding_box_height found.")
    end

    #CHECK THAT NO holding WILL BE LOCATED ON AN EDGE OR OUTSIDE THE BOUNDING BOX
    #Locations on bottom edge and left edge are okay, such holding will be allocated a grid in the grid configuration.
    #holding on top edge and right edge would be missed!
    error_flag = 0 #Initialise variable to determine whether error should be thrown or not
    if (sum(holding_loc_X_vals .<= xMin) > 0)
        println("Incompatible boundary. holding located on or outside left edge of landscape box.")
        error_flag = 1
    end

    if (sum(holding_loc_X_vals .>= xMax) > 0)
        println("Incompatible boundary. holding located on or outside right edge of landscape box.")
        error_flag = 1
    end

    if (sum(holding_loc_Y_vals .<= yMin) > 0)
        println("Incompatible boundary. holding located on or below bottom edge of landscape box.")
        error_flag = 1
    end

    if (sum(holding_loc_Y_vals .>= yMax) > 0)
        println("Incompatible boundary. holding located on or above top edge of landscape box.")
        error_flag = 1
    end

    if error_flag == 1
        error("Incompatible landscape bounding box. Amend boundary values. Exiting programme.")
    end

    return
end
#-------------------------------------------------------------------------------
### GET LANDSCAPE SIZE ATTRIBUTES (& CHECK VALUES ARE VALID)
#-------------------------------------------------------------------------------
"""
    get_landscape_size(landscape_grid_data::GridData,
                        holding_loc_X_vals::Array{Float64,1},
                        holding_loc_Y_vals::Array{Float64,1})

Get landscape attributes and check the attributes are valid.

Inputs:
- `landscape_grid_data::GridData`: Composite type with info such as co-ordinate system type & bounding box information [Min_x,Max_x,Min_y,Max_y]
- `holding_loc_X_vals::Array{Float64,1}`: East-west plane co-ordinate (per holding).
- `holding_loc_Y_vals::Array{Float64,1}`: North-south plane co-ordinate (per holding).

Outputs:
- `bounding_box_width::Float64`: Width of the bounding box for the landscape.
- `bounding_box_height::Float64`: Height of the bounding box for the landscape.
- `longest_landscape_edge::Float64`: Longest edge of the bounding box for the landscape.
- `max_dist_within_landscape::Int64`: Distance of diagonal across landscape (maximum separation between two points within the landscape).

Location: grid\\_fns.jl
"""
function get_landscape_size(landscape_grid_data::GridData,
                            holding_loc_X_vals::Array{Float64,1},
                            holding_loc_Y_vals::Array{Float64,1})

    #Calculate length of longest edge
    bounding_box_width = landscape_grid_data.bounding_box_vals[2] - landscape_grid_data.bounding_box_vals[1]
    bounding_box_height = landscape_grid_data.bounding_box_vals[4] - landscape_grid_data.bounding_box_vals[3]
    longest_landscape_edge = max(bounding_box_width,bounding_box_height)
    println("longest_landscape_edge is $longest_landscape_edge. CHECK THIS IS CORRECT!")

    #Calculate distance (in metres) along diagonal for landscape
    if landscape_grid_data.coord_type == 1 #Cartesian co-ordinate system. Gives distance between nodes in metres.
        max_dist_within_landscape_float = sqrt(longest_landscape_edge*longest_landscape_edge + longest_landscape_edge*longest_landscape_edge)
    elseif landscape_grid_data.coord_type == 2 #Cartesian co-ordinate system. Gives distance between nodes in km. Convert to metres.
        max_dist_within_landscape_km = sqrt(longest_landscape_edge*longest_landscape_edge + longest_landscape_edge*longest_landscape_edge)
        max_dist_within_landscape_float = max_dist_within_landscape_km*1000
    else #LatLong co-ordinate system. Output in km. Convert to metres.
        max_dist_within_landscape_float = GreatCircleDistance(landscape_grid_data.bounding_box_vals[1]+longest_landscape_edge, landscape_grid_data.bounding_box_vals[1],
                                                            landscape_grid_data.bounding_box_vals[2]+longest_landscape_edge, landscape_grid_data.bounding_box_vals[2])
    end

    #Take ceiling of float distance and convert to integer
    max_dist_within_landscape = convert(Int64,ceil(max_dist_within_landscape_float))
    println("max_dist_within_landscape is $max_dist_within_landscape. CHECK THIS IS CORRECT!")

    return bounding_box_width::Float64,
            bounding_box_height::Float64,
            longest_landscape_edge::Float64,
            max_dist_within_landscape::Int64
end

#-------------------------------------------------------------------------------
### RETURN CELL LEVEL ATTRIBUTES (SUSCEPTIBILITY, TRANSMISSIBILITY, NUMBER OF holding)
#-------------------------------------------------------------------------------
"""
    get_cell_attributes(holding_cell_IDs::Array{Int64},
                        holding_loc_X_and_Y_vals::Array{Float64,2},
                        holding_suscept::Vector{Float64},
                        holding_transmiss::Vector{Float64})

Return cell level attributes (susceptibility, transmissibility, number of holding).

Inputs:
- `holding_cell_IDs::Array{Int64}`: ID for each holding.
- `holding_loc_X_and_Y_vals::Array{Float64,2}`: Coordinates for each holding.
- `holding_suscept::Vector{Float64}`: holding-level susceptibility.
- `holding_transmiss::Vector{Float64}`: holding-level transmissibility.

Outputs:
- `n_holdings_per_cell::Vector{Int64}`: Number of holdings in each cell of the grid.
- `holding_in_cell_IDs::Vector{Vector{Int64}}`: For each cell, vector of holding IDs that reside in that cell.
- `max_suscept_per_cell::Vector{Float64}`: Maximum susceptibility per cell of the grid.
- `max_trans_per_cell::Vector{Float64}`: Maximum transmissibility per cell of the grid.

Location: grid\\_fns.jl
"""
function get_cell_attributes(holding_cell_IDs::Array{Int64},
                            holding_loc_X_and_Y_vals::Array{Float64,2},
                            holding_suscept::Vector{Float64},
                            holding_transmiss::Vector{Float64})

    #Get largest grid ID value
    holding_cell_IDs_max_val = maximum(holding_cell_IDs)

    #Initialise vectors for grid attributes
    n_holdings_per_cell = zeros(Int64,holding_cell_IDs_max_val) #number of holding per grid

    #Initialise tuple for storing holding IDs residing in each cell
    holding_in_cell_IDs = Array{Array{Int64,1},1}(undef,holding_cell_IDs_max_val)

    #Initialise grid-level susceptibility and transmissibility vectors
    max_suscept_per_cell = zeros(holding_cell_IDs_max_val)
    max_trans_per_cell = zeros(holding_cell_IDs_max_val)

    #Iterate over grid IDs
    for ii = 1:holding_cell_IDs_max_val
        m = findall(isequal(ii),holding_cell_IDs)::Array{Int64,1}
        holding_in_cell_IDs[ii] = m
        n_holdings_per_cell[ii] = length(m) #Get number of holding in grid ii
        if n_holdings_per_cell[ii]>0
            max_suscept_per_cell[ii] = maximum(holding_suscept[m])
            max_trans_per_cell[ii] = maximum(holding_transmiss[m])
        else
            max_suscept_per_cell[ii] = 0.
            max_trans_per_cell[ii] = 0.
        end
    end

    return  n_holdings_per_cell::Vector{Int64},
            holding_in_cell_IDs::Vector{Vector{Int64}},
            max_suscept_per_cell::Vector{Float64},
            max_trans_per_cell::Vector{Float64}
end

#-------------------------------------------------------------------------------
### CONSTRUCT BINOMIAL RNG ARRAY
#-------------------------------------------------------------------------------
"""
    construct_binomial_RNG(n_holdings_per_cell::Array{Int64,1},
                            P_CS::Array{Float64,1})

Construct a collection of Binomial RNGs (random number generators).

Inputs:
- `n_holdings_per_cell::Array{Int64,1}`: Vector with entry for number of holding per grid cell.
- `P_CS::Array{Float64,1})`: Preset probabilities to initialise RNG with.

Outputs:
- `binomial_RNG_array::Array{Binomial{Float64},2}`: Array of Binomial RNGs. Row per susceptible number, column per probability threshold.

Location: grid\\_fns.jl
"""
function construct_binomial_RNG(n_holdings_per_cell::Array{Int64,1},
                                P_CS::Array{Float64,1})

    #Get maximum possible number of susceptibles per grid.
    MaxPremPerGrid = maximum(n_holdings_per_cell)

    #Get quantity of preset probabilities to initialise RNG with
    PresetProbThresholds = length(P_CS) - 1 #Ignore first zero entry

    #Construct the Binomial RNG
    #Row per susceptible number, column per probability threshold
    binomial_RNG_array = Array{Binomial{Float64},2}(undef,MaxPremPerGrid,PresetProbThresholds) #Initialise array of Binomial RNG
    for ii =1:PresetProbThresholds
        for jj = 1:MaxPremPerGrid
            binomial_RNG_array[jj,ii] = Binomial(jj, P_CS[ii+1])  #First entry of P_CS is zero. Skip that entry.
        end
    end

    return binomial_RNG_array::Array{Binomial{Float64},2}
end

#-------------------------------------------------------------------------------
### PRECALCULATE CELL-TO-CELL MAXIMUM TRANSMISSION PROBABILITIES
#-------------------------------------------------------------------------------
"""
    max_cell_prob_calc(cell_lim_array::Array{Float64,2},
                        n_holdings_per_cell::Array{Int64,1},
                        coord_type::Int64,
                        kernel_lookup_vec::Array{Float64,1},
                        max_suscept_per_cell::Array{Float64,1},
                        max_trans_per_cell::Array{Float64,1},
                        delta_t::Float64)

Precalculate cell-to-cell maximum transmission probabilities.

Inputs:
- `cell_lim_array::Array{Float64,2}`: Vector with limits of each cell (collapsed array). Each block of four elements corresponds to [cell_xMin cell_xMax cell_yMin cell_yMax].
- `n_holdings_per_cell::Array{Int64,1}`: Vector with entry for number of holding per grid cell.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `max_suscept_per_cell::Array{Float64,1}`: Maximum susceptibility per cell of the grid.
- `max_trans_per_cell::Array{Float64,1}`: Maximum transmissibility per cell of the grid.
- `delta_t::Float64`: Timestep increment.

Outputs:
- `max_cell_prob::Array{Float64,2}`: Array of grid to grid transmission probabilities.

Location: grid\\_fns.jl
"""
function max_cell_prob_calc(cell_lim_array::Array{Float64,2},
                            n_holdings_per_cell::Array{Int64,1},
                            coord_type::Int64,
                            kernel_lookup_vec::Array{Float64,1},
                            max_suscept_per_cell::Array{Float64,1},
                            max_trans_per_cell::Array{Float64,1},
                            delta_t::Float64)

    #Get number of grids in use
    holding_cell_IDs_max_val = length(n_holdings_per_cell)

    #Populate max_cell_to_cell_transmission_rate
    max_cell_prob = zeros(holding_cell_IDs_max_val,holding_cell_IDs_max_val)
    for jj = 1:holding_cell_IDs_max_val
        for ii = 1:holding_cell_IDs_max_val
            if ii==jj  #Within-cell event
                max_cell_prob[ii,jj] = 1.
            elseif n_holdings_per_cell[ii]==0 || n_holdings_per_cell[jj]==0 #no holding present in grid
                max_cell_prob[ii,jj] = 0.
            else  #Calculate distance based on type of location data
                distance = 0.0 #Initialise distance variable
                if coord_type == 1 #Cartesian co-ords (metres)
                    distance = eucl_distance_between_cells(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                            cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                            cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                            cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax

                elseif coord_type == 2 #Cartesian co-ords (km)
                    distance = eucl_distance_between_cells_convert_to_metres(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                                cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                                cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                                cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax
                elseif coord_type == 3 #Lat/Long co-ords

                    #Get distance by calling great circle distance function
                    distance = great_circle_distance_between_cells(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                            cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                            cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                            cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax
                end

                #Compute over estimated grid-to-grid transmission rate
                dist_idx = return_dist_idx_for_kernel(distance)
                max_cell_to_cell_transmission_rate = max_trans_per_cell[ii]*max_suscept_per_cell[jj]*kernel_lookup_vec[dist_idx]

                #Computer over-estimated grid-to-grid transmission rate
                max_cell_prob[ii,jj] = oneMinusExp(-max_cell_to_cell_transmission_rate*delta_t)

                if ii==1 && jj==10
                    println("cell_ii, cell_jj, distance: $ii, $jj, $distance")
                    println(kernel_lookup_vec[dist_idx])
                end

                if ii==1 && jj==2
                    println("cell_ii, cell_jj, distance: $ii, $jj, $distance")
                    println(kernel_lookup_vec[dist_idx])
                end
            end
        end
    end

    println("Check max_cell_prob")

    return max_cell_prob::Array{Float64,2}
end

#-------------------------------------------------------------------------------
### PRECALCULATE CELL-TO-CELL MAXIMUM TRANSMISSION RATE
#-------------------------------------------------------------------------------
"""
    max_cell_to_cell_transmission_rate_calc(cell_lim_array::Array{Float64,2},
                    n_holdings_per_cell::Array{Int64,1},
                    coord_type::Int64,
                    kernel_lookup_vec::Array{Float64,1},
                    max_suscept_per_cell::Array{Float64,1},
                    max_trans_per_cell::Array{Float64,1})

Precalculate cell-to-cell maximum transmission rates.

Inputs:
- `cell_lim_array::Array{Float64,2}`: Vector with limits of each cell (collapsed array). Each block of four elements corresponds to [cell_xMin cell_xMax cell_yMin cell_yMax].
- `n_holdings_per_cell::Array{Int64,1}`: Vector with entry for number of holding per grid cell.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").
- `kernel_lookup_vec::Array{Float64,1}`: Profile of infection risk against distance. Entry for each one metre increment.
- `max_suscept_per_cell::Array{Float64,1}`: Maximum susceptibility per cell of the grid.
- `max_trans_per_cell::Array{Float64,1}`: Maximum transmissibility per cell of the grid.

Outputs:
- `max_cell_to_cell_transmission_rate::Array{Float64,2}`: Array of ccell-to-cell transmission rates.

Location: grid\\_fns.jl
"""
function max_cell_to_cell_transmission_rate_calc(cell_lim_array::Array{Float64,2},
                            n_holdings_per_cell::Array{Int64,1},
                            coord_type::Int64,
                            kernel_lookup_vec::Array{Float64,1},
                            max_suscept_per_cell::Array{Float64,1},
                            max_trans_per_cell::Array{Float64,1})

    #Get number of grids in use
    holding_cell_IDs_max_val = length(n_holdings_per_cell)

    #Populate max_cell_to_cell_transmission_rate
    max_cell_to_cell_transmission_rate = zeros(holding_cell_IDs_max_val,holding_cell_IDs_max_val)
    for jj = 1:holding_cell_IDs_max_val
        for ii = 1:holding_cell_IDs_max_val
            if ii==jj || n_holdings_per_cell[ii]==0 || n_holdings_per_cell[jj]==0 #Within-grid event or no holding present in grid
                max_cell_to_cell_transmission_rate[ii,jj] = Inf
            else  #Calculate distance based on type of location data
                distance = 0.0 #Initialise distance variable
                if coord_type == 1 #Cartesian co-ords (metres)
                    distance = eucl_distance_between_cells(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                            cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                            cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                            cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax
                elseif coord_type == 2 #Cartesian co-ords (km)
                    distance = eucl_distance_between_cells_convert_to_metres(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                                cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                                cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                                cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax
                elseif coord_type == 3 #Lat/Long co-ords

                    #Get distance by calling great circle distance function
                    distance = great_circle_distance_between_cells(cell_lim_array[ii,1],cell_lim_array[ii,2], #Cell_A_xMin, Cell_A_xMax
                                                            cell_lim_array[ii,3],cell_lim_array[ii,4], #Cell_A_yMin, Cell_A_yMax
                                                            cell_lim_array[jj,1],cell_lim_array[jj,2], #Cell_B_xMin, Cell_B_xMax
                                                            cell_lim_array[jj,3],cell_lim_array[jj,4]) #Cell_B_yMin, Cell_B_yMax
                end

                #Compute over estimated grid-to-grid transmission rate
                dist_idx = return_dist_idx_for_kernel(distance)
                max_cell_to_cell_transmission_rate[ii,jj] = max_trans_per_cell[ii]*max_suscept_per_cell[jj]*kernel_lookup_vec[dist_idx]

                if ii==1 && jj==10
                    println("cell_ii, cell_jj, distance: $ii, $jj, $distance")
                    println(kernel_lookup_vec[dist_idx])
                end

                if ii==1 && jj==2
                    println("cell_ii, cell_jj, distance: $ii, $jj, $distance")
                    println(kernel_lookup_vec[dist_idx])
                end
            end
        end
    end

    println("Check max_cell_to_cell_transmission_rate")

    return max_cell_to_cell_transmission_rate::Array{Float64,2}
end


#-------------------------------------------------------------------------------
### CONSTRUCT SUSCEPTIBLE & INFECTIOUS PREM ID & LOCATION TUPLES
#-------------------------------------------------------------------------------
"""
    construct_holding_by_cell_tuples(n_cells::Int64,
                                    holding_in_cell_IDs::Vector{Vector{Int64}},
                                    holding_status::Vector{Float64},
                                    holding_infectious_flag_vec::Array{Bool,1},
                                    holding_suscep_flag_vec::BitArray{1})

Construct susceptible and infectious holding ID & location tuples.

Inputs:
- `n_cells::Int64`: Number of cells in the entire grid.
- `holding_in_cell_IDs::Vector{Vector{Int64}}`: Vector of vectors. Vector per cell, with a list of holding IDs that reside in that cell.
- `holding_status::Vector{Float64}`: holding-level disease status. Entry per holding. Susceptible: 0. Infected: >0. Increases by 1 each day infected until culled.
- `holding_infectious_flag_vec::Array{Bool,1},`: Indicator of the infectious status of each holding.
- `holding_suscep_flag_vec::BitArray{1}`: Flag vector, denoting whether holding is susceptible (value 1) or not susceptible (value 0).

Outputs:
- `suscep_holding_by_cell_IDs::Vector{Vector{Int64}}`: Vector of vectors. Vector per cell, with a vector of susceptible holding in that cell.
- `infected_holding_by_cell_IDs::Vector{Vector{Int64}}`: Vector of vectors. Vector per cell, with a vector of infected holding in that cell.

Location: grid\\_fns.jl
"""
function construct_holding_by_cell_tuples(n_cells::Int64,
                                          holding_in_cell_IDs::Vector{Vector{Int64}},
                                          holding_status::Vector{Float64},
                                          holding_infectious_flag_vec::Array{Bool,1},
                                          holding_suscep_flag_vec::BitArray{1})
    #Initialise tuples
    suscep_holding_by_cell_IDs = Vector{Vector{Int64}}(undef,n_cells)
    infected_holding_by_cell_IDs = Vector{Vector{Int64}}(undef,n_cells)

     for cell_idx_itr = 1:n_cells #Populate tuples

         #Get holding in current cell of interest
         current_cell_holding_IDs = holding_in_cell_IDs[cell_idx_itr]::Array{Int64,1}

         #Get IDs of holding that are susceptible in each cell
         temp_logic_vec = current_cell_holding_IDs.*holding_suscep_flag_vec[current_cell_holding_IDs] #Multiply susceptible flag by the node ID
                                                                                 #Non-zero entries correspond to those nodes that are susceptible
         suscep_holding_by_cell_IDs[cell_idx_itr] = filter(!iszero,temp_logic_vec) #Remove zero entries.

         #Get IDs of holding that are infectious in each cell
         temp_logic_vec .= current_cell_holding_IDs.*holding_infectious_flag_vec[current_cell_holding_IDs] #Multiply infectious flag by the grid ID each node resides in.
                                                                                #Non-zero entries correspond to those nodes that are susceptible
                                                                                        #Non-zero entries correspond to those nodes that are infectious
         infected_holding_by_cell_IDs[cell_idx_itr] = filter(!iszero,temp_logic_vec) #Remove zero entries.

     end

    return suscep_holding_by_cell_IDs::Vector{Vector{Int64}},
                infected_holding_by_cell_IDs::Vector{Vector{Int64}}

end

#-------------------------------------------------------------------------------
### STATIC GRID CONSTRUCTION. PAIR OF FUNCTIONS TO COMPLETE THIS.
#-------------------------------------------------------------------------------
"""
    which_cell_fn(holding_loc_X_vals::Array{Float64,1},
                holding_loc_Y_vals::Array{Float64,1},
                bounding_box_width::Float64,
                bounding_box_height::Float64,
                n_cells_per_side::Int64)

Static grid construction: Return the grid a particular location is in.

Inputs:
- `holding_loc_X_vals::Array{Float64,1}`: East-west plane co-ordinate (per holding).
- `holding_loc_Y_vals::Array{Float64,1}`: North-south plane co-ordinate (per holding).
- `bounding_box_width::Float64`: Width of the landscape bounding box.
- `bounding_box_height::Float64`: Height of the landscape bounding box.
- `n_cells_per_side::Int64`: Grid squares need to span longest dimension of bounding box.

Outputs:
- `G::Int64`: Grid/cell index the location resides in.

Location: grid\\_fns.jl
"""
function which_cell_fn(holding_loc_X_vals::Array{Float64,1},
                        holding_loc_Y_vals::Array{Float64,1},
                        bounding_box_width::Float64,
                        bounding_box_height::Float64,
                        n_cells_per_side::Int64)

    # Get overall grid index for each holding
    #Column-major order!
    G = (floor.(holding_loc_X_vals.*n_cells_per_side/bounding_box_width)).*n_cells_per_side .+
            floor.(holding_loc_Y_vals.*n_cells_per_side/bounding_box_height) .+ 1

    return convert(Array{Int64},G)  #Ensure returned values are integers rather than floats
end

"""
    static_grid_construction(n_cells_per_side::Int64,
                            bounding_box_var::Array{Float64,1},
                            holding_loc_X_vals::Array{Float64,1},
                            holding_loc_Y_vals::Array{Float64,1})

Run static grid construction.

Inputs:
- `n_cells_per_side::Int64`: Grid squares need to span longest dimension of bounding box.
- `bounding_box_var::Array{Float64,1}`: Vector with limits for landscape bounding box. Entries [Min_x,Max_x,Min_y,Max_y].
- `holding_loc_X_vals::Array{Float64,1}`: East-west plane co-ordinate (per holding).
- `holding_loc_Y_vals::Array{Float64,1}`: North-south plane co-ordinate (per holding).

Outputs:
- `holding_cell_IDs::Array{Int64,1}`: For each holding, the ID of the cell it resides in within the final grid configuration.
- `cell_xMin::Array{Float64,1}`,`cell_yMin::Array{Float64,1}`,`cell_xMax::Array{Float64,1}`,`cell_yMax::Array{Float64,1}`: Limits for each grid cell (stored separately horizontal and vertical dimensions).

Location: grid\\_fns.jl
"""
function static_grid_construction(n_cells_per_side::Int64,
                                    bounding_box_var::Array{Float64,1},
                                    holding_loc_X_vals::Array{Float64,1},
                                    holding_loc_Y_vals::Array{Float64,1})

    #Get lengths of bounding box
    bounding_box_width = bounding_box_var[2] - bounding_box_var[1]
    bounding_box_height = bounding_box_var[4] - bounding_box_var[3]

    #Construct n_cells_per_side^2 rectangular cells, with side longest_landscape_edge/n_cells_per_side
    #Calculates which grid a particular location is in
    holding_cell_IDs = which_cell_fn(holding_loc_X_vals,holding_loc_Y_vals,bounding_box_width,bounding_box_height,n_cells_per_side)

    println("CHECK OUTPUTS OF which_cell_fn!")
    println("holding_cell_IDs: $holding_cell_IDs")

    #Get all possible max/min grid boundary values
    X_cell_edge_vals = bounding_box_var[1]:bounding_box_width/n_cells_per_side:bounding_box_var[2]
    Y_cell_edge_vals = bounding_box_var[3]:bounding_box_height/n_cells_per_side:bounding_box_var[4]
        #bounding_box_width/n_cells_per_side: Length of grid side in horizontal direction
        #bounding_box_height/n_cells_per_side: Length of grid side in vertical direction

    cell_edge_vals_xMin = X_cell_edge_vals[1:end-1]; cell_edge_vals_xMax = X_cell_edge_vals[2:end];
    cell_edge_vals_yMin  = Y_cell_edge_vals[1:end-1]; cell_edge_vals_yMax = Y_cell_edge_vals[2:end];

    #Iterate over each grid.
    #Assign cell_xMin, cell_xMax, cell_yMin, cell_yMax
    cell_xMin = zeros(maximum(holding_cell_IDs)); cell_xMax = zeros(maximum(holding_cell_IDs));
    cell_yMin = zeros(maximum(holding_cell_IDs)); cell_yMax = zeros(maximum(holding_cell_IDs));
    for ii = 1:maximum(holding_cell_IDs)
        Xgrid = (floor(Int64,(ii-1)/n_cells_per_side)) + 1
        Ygrid = mod(ii,n_cells_per_side)

        if Ygrid == 0 #Deal with case where modulo fn sets value to 0.
            Ygrid = n_cells_per_side
        end

        #Assign x boundary values to storage vectors
        cell_xMin[ii] = cell_edge_vals_xMin[Xgrid]
        cell_xMax[ii] = cell_edge_vals_xMax[Xgrid]

        #Assign y boundary values to storage vectors
        cell_yMin[ii] = cell_edge_vals_yMin[Ygrid]
        cell_yMax[ii] = cell_edge_vals_yMax[Ygrid]
    end

    return holding_cell_IDs::Array{Int64,1},
            cell_xMin::Array{Float64,1}, 
            cell_xMax::Array{Float64,1}, 
            cell_yMin::Array{Float64,1}, 
            cell_yMax::Array{Float64,1}
end
