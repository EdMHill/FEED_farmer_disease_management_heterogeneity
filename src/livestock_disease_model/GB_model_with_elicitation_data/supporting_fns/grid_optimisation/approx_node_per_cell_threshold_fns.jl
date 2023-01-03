#=
Purpose:
Functions to approximate the optimum number of nodes per cell to optimise the efficiency
of the spatial simulation

Described in Sellman et al (2018) "Need for sped; An optimized gridding approach for spatially
explicit disease simulations" Plos Comp Bio.
=#

#-------------------------------------------------------------------------------
#COMPUTE DISTANCES BETWEEN PAIRS OF GRID CELLS
#-------------------------------------------------------------------------------
"""
    get_distance_between_cell_pairs(A_x::Int64,A_y::Int64,
                                    B_x::Int64,B_y::Int64,
                                    longest_side::Float64,
                                    grid_resolution::Int64,
                                    coord_type::Int64)

Function used in approximation of optimum number of nodes per cell: Compute distances between pairs of grid cells.

Inputs:
- `A_x::Int64,A_y::Int64`: Position for cell A.
- `B_x::Int64,B_y::Int64`: Position for cell B.
- `longest_side::Float64`: Length of landscape.
- `grid_resolution::Int64`: Number of cells per side.
- `coord_type::Int64`: Value 1 ("Cartesian", metres), 2 ("Cartesian", km) or 3 ("LatLong").

Outputs:
- `cell_pair_dist::Float64`: Threshold number of nodes per cell. To be used in adaptive grid configuration construction.

Location: approx\\_node\\_per\\_cell\\_threshold\\_fns.jl
"""
function get_distance_between_cell_pairs(A_x::Int64,A_y::Int64,
                                            B_x::Int64,B_y::Int64,
                                            longest_side::Float64,
                                            grid_resolution::Int64,
                                            coord_type::Int64)

    # Set up flag values for magnitude of x,y position of cell B relative to
    # cell A
    a = 0
    b = 0

    if (A_x == B_x)
        a = 0
    elseif (A_x < B_x)
        a = 1
    elseif (A_x > B_x)
        a = -1
    end

    if (A_y == B_y)
        b = 0
    elseif (A_y < B_y)
        b = 1
    elseif (A_y > B_y)
        b = -1
    end

    #Compute distance based on co-ordinate system in use
    if coord_type == 1 #Cartesian co-ordinates (units of metres)
        t1 = ((A_x + a)*longest_side - B_x * longest_side) / grid_resolution
        t2 = ((A_y + b)*longest_side - B_y * longest_side) / grid_resolution
        d_sq = t1.^2 + t2.^2
        cell_pair_dist = sqrt(d_sq)
    elseif coord_type == 2 #Cartesian co-ordinates (units of km, so multiply output by 1000 so output returned is in metres)
        t1 = ((A_x + a)*longest_side - B_x * longest_side) / grid_resolution
        t2= ((A_y + b)*longest_side - B_y * longest_side) / grid_resolution
        d_sq = t1.^2 + t2.^2
        cell_pair_dist = (sqrt(d_sq))*1000
    elseif coord_type == 3 #LatLong co-ordinates
        cell_length::Float64 = longest_side/grid_resolution
        lat1 = (A_y + b)*cell_length; lat2 = B_y*cell_length;
        lon1 = (A_x + a)*cell_length; lon2 = B_x*cell_length;

        cell_pair_dist = great_circle_distance(lat1, lon1, lat2, lon2) #Will return distance in metres
    end

    return cell_pair_dist::Float64
end

#-------------------------------------------------------------------------------
# PRIMARY FUNCTION
#-------------------------------------------------------------------------------
"""
    approx_node_per_cell_threshold_fn(L::LandscapeData,
                                        kernel_lookup_vec::Array{Float64,1},
                                        delta_t::Float64,
                                        approx_node_per_cell_grid_size_threshold::Int64)

Approximate the optimum number of nodes per cell, for use in adaptive grid construction algorithm.

Inputs:
- `L::LandscapeData`: For landscape, contains boundary information, maxium transmissibility and median susceptibility.
- `kernel_lookup_vec::Array{Float64,1}`: Profile to specify risk of infection with respect to distance.
- `delta_t::Float64`: Timestep length
- `approx_node_per_cell_grid_size_threshold::Int64`: Maximum number of cells per grid side to test.

Outputs:
- `approx_node_per_cell_threshold::Float64`: Threshold number of nodes per cell. To be used in adaptive grid configuration construction.

Location: approx\\_node\\_per\\_cell\\_threshold\\_fns.jl
"""
function approx_node_per_cell_threshold_fn(L::LandscapeData,
                                        kernel_lookup_vec::Array{Float64,1},
                                        delta_t::Float64,
                                        approx_node_per_cell_grid_size_threshold::Int64)

    #Set up cell resolutions to be tested
    k_to_test = 1:1:approx_node_per_cell_grid_size_threshold

    #Initialise variables used in main loop
    k_element = 1 #Increment variable
    k_calls_vector = zeros(length(k_to_test))
    avg_k_calls_vector = zeros(length(k_to_test))
    nodes_per_cell_vec = zeros(length(k_to_test))

    #Iterate over the different grid configuration
    for grid_resolution_itr = 1:length(k_to_test)

        #Get number of cells per side of grid
        k::Int64 = k_to_test[grid_resolution_itr]

        #Find number of nodes per cell, assign to storage vector
        nodes_per_cell = L.get_n_nodes/(k*k)
        nodes_per_cell_vec[grid_resolution_itr] = nodes_per_cell

        k_calls_this_k = 0 #number of kernel calls in total for this grid conf.
        if (k == 1) #No gridding (one large cell). Only pairwise.
            k_calls_this_k = k_calls_this_k + (nodes_per_cell - 1)
        else #More than one cell present

            #To save time we only calculate the number of kernel calls for one quadrant of the cells to
            #all other cells and multiply with 4. If k is uneven we do that for a rectangle and add on
            #the number of kernel calls for the middle cell to all other cells at the end.

            #Get indexes for each grid based on value of k
            x_val_all = 1:1:k #x dimension
            y_val_all = 1:1:k #y dimension

            #Get coordinates for cells in first quadrant
            if mod(k,2) == 0  #when k is even
                quadrant_width = convert(Int64,k/2)
                quadrant_height = convert(Int64,k/2)
            else #when k is odd
                quadrant_width = convert(Int64,(k+1)/2)
                quadrant_height = convert(Int64,(k-1)/2)
            end

            x_val_subset = 1:1:quadrant_width
            y_val_subset = 1:1:quadrant_height

            #Loop over all pairs of cells (first cell A x and y, and then cell B x and y).
            for A_x_itr = 1:length(x_val_subset) #Cell A_x coordinates.
                for A_y_itr = 1:length(y_val_subset) #Cell A_y coordinates.
                    A_x = x_val_subset[A_x_itr]
                    A_y = y_val_subset[A_y_itr]
                    #For this pair of Ax and Ay do all B cells (including A itself).
                    k_calls_A = (nodes_per_cell - 1) #Pairwise within the cell itself.
                    for B_x_itr = 1:length(x_val_all) #Cell B_x coordinates.
                        for B_y_itr = 1:length(y_val_all) #Cell B_x coordinates.
                            B_x = x_val_all[B_x_itr]
                            B_y = y_val_all[B_y_itr]
                            if (!(A_x == B_x && A_y == B_y)) #Only if A and B are not the same cell.
                                d = get_distance_between_cell_pairs(A_x, A_y, B_x, B_y, L.longest_side, k, L.coord_type)
                                distIdx = return_dist_idx_for_kernel(d)
                                kernel_value = kernel_lookup_vec[distIdx]
                                exponent = -L.get_trans * L.get_susc * kernel_value * delta_t
                                #                               double exponent = -L.get_trans() * L.get_susc() * distance_kernel(d);
                                p_over = oneMinusExp(exponent)
                                #                                p_over = 1 - std::pow(1 - p_over, nodes_per_cell);
                                k_calls_this_pair = (nodes_per_cell * p_over) + 1 #Expected n kernel calls for Cell(Ax, Ay)->Cell(Bx, By) plus one for the calculation of p_over.
                                k_calls_A += k_calls_this_pair  #Add the kernel calls between cells A and B to A's total.
                            end
                        end
                    end
                    k_calls_this_k += 4*k_calls_A; #Add A's total kernel calls to the grand total (cells equivalent to A occur 4 times in the landscape).
                end
            end

            #When all pairs have been looped over and k is odd there still remains the central cell in the configuration.
            if (mod(k,2) == 1)
                #Get grid index for central grid
                A_x = convert(Int64,ceil(k/2))
                A_y = convert(Int64,ceil(k/2))
                k_calls_A = (nodes_per_cell - 1) #Pairwise within the cell itself.
                for B_x_itr = 1:length(x_val_all) #Cell B_x coordinates.
                    for B_y_itr = 1:length(y_val_all) #Cell B_x coordinates.
                        B_x = x_val_all[B_x_itr]
                        B_y = y_val_all[B_y_itr]
                        if (!(A_x == B_x && A_y == B_y)) #Only if A and B are not the same cell.
                            d = get_distance_between_cell_pairs(A_x, A_y, B_x, B_y, L.longest_side, k, L.coord_type)
                            distIdx = return_dist_idx_for_kernel(d)
                            kernel_value = kernel_lookup_vec[distIdx]
                            exponent = -L.get_trans * L.get_susc * kernel_value * delta_t
                            #                               double exponent = -L.get_trans() * L.get_susc() * distance_kernel(d);
                            p_over = oneMinusExp(exponent)
                            #                                p_over = 1 - std::pow(1 - p_over, nodes_per_cell);
                            k_calls_this_pair = (nodes_per_cell * p_over) + 1 #Expected n kernel calls for Cell(Ax, Ay)->Cell(Bx, By) plus one for the calculation of p_over.
                            k_calls_A += k_calls_this_pair  #Add the kernel calls between cells A and B to A's total
                        end
                    end
                end
                k_calls_this_k += k_calls_A
            end
        end

        #Assign kernel calls for this k value to storage vector
        k_calls_vector[k_element] = k_calls_this_k

        #Assign average number of expected kernel function call per cell to storage vector
        #Total number of cells is k*k
        avg_k_calls_vector[k_element] = k_calls_this_k/(k*k)

        #Increment k_element
        k_element = k_element + 1
    end

    #Find configuration that returned the lowest number of kernel calls
    k_min_idx::Int64 = argmin(avg_k_calls_vector)
    println("nodes_per_cell_vec: $nodes_per_cell_vec")
    println("k_calls_vector: $k_calls_vector")
    println("avg_k_calls_vector: $avg_k_calls_vector")
    println("k_min_idx: $k_min_idx")

    #Get expected number of nodes per cell, given selected configuration
    approx_node_per_cell_threshold = nodes_per_cell_vec[k_min_idx]

    return approx_node_per_cell_threshold::Float64

end
