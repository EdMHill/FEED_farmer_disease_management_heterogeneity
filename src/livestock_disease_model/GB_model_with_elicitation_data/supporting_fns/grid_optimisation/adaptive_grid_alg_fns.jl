#=
Purpose:
Functions to perform adaptive grid construction algorithm

Described in Sellman et al (2018) "Need for sped; An optimized gridding approach for spatially
explicit disease simulations" Plos Comp Bio.
=#

#-------------------------------------------------------------------------------
# NON-PRIMARY FUNCTIONS THAT ARE CALLED WITHIN THE CORE ADAPTIVE GRID ALGORITHM FN
#-------------------------------------------------------------------------------
"""
    subdivide_quad(parent_cell_limits::Array{Float64,1})

Adaptive grid algorithm function: Split index grid into four, equal grids

Inputs:
- `parent_cell_limits::Array{Float64,1}`: Vector with entries [Xmin,Xmax,Ymin,Ymax].

Outputs:
- `child_cell_limits::Array{Float64,2}`: Array with four rows. Per row, boundary limits [Xmin,Xmax,Ymin,Ymax] for a cell quadrant.

Location: adaptive//_grid//_alg//_fns.jl
"""
function subdivide_quad(parent_cell_limits::Array{Float64,1})

    #Disaggregate parent_cell_limits
    parent_cell_Xmin = parent_cell_limits[1]::Float64
    parent_cell_Xmax = parent_cell_limits[2]::Float64
    parent_cell_Ymin = parent_cell_limits[3]::Float64
    parent_cell_Ymax = parent_cell_limits[4]::Float64

    #Check size of cell
    width = parent_cell_Xmax - parent_cell_Xmin
    height = parent_cell_Ymax - parent_cell_Ymin

    #Tell a cell to subdivide into 4 smaller cells if number of nodes
    #is above a certain value.
    #Only allowed for quadratic cells. Throw error if check failed.
    if (abs(width - height) > 1e-8 )  #Throw error if discrepancy exceeds floating point precision
        error("Attempt to subdivide a cell that is not quadratic into 4 children.
                    Height: $height. Width: $width. Difference: $(abs(width - height)).")
    end

    #Get length of subdivided grids
    half_side = width * 0.5

    #Generate new grids
    child_cell_BL = [parent_cell_Xmin parent_cell_Xmin+half_side parent_cell_Ymin parent_cell_Ymin+half_side] #Bottom-left quadrant
    child_cell_TL = [parent_cell_Xmin parent_cell_Xmin+half_side parent_cell_Ymin+half_side parent_cell_Ymax] #Top-left quadrant
    child_cell_BR = [parent_cell_Xmin+half_side parent_cell_Xmax parent_cell_Ymin parent_cell_Ymin+half_side] #Bottom-right quadrant
    child_cell_TR = [parent_cell_Xmin+half_side parent_cell_Xmax parent_cell_Ymin+half_side parent_cell_Ymax] #Top-right quadrant

    #Assign newly generated child cells to output variable
    child_cell_limits = [child_cell_BL;child_cell_TL;child_cell_BR;child_cell_TR]

    return child_cell_limits::Array{Float64,2}
end

"""
    nodes_within_check(cell_limits,
                        prem_X_loc,
                        prem_Y_loc)

Adaptive grid algorithm function: Check which nodes are within the given grid cell.

Inputs:
- `cell_limits`: Vector with cell boundary limits [Xmin,Xmax,Ymin,Ymax].
- `prem_X_loc`,`prem_Y_loc`: Location information for each premises.

Outputs:
- `in_cell_flag::BitArray{1}`: Indicator vector. Entry per premises. 1 = in cell; 0 = not in cell.

Location: adaptive//_grid//_alg//_fns.jl
"""
function nodes_within_check(cell_limits,
                            prem_X_loc,
                            prem_Y_loc)

    #Disaggregate cell_limits
    cell_Xmin = cell_limits[1]::Float64
    cell_Xmax = cell_limits[2]::Float64
    cell_Ymin = cell_limits[3]::Float64
    cell_Ymax = cell_limits[4]::Float64

    #Find premises that lie within xLimits & yLimits
    in_cell_flag_x = cell_Xmin .<= prem_X_loc .< cell_Xmax
    in_cell_flag_y = cell_Ymin .<= prem_Y_loc .< cell_Ymax

    #Multiply two flag variables. Will return one if both conditions satisfied, zero otherwise
    in_cell_flag = in_cell_flag_x.*in_cell_flag_y

    return in_cell_flag::BitArray{1}
end

"""
    nodes_within_count(cell_limits,
                        prem_X_loc,
                        prem_Y_loc)

Adaptive grid algorithm function: Count number of nodes within the given grid cell.

Inputs:
- `cell_limits`: Vector with cell boundary limits [Xmin,Xmax,Ymin,Ymax].
- `prem_X_loc`,`prem_Y_loc`: Location information for each premises.

Outputs:
- `node_count_in_cell::Int64`: Number of nodes in given cell.

Location: adaptive//_grid//_alg//_fns.jl
"""
function nodes_within_count(cell_limits,
                            prem_X_loc,
                            prem_Y_loc)

    in_cell_flag = nodes_within_check(cell_limits,
                                    prem_X_loc,
                                    prem_Y_loc)

    #Get number of nodes that reside in the cell
    node_count_in_cell = sum(in_cell_flag)

    return node_count_in_cell::Int64
end

"""
    nodes_within_get_holding_ID(cell_limits,
                        prem_X_loc,
                        prem_Y_loc)

Adaptive grid algorithm function: Count number of nodes within the given grid cell.

Inputs:
- `cell_limits`: Vector with cell boundary limits [Xmin,Xmax,Ymin,Ymax].
- `prem_X_loc`,`prem_Y_loc`: Location information for each premises.

Outputs:
- `holding_in_cell_idx::Array{Int64,1}`: IDs of premises that reside in that cell.

Location: adaptive//_grid//_alg//_fns.jl
"""
function nodes_within_get_holding_ID(cell_limits,
                                prem_X_loc,
                                prem_Y_loc)

    in_cell_flag = nodes_within_check(cell_limits,
                                    prem_X_loc,
                                    prem_Y_loc)

    #Get the premises ID that reside in that cell
    holding_in_cell_idx = findall(in_cell_flag .== 1)

    return holding_in_cell_idx::Array{Int64,1}
end


#-------------------------------------------------------------------------------
# PRIMARY FUNCTION
#-------------------------------------------------------------------------------

"""
    adaptive_grid_construct_fn(avg_node_per_grid_target::Float64,
                            bounding_box_var::Array{Float64,1},
                            longest_landscape_edge::Float64,
                            holding_loc_X_vals::Array{Float64,1},
                            holding_loc_Y_vals::Array{Float64,1})

Construct grid over the landscape using adaptive algorithm with the final grid configuration having (on average) avg_node_per_grid_target nodes per grid cell.

Inputs:
- `avg_node_per_grid_target::Float64`: The algorithm aims for final grid configuration to get (on average) this number of nodes per grid cell.
- `bounding_box_var::Array{Float64,1}`: Vector [x_min, x_max, y_min, y_max] for landscape box.
- `longest_landscape_edge::Float64`: Dimension of landscape of highest value. Used to construct initial parent cell square.
- `holding_loc_X_vals::Array{Float64,1}`: Premises location data, x co-ordinate.
- `holding_loc_Y_vals::Array{Float64,1}`: Premises location data, y co-ordinate.

Outputs:
- `holding_cell_IDs::Array{Int64,1}`: For each premises, the ID of the cell it resides in within the final grid configuration.
- `grid_xMin::Array{Float64,1}`: For each grid, left edge value.
- `grid_xMax::Array{Float64,1}`: For each grid, right edge value.
- `grid_yMin::Array{Float64,1}`: For each grid, bottom edge value.
- `grid_yMax::Array{Float64,1}`: For each grid, top edge value.

Location: adaptive//_grid//_alg//_fns.jl
"""
function  adaptive_grid_construct_fn(avg_node_per_grid_target::Float64,
                                    bounding_box_var::Array{Float64,1},
                                    longest_landscape_edge::Float64,
                                    holding_loc_X_vals::Array{Float64,1},
                                    holding_loc_Y_vals::Array{Float64,1})
#Principle. List of grids to test subdivide.
#Iterate through. For those accepted, add to list to be tested in next generation
#Keep subdividing until get to granularity that minimises the squared difference between log-number of nodes per cell
#and the specified threshold log-number of nodes per cell.

    #Get total number of holdings in landscape
    n_holdings = length(holding_loc_X_vals)

    #Take landscsape bounds.
    #Construct square from bottom left corner. Side length of longest_landscape_edge
    #Set up initial grid square. Put into a 2D array
    parent_cell_square = [bounding_box_var[1] bounding_box_var[1]+longest_landscape_edge bounding_box_var[3] bounding_box_var[3]+longest_landscape_edge]

    #Initialise variables
    cell_subdivide_check_val = [1] #Initialise as a 1D array. Ensure logical indexing works during first iteration!
    parent_cell_node_count = [n_holdings] #Initialise as a 1D array. Ensure logical indexing works during first iteration!

    #Initialise conidition variable. Will be zero when algorithm can stop
    subdivide_check_flag = 1

    #Enter algorithm loop
    while subdivide_check_flag > 0

        #Assign logic vectors to variables
        grid_unchecked_logic_vec = cell_subdivide_check_val.==0
        grid_checked_logic_vec = cell_subdivide_check_val.==1

        #Get those grid cells that do not need checking for subdivision, initialise revised variables
        parent_cell_square_unchecked = parent_cell_square[grid_unchecked_logic_vec,:]
        revised_parent_cell_square = parent_cell_square_unchecked #Initialise parent grid cells variable
        revised_cell_subdivide_check_val = zeros(Int64,size(parent_cell_square_unchecked,1))
        revised_parent_cell_node_count = parent_cell_node_count[grid_unchecked_logic_vec]

        #Get boundary limits of parent cells for which subdivision test will be performed
        parent_cell_to_check = parent_cell_square[grid_checked_logic_vec,:]

        #Get ID of parent_cell_square to be checked
        parent_cell_to_check_IDs = findall(grid_checked_logic_vec)

        #For parent cells to be tested, get "closeness" to target number of nodes per grid
        parent_cell_closeness = (log.(parent_cell_node_count[grid_checked_logic_vec]) .- log.(avg_node_per_grid_target)).^2

        #Iterate over grids to test subdivision
        for ii = 1:length(parent_cell_to_check_IDs)

            #Get ID of parent cell under consideration in present iteration
            parent_cell_ID = parent_cell_to_check_IDs[ii]

            #Input all premises. Check which lie within the current cell of interest
            in_cell_holding_ID = nodes_within_get_holding_ID(parent_cell_to_check[ii,:],
                                                 holding_loc_X_vals,
                                                 holding_loc_Y_vals)

            #Get premises located within current parent cell
            parent_holding_X_loc = holding_loc_X_vals[in_cell_holding_ID]
            parent_holding_Y_loc = holding_loc_Y_vals[in_cell_holding_ID]

            #For current grid, subdivide into four grids
            child_cell_limits = subdivide_quad(parent_cell_to_check[ii,:]::Array{Float64,1})

            #Get counts per child grid
            node_per_child_vec = zeros(Int64,4)
            for child_grid_idx = 1:4
                node_per_child_vec[child_grid_idx] = nodes_within_count(child_cell_limits[child_grid_idx,:],
                                                                        parent_holding_X_loc,
                                                                        parent_holding_Y_loc)
            end

            #Assign logic vector of child cells containing nodes to variable
            child_with_nodes_logic_vec = node_per_child_vec.>0

            #Number of child grids actually containing nodes
            child_with_node_num = sum(child_with_nodes_logic_vec)

            #Calculate new and original variance from the target value
            parent_var = parent_cell_closeness[ii] #Parent cell term
            child_var = (sum((log.(node_per_child_vec[child_with_nodes_logic_vec]) .- log(avg_node_per_grid_target)).^2))/child_with_node_num #Child cell term
            if  child_var < parent_var
                #Accepted, move forward with subdivide retained
                #Only keep those gird cells containing nodes
                current_cell_grid_next_gen = child_cell_limits[child_with_nodes_logic_vec,:]
                current_cell_node_count = node_per_child_vec[child_with_nodes_logic_vec]

                #Set value to be appended to cell_subdivide_check_val
                current_cell_subdivide_check_val = ones(Int64,child_with_node_num)

                #Construct revised parent_cell_square & cell_subdivide_check_val lists
                revised_parent_cell_square = [revised_parent_cell_square;current_cell_grid_next_gen]
                revised_cell_subdivide_check_val = [revised_cell_subdivide_check_val;current_cell_subdivide_check_val]
                revised_parent_cell_node_count = [revised_parent_cell_node_count;current_cell_node_count]
            else
                #Not accepted. Retain parent cell.
                current_cell_grid_next_gen = parent_cell_to_check[ii,:]
                current_cell_node_count = parent_cell_node_count[parent_cell_ID]

                #Set value to be appended to cell_subdivide_check_val
                current_cell_subdivide_check_val = 0

                #Construct revised parent_cell_square & cell_subdivide_check_val lists
                revised_parent_cell_square = vcat(revised_parent_cell_square,current_cell_grid_next_gen') #Transpose current_cell_grid_next_gen so it is a row vector
                revised_cell_subdivide_check_val = [revised_cell_subdivide_check_val;current_cell_subdivide_check_val]
                revised_parent_cell_node_count = [revised_parent_cell_node_count;current_cell_node_count]
            end
        end

        #End of generation, set parent_cell_square, cell_subdivide_check_val, parent_cell_node_count
        parent_cell_square = revised_parent_cell_square
        cell_subdivide_check_val = revised_cell_subdivide_check_val
        parent_cell_node_count = revised_parent_cell_node_count

        #Check if there remains any parent cells that need to be checked for subdivision
        subdivide_check_flag = sum(cell_subdivide_check_val) #Set condition variable to 0 to leave while loop
    end

    #FInal configuration obtained
    complete_cell_config = parent_cell_square

    #For each premises, get the cell ID in which the premises resides
    complete_cell_Num = size(complete_cell_config,1) #Number of cells in final config, matches row count in parent_cell_square
    holding_cell_IDs = zeros(Int64,n_holdings) #Initialise storage vector. Ensure length matches number of premises
    for cell_IdxItr = 1:complete_cell_Num
       #Input all premises. Check which lie within the current cell of interest
       in_cell_holding_ID = nodes_within_get_holding_ID(complete_cell_config[cell_IdxItr,:],
                                            holding_loc_X_vals,
                                            holding_loc_Y_vals)
        #Assign the grid ID
        holding_cell_IDs[in_cell_holding_ID] .= cell_IdxItr
    end

    #Assign cell boundary limits to variables
    grid_xMin = complete_cell_config[:,1]
    grid_xMax = complete_cell_config[:,2]
    grid_yMin = complete_cell_config[:,3]
    grid_yMax = complete_cell_config[:,4]

    return holding_cell_IDs::Array{Int64,1},
            grid_xMin::Array{Float64,1}, grid_xMax::Array{Float64,1}, grid_yMin::Array{Float64,1}, grid_yMax::Array{Float64,1}
end
