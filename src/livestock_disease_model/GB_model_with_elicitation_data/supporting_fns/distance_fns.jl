#=
Purpose:
File to house functions to compute distances between point locations
Used in spatial disease outbreak simulations

Formats:
 - Euclidean distance squared (node-to-node)
 - Euclidean distance (node-to-node)
 - Euclidean distance squared (cell-to-cell)
 - Euclidean distance (cell-to-cell)
 - Great circle distance (node-to-node)
 - Great circle distance (cell-to-cell)

File also contains the function 'return_dist_idx_for_kernel', returns a distance value
used in tranmission kernel lookup vector.
=#

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE SQUARED (BETWEEN NODES)
#-------------------------------------------------------------------------------
"""
    sq_distance(Loc_A_xVal::Float64,
                Loc_A_yVal::Float64,
                Loc_B_xVal::Float64,
                Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function sq_distance(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    delta_x::Float64 = Loc_A_xVal - Loc_B_xVal #Find difference in east-west plane
    delta_y::Float64 = Loc_A_yVal - Loc_B_yVal #Find difference in south-north plane

    return (delta_x*delta_x + delta_y*delta_y)::Float64
end

"""
    sq_distance_convert_to_km(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function sq_distance_convert_to_km(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    return ((sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*0.001)::Float64
end

"""
    sq_distance_convert_to_metres(Loc_A_xVal::Float64,
                                Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64,
                                Loc_B_yVal::Float64)

Return squared euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: distance_fns.jl
"""
function sq_distance_convert_to_metres(Loc_A_xVal::Float64,
                            Loc_A_yVal::Float64,
                            Loc_B_xVal::Float64,
                            Loc_B_yVal::Float64)

    return ((sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*1000.)::Float64
end

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE (BETWEEN NODES)
#-------------------------------------------------------------------------------
"""
    eucl_distance(Loc_A_xVal::Float64,
                Loc_A_yVal::Float64,
                Loc_B_xVal::Float64,
                Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function eucl_distance(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    # Take square root of squared distance between nodes A & B
    return (sqrt(sq_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal)))::Float64

end


"""
    eucl_distance_convert_to_km(Loc_A_xVal::Float64,
                                Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64,
                                Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function eucl_distance_convert_to_km(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    return ((eucl_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*0.001)::Float64

end

"""
    eucl_distance_convert_to_metres(Loc_A_xVal::Float64,
                                    Loc_A_yVal::Float64,
                                    Loc_B_xVal::Float64,
                                    Loc_B_yVal::Float64)

Return euclidean distance between two nodes, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: distance_fns.jl
"""
function eucl_distance_convert_to_metres(Loc_A_xVal::Float64, Loc_A_yVal::Float64,
                                Loc_B_xVal::Float64, Loc_B_yVal::Float64)


    return ((eucl_distance(Loc_A_xVal, Loc_A_yVal,
                                    Loc_B_xVal, Loc_B_yVal))*1000.)::Float64

end

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE SQUARED (BETWEEN CELLS)
#-------------------------------------------------------------------------------
"""
    sq_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                            cell_A_yMin::Float64, cell_A_yMax::Float64,
                            cell_B_xMin::Float64, cell_B_xMax::Float64,
                            cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest squared euclidean distance between two cells, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function sq_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    #If an edge of cell_B overlaps with an edge in cell_A, the two points
    #between which the shortest distance is found is on a straight horizontal or
    #vertical line.

    #Otherwise calculate the euclidean distance.

    if ( ((cell_A_yMin >= cell_B_yMin) && (cell_A_yMax <= cell_B_yMax)) || ((cell_B_yMin >= cell_A_yMin) && (cell_B_yMax <= cell_A_yMax)) )
        #Check for horizontal alignment first. i.e. Shortest distance is found on a straight horizontal line

        #Get shortest horizontal distance between cells A & B
        distances_x = abs.([cell_A_xMin - cell_B_xMin,
                            cell_A_xMin - cell_B_xMax,
                            cell_A_xMax - cell_B_xMin,
                            cell_A_xMax - cell_B_xMax])

        delta_x::Float64 = minimum(distances_x)

        return (delta_x*delta_x)::Float64

    elseif ( ((cell_A_xMin >= cell_B_xMin) && (cell_A_xMax <= cell_B_xMax)) || ((cell_B_xMin >= cell_A_xMin) && (cell_B_xMax <= cell_A_xMax)) )
        #Check for vertical alignment. i.e. Shortest distance is found on a straight vertical line

        #Get shortest vertical distance between cells A & B
        distances_y = abs.([cell_A_yMin - cell_B_yMin,
                            cell_A_yMin - cell_B_yMax,
                            cell_A_yMax - cell_B_yMin,
                            cell_A_yMax - cell_B_yMax])

        delta_y::Float64 = minimum(distances_y)    

        return (delta_y*delta_y)::Float64

    else
        ##No vertical or horizontal overlap, calculate euclidean distance


        #Get shortest horizontal distance between cells A & B
        #Ensure all distances are positive, use absolute value fn
        distances_x = abs.([cell_A_xMin - cell_B_xMin,  #cell_A.xVal_min - cell_B.xVal_min
                            cell_A_xMin - cell_B_xMax,  #cell_A.xVal_min - cell_B.xVal_max
                            cell_A_xMax - cell_B_xMin,  #cell_A.xVal_max - cell_B.xVal_min
                            cell_A_xMax - cell_B_xMax])  #cell_A.xVal_max - cell_B.xVal_max

        delta_x = minimum(distances_x)

        #Get shortest vertical distance between cells A & B
        #Ensure all distances are positive, use absolute value fn
        distances_y = abs.([cell_A_yMin - cell_B_yMin,  #cell_A.yVal_min - cell_B.yVal_min
                            cell_A_yMin - cell_B_yMax,  #cell_A.yVal_min - cell_B.yVal_max
                            cell_A_yMax - cell_B_yMin,  #cell_A.yVal_max - cell_B.yVal_min
                            cell_A_yMax - cell_B_yMax])  #cell_A.yVal_max - cell_B.yVal_max

        delta_y = minimum(distances_y)

        return (delta_x*delta_x + delta_y*delta_y)::Float64
    end

end


"""
    sq_distance_between_cells_convert_to_km(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                        cell_A_yMin::Float64, cell_A_yMax::Float64,
                                        cell_B_xMin::Float64, cell_B_xMax::Float64,
                                        cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest squared euclidean distance between two cells, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function sq_distance_between_cells_convert_to_km(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    # Take square root of squared distance between cells A & B
    return ((sq_distance_between_cells(cell_A_xMin, cell_A_xMax,
                                        cell_A_yMin, cell_A_yMax,
                                        cell_B_xMin, cell_B_xMax,
                                        cell_B_yMin, cell_B_yMax))*0.001)::Float64

end

"""
    sq_distance_between_cells_convert_to_metres(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                                cell_A_yMin::Float64, cell_A_yMax::Float64,
                                                cell_B_xMin::Float64, cell_B_xMax::Float64,
                                                cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest squared euclidean distance between two cells, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: distance_fns.jl
"""
function sq_distance_between_cells_convert_to_metres(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    # Take square root of squared distance between cells A & B
    return ((sq_distance_between_cells(cell_A_xMin, cell_A_xMax,
                                        cell_A_yMin, cell_A_yMax,
                                        cell_B_xMin, cell_B_xMax,
                                        cell_B_yMin, cell_B_yMax))*1000.)::Float64

end

#-------------------------------------------------------------------------------
### EUCLIDEAN DISTANCE (BETWEEN CELLS)
#-------------------------------------------------------------------------------
"""
    eucl_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                cell_A_yMin::Float64, cell_A_yMax::Float64,
                                cell_B_xMin::Float64, cell_B_xMax::Float64,
                                cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest euclidean distance between two cells, in metres.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function eucl_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    # Take square root of squared distance between cells A & B
    return (sqrt(sq_distance_between_cells(cell_A_xMin, cell_A_xMax,
                                        cell_A_yMin, cell_A_yMax,
                                        cell_B_xMin, cell_B_xMax,
                                        cell_B_yMin, cell_B_yMax)))::Float64

end

"""
    eucl_distance_between_cells_convert_to_km(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                            cell_A_yMin::Float64, cell_A_yMax::Float64,
                                            cell_B_xMin::Float64, cell_B_xMax::Float64,
                                            cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest euclidean distance between two cells, in km.\n
Used when inputs using a co-ordinate system with distances in units of metre.

Location: distance_fns.jl
"""
function eucl_distance_between_cells_convert_to_km(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    # Take square root of squared distance between cells A & B
    return ((eucl_distance_between_cells(cell_A_xMin, cell_A_xMax,
                                        cell_A_yMin, cell_A_yMax,
                                        cell_B_xMin, cell_B_xMax,
                                        cell_B_yMin, cell_B_yMax))*0.001)::Float64

end

"""
    eucl_distance_between_cells_convert_to_metres(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                                cell_A_yMin::Float64, cell_A_yMax::Float64,
                                                cell_B_xMin::Float64, cell_B_xMax::Float64,
                                                cell_B_yMin::Float64, cell_B_yMax::Float64)

Return shortest euclidean distance between two cells, in metres.\n
Used when inputs using a co-ordinate system with distances in units of km.

Location: distance_fns.jl
"""
function eucl_distance_between_cells_convert_to_metres(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                    cell_A_yMin::Float64, cell_A_yMax::Float64,
                                    cell_B_xMin::Float64, cell_B_xMax::Float64,
                                    cell_B_yMin::Float64, cell_B_yMax::Float64)

    # Take square root of squared distance between cells A & B
    return ((eucl_distance_between_cells(cell_A_xMin, cell_A_xMax,
                                        cell_A_yMin, cell_A_yMax,
                                        cell_B_xMin, cell_B_xMax,
                                        cell_B_yMin, cell_B_yMax))*1000.)::Float64

end

#-------------------------------------------------------------------------------
### GREAT CIRCLE DISTANCE (BY DEFAULT IN KM, SCALE UP TO METRES)
#-------------------------------------------------------------------------------

"""
    great_circle_distance(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)

Use the haversine formula to give the great-circle distances between two points on a sphere from their longitudes and latitudes.

Note, unit of output is metres.

Location: distance_fns.jl
"""
function great_circle_distance(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)

    #Compute haversine formula
    lat_long_dist = 2 * 6371 * asin(sqrt(sind((lat2 - lat1)*0.5)*sind((lat2 - lat1)*0.5) +
                                            cosd(lat1) * cosd(lat2) * sind((lon2 - lon1)*0.5)* sind((lon2 - lon1)*0.5)))

    #Return distance
    return lat_long_dist*1000.::Float64

end

"""
    great_circle_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                        cell_A_yMin::Float64, cell_A_yMax::Float64,
                                        cell_B_xMin::Float64, cell_B_xMax::Float64,
                                        cell_B_yMin::Float64, cell_B_yMax::Float64)

Use the haversine formula to give the shortest great-circle distances between two cells on a sphere from their longitudes and latitudes.

Note, unit of output is metres.

Location: distance_fns.jl
"""
function great_circle_distance_between_cells(cell_A_xMin::Float64, cell_A_xMax::Float64,
                                            cell_A_yMin::Float64, cell_A_yMax::Float64,
                                            cell_B_xMin::Float64, cell_B_xMax::Float64,
                                            cell_B_yMin::Float64, cell_B_yMax::Float64)

    if ( ((cell_A_yMin >= cell_B_yMin) && (cell_A_yMax <= cell_B_yMax)) || ((cell_B_yMin >= cell_A_yMin) && (cell_B_yMax <= cell_A_yMax)) )
        #Check for horizontal alignment first. i.e. Shortest distance is found on a constant latitude line

        #Get shortest horizontal distance between cells A & B
        cell_distances_x = abs.([cell_A_xMin - cell_B_xMin,
                                cell_A_xMin - cell_B_xMax,
                                cell_A_xMax - cell_B_xMin,
                                cell_A_xMax - cell_B_xMax])

        min_grid_dist_x_idx::Int64 = findmin(cell_distances_x)[2] #Second entry, index of the minimum over the given dimensions.

        #Will take value 1, 2, 3 or 4. Throw error if another value
        #Lookup table based on index value obtained
        #Assign longitude coords
        if min_grid_dist_x_idx == 1 #cell_A min, cell_B min
            lon1 = cell_A_xMin
            lon2 = cell_B_xMin
        elseif min_grid_dist_x_idx == 2 #cell_A min, cell_B max
            lon1 = cell_A_xMin
            lon2 = cell_B_xMax
        elseif min_grid_dist_x_idx == 3 #cell_A max, cell_B min
            lon1 = cell_A_xMax
            lon2 = cell_B_xMin
        elseif min_grid_dist_x_idx == 4 #cell_A max, cell_B max
            lon1 = cell_A_xMax
            lon2 = cell_B_xMax
        else #Unexpected value for min_grid_dist_x_idx, throw error
            error("min_grid_dist_x_idx has value $min_grid_dist_x_idx. min_grid_dist_x_idx must have value 1, 2, 3 or 4.")
        end

        #Get latitude midpoint of smaller grid cell
        if ((cell_A_yMin >= cell_B_yMin) && (cell_A_yMax <= cell_B_yMax)) 
            #Cell A is vertically within top and bottom boundaries of cell B
            lat1 = (cell_A_yMin + cell_A_yMax)/2
            lat2 = (cell_A_yMin + cell_A_yMax)/2
        else
            #Cell B is vertically within top and bottom boundaries of cell A
            lat1 = (cell_B_yMin + cell_B_yMax)/2
            lat2 = (cell_B_yMin + cell_B_yMax)/2
        end

        #Return great circle distance
        return (great_circle_distance(lat1, lon1, lat2, lon2))::Float64

    elseif ( ((cell_A_xMin >= cell_B_xMin) && (cell_A_xMax <= cell_B_xMax)) || ((cell_B_xMin >= cell_A_xMin) && (cell_B_xMax <= cell_A_xMax)) )
        #Check for vertical alignment. i.e. Shortest distance is found on a constant longitude line

        #Get shortest vertical distance between cells A & B
        cell_distances_y = abs.([cell_A_yMin - cell_B_yMin,
                                cell_A_yMin - cell_B_yMax,
                                cell_A_yMax - cell_B_yMin,
                                cell_A_yMax - cell_B_yMax])

        min_grid_dist_y_idx = findmin(cell_distances_y)[2] #Second entry, index of the minimum over the given dimensions.

        #Will take value 1, 2, 3 or 4. Throw error if another value
        #Lookup table based on index value obtained
        #Assign latitude coords
        if min_grid_dist_y_idx == 1 #cell_A min, cell_B min
            lat1 = cell_A_yMin
            lat2 = cell_B_yMin
        elseif min_grid_dist_y_idx == 2 #cell_A min, cell_B max
            lat1 = cell_A_yMin
            lat2 = cell_B_yMax
        elseif min_grid_dist_y_idx == 3 #cell_A max, cell_B min
            lat1 = cell_A_yMax
            lat2 = cell_B_yMin
        elseif min_grid_dist_y_idx == 4 #cell_A max, cell_B max
            lat1 = cell_A_yMax
            lat2 = cell_B_yMax
        else #Unexpected value for min_grid_dist_x_idx, throw error
            error("min_grid_dist_y_idx has value $min_grid_dist_y_idx. min_grid_dist_y_idx must have value 1, 2, 3 or 4.")
        end

        #Get longitude midpoint of smaller grid cell
        if ((cell_A_xMin >= cell_B_xMin) && (cell_A_xMax <= cell_B_xMax)) 
            #Cell A is horizontally within left and right boundaries of cell B
            lon1 = (cell_A_xMin + cell_A_xMax)/2
            lon2 = (cell_A_xMin + cell_A_xMax)/2
        else
            #Cell B is horizontally within left and right boundaries of cell A
            lon1 = (cell_B_xMin + cell_B_xMax)/2
            lon2 = (cell_B_xMin + cell_B_xMax)/2
        end

        #Return great circle distance
        return (great_circle_distance(lat1, lon1, lat2, lon2))::Float64

    else
    ##No vertical or horizontal overlap, calculate great circle distance

        #Get shortest horizontal distance between cells A & B
        #Ensure all distances are positive, use absolute value fn
        cell_distances_x = abs.([cell_A_xMin - cell_B_xMin,
                                cell_A_xMin - cell_B_xMax,
                                cell_A_xMax - cell_B_xMin,
                                cell_A_xMax - cell_B_xMax])

        min_grid_dist_x_idx = findmin(cell_distances_x)[2] #Second entry, index of the minimum over the given dimensions.

        #Will take value 1, 2, 3 or 4. Throw error if another value
        #Lookup table based on index value obtained
        #Assign longitude coords
        if min_grid_dist_x_idx == 1 #cell_A min, cell_B min
            lon1 = cell_A_xMin
            lon2 = cell_B_xMin
        elseif min_grid_dist_x_idx == 2 #cell_A min, cell_B max
            lon1 = cell_A_xMin
            lon2 = cell_B_xMax
        elseif min_grid_dist_x_idx == 3 #cell_A max, cell_B min
            lon1 = cell_A_xMax
            lon2 = cell_B_xMin
        elseif min_grid_dist_x_idx == 4 #cell_A max, cell_B max
            lon1 = cell_A_xMax
            lon2 = cell_B_xMax
        else #Unexpected value for min_grid_dist_x_idx, throw error
            error("min_grid_dist_x_idx has value $min_grid_dist_x_idx. min_grid_dist_x_idx must have value 1, 2, 3 or 4.")
        end


        #Get shortest vertical distance between cells A & B
        #Ensure all distances are positive, use absolute value fn
        cell_distances_y = abs.([cell_A_yMin - cell_B_yMin,
                                cell_A_yMin - cell_B_yMax,
                                cell_A_yMax - cell_B_yMin,
                                cell_A_yMax - cell_B_yMax])

        min_grid_dist_y_idx = findmin(cell_distances_y)[2] #Second entry, index of the minimum over the given dimensions.

        #Will take value 1, 2, 3 or 4. Throw error if another value
        #Lookup table based on index value obtained
        #Assign latitude coords
        if min_grid_dist_y_idx == 1 #cell_A min, cell_B min
            lat1 = cell_A_yMin
            lat2 = cell_B_yMin
        elseif min_grid_dist_y_idx == 2 #cell_A min, cell_B max
            lat1 = cell_A_yMin
            lat2 = cell_B_yMax
        elseif min_grid_dist_y_idx == 3 #cell_A max, cell_B min
            lat1 = cell_A_yMax
            lat2 = cell_B_yMin
        elseif min_grid_dist_y_idx == 4 #cell_A max, cell_B max
            lat1 = cell_A_yMax
            lat2 = cell_B_yMax
        else #Unexpected value for min_grid_dist_x_idx, throw error
            error("min_grid_dist_y_idx has value $min_grid_dist_y_idx. min_grid_dist_y_idx must have value 1, 2, 3 or 4.")
        end

        #Return great circle distance
        return (great_circle_distance(lat1, lon1, lat2, lon2))::Float64
    end
end


#-------------------------------------------------------------------------------
### RETURN DISTANCE INDEX, USED TO LOOK UP KERNEL VALUE
#-------------------------------------------------------------------------------
"""
    return_dist_idx_for_kernel(d::Float64)

Return distance index, used to look up kernel value.

Inputs:
- `d::Float64`: Distance between epidemiological units of interest.

Outputs:
- `distIdx::Int64`: Array index to use in kernel lookup array.

Location: Grid.jl
"""
function return_dist_idx_for_kernel(d::Float64)

    if d <= 0.5
        distIdx = 1
    else  #For distances of 0m, use first entry in vector.
        distIdx = round(Int, d)
    end

    return distIdx::Int64
end
