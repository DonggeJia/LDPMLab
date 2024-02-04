#Set the Boundary Conditions on Restrained Nodes
# 0 = Free
# 1 = Fixed
# 2 = Fixed (Used to Compute the FD Curve, must be in the direction of the applied load)
# 3 = Imposed Displacement

function BoundaryConditions(Positions)

    width = Geometry_parameters[1,1];     
    #depth = Geometry_parameters[2,1];                                                             ### depth of the specimen (Y-dir)                                                       
    mean_d_x = Geometry_parameters[8,1];                                                          ### Mean Distance Between Particles to Be Placed on the External Surfaces (X-dir)
    #mean_d_y = Geometry_parameters[9,1];                                                          ### Mean Distance Between Particles to Be Placed on the External Surfaces (Y-dir)
    #mean_d_z = Geometry_parameters[10,1]; 
    support_start= Geometry_parameters[12]
    support_end  = width-Geometry_parameters[12]       
    Boun = zeros(Int64, gdl);

    Point_loading=Geometry_parameters[13,1];

######## bottom nodes, simply-support beam
    bottom_all= findall(Positions[:,3] .== 0 ); # Bottom Nodes of the Cube
    bottom_left = []
    bottom_right = []
    for i in bottom_all
        if support_start-mean_d_x .<= Positions[i,1].<= support_start+mean_d_x
            append!(bottom_left,i)
        elseif support_end-mean_d_x .<= Positions[i,1].<= support_end+mean_d_x
            append!(bottom_right,i)
        end
    end
    # for i in bottom_all
    #     if  Positions[i,1].== support_start+mean_d_x
    #         append!(bottom_left,i)
    #     elseif  Positions[i,1].== support_end-mean_d_x
    #         append!(bottom_right,i)
    #     end
    # end   



    @inbounds for i=1:size(bottom_left,1) # Nodes at the Bottom left are Restrained in UX, UY UZ, and UR1
        id_boun = bottom_left[i,1];
        Boun[gdl_n*(id_boun-1)+1:gdl_n*(id_boun-1)+3] .= 1; ## bottom left nodes are restrained at x y, z, and ROTX direction
        Boun[gdl_n*id_boun-3] = 2; # Restrained in Z-dir (Used to check for Reaction force)
    end


    @inbounds for i=1:size(bottom_right,1) # Nodes at the Bottom Right are Restrained in Y and Z-dir
        id_boun = bottom_right[i,1];
        Boun[gdl_n*(id_boun-1)+2:gdl_n*(id_boun-1)+3] .= 1; # # bottom rignt nodes are restrained at y, z, and ROTX direction
        Boun[gdl_n*id_boun-3] = 2; # Restrained in Z-dir (Used to check for Reaction force)
    end


######## top nodes, only restrain z_direction (apply disp on z_direction)
    top_all= findall(Positions[:,3] .== Geometry_parameters[3,1] ); 
    top_middle = []
    for i in top_all
        ## if  Positions[i,1].>= Geometry_parameters[1,1]/2-12.5 && Positions[i,1].<= Geometry_parameters[1,1]/2+12.5


        if  (Point_loading-mean_d_x/2 .<= Positions[i,1].<= Point_loading+mean_d_x/2 || width-Point_loading-mean_d_x/2 .<= Positions[i,1].<= width-Point_loading+mean_d_x/2  )
            append!(top_middle,i)
        end


    end
    
    @inbounds for i=1:size(top_middle,1); # apply disp to Top_middle Nodes 
        id_boun = top_middle[i,1];

        Boun[gdl_n*id_boun-3] = 3; # apply disp in Z-dir (Used to check for Reaction force)
    end
    
    return Boun, bottom_left, bottom_right,top_middle 
    end

    Boun, bottom_left, bottom_right,top_middle  = BoundaryConditions(Positions);

    # println("Number of restrained DOFs = ", size(findall(Boun.==1),1))
    # println("Number of imposed displacement DOFs = ", size(findall(Boun.>1),1))


