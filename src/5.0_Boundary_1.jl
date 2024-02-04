#Set the Boundary Conditions on Restrained Nodes
# 0 = Free
# 1 = Fixed
# 2 = Fixed (Used to Compute the FD Curve, must be in the direction of the applied load)
# 3 = Imposed Displacement

function BoundaryConditions(Positions)
    
    Boun = zeros(Int64, gdl);

######## bottom nodes, simply-support beam
    bottom_all= findall(Positions[:,3] .== 0 ); # Bottom Nodes of the Cube
    bottom_left = []
    bottom_right = []
    for i in bottom_all
        if Positions[i,1].== 25
            append!(bottom_left,i)
        elseif Positions[i,1].==475
            append!(bottom_right,i)
        end
    end
    
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
        if  Positions[i,1].== Geometry_parameters[1,1]/2
            append!(top_middle,i)
        end
    end
    
    @inbounds for i=1:size(top_middle,1); # apply disp to Top_middle Nodes 
        id_boun = top_middle[i,1];

        Boun[gdl_n*id_boun-3] = 3; # apply disp in Z-dir (Used to check for Reaction force)
    end
    
    return Boun
    end
    Boun = BoundaryConditions(Positions);

    # println("Number of restrained DOFs = ", size(findall(Boun.==1),1))
    # println("Number of imposed displacement DOFs = ", size(findall(Boun.>1),1))



# for j in bottom_left
#     print(Positions[j,:])
# end
# for j in bottom_right
#     print(Positions[j,:])
# end
# for j in top_middle
#     print(Positions[j,:])
# end


# function BoundaryConditions(Positions)
#     Boun = zeros(Int64, gdl);
#     top = findall(Positions[:,3] .>= 145); # Top Nodes of the Cube
#     bot = findall(Positions[:,3] .<= 5); # Bottom Nodes of the Cube
#     # bot_2 = find(Positions(:,3)==0 & Positions(:,1)>=390); % Bottom Nodes 2
#     # left = find(Positions(:,2)<=0);
#     # right = find(Positions(:,2)>= 100);
#     # top = find(Positions(:,3)>= 45 & Positions(:,2)>= 200 & Positions(:,2)<=300); % Restrained Top Nodes
#     # bot_1 = find(Positions(:,3)<=10 & Positions(:,2)<= 25); % Bottom Nodes 1
#     # bot_2 = find(Positions(:,3)<=10 & Positions(:,2)>=475); % Bottom Nodes 2
#     @inbounds for i=1:size(bot,1) # Nodes at the Bottom are Restrained in X, Y and Z-dir
#         id_boun = bot[i,1];
#         Boun[gdl_n*(id_boun-1)+1:gdl_n*id_boun] .= 1; # Fully Restrained
#         Boun[gdl_n*id_boun-3] = 2; # Restrained in Z-dir (Used to check for Reaction force)
#     end
    
#     @inbounds for i=1:size(top,1); # Top Nodes are Fully Restrained
#         id_boun = top[i,1];
#         Boun[gdl_n*(id_boun-1)+1:gdl_n*id_boun] .= 1; # Fully Restrained
#         Boun[gdl_n*id_boun-3] = 3; # Restrained in Z-dir (Used to check for Reaction force)
#     end
    
#     return Boun
#     end
#     Boun = BoundaryConditions(Positions);
#     println("Number of restrained DOFs = ", size(findall(Boun.==1),1))
#     println("Number of imposed displacement DOFs = ", size(findall(Boun.>1),1))
    


    

