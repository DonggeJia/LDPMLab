#Set the Boundary Conditions on Restrained Nodes
# 0 = Free
# 1 = Fixed
# 3 = Imposed Displacement

    #width = Geometry_parameters[1,1];     
    #depth = Geometry_parameters[2,1];                                                             ### depth of the specimen (Y-dir)                                                       
    #mean_d_x = Geometry_parameters[8,1];                                                          ### Mean Distance Between Particles to Be Placed on the External Surfaces (X-dir)
    #mean_d_y = Geometry_parameters[9,1];                                                          ### Mean Distance Between Particles to Be Placed on the External Surfaces (Y-dir)
    #mean_d_z = Geometry_parameters[10,1]; 
    #support_start= Geometry_parameters[12]
    #support_end  = width-Geometry_parameters[12]       
    Boun = zeros(gdl);
    Boun1 = zeros(Int64, gdl);
    #Point_loading=Geometry_parameters[13,1];
    #fixed_region = [[[0 10;0 200;65 70],[3,4,5]],[[100 110;0 200;65 70],[3,4,5]],[[190 200;0 200;65 70],[3,4,5]]]
    for region in fixed_region
        fixed_points = intersect(findall(region[1][1,1] .<= Positions[:,1] .<= region[1][1,2]), findall(region[1][2,1] .<= Positions[:,2] .<= region[1][2,2]), findall(region[1][3,1] .<= Positions[:,3] .<= region[1][3,2]));
        for fixed_point in fixed_points # Nodes at the Bottom left are Restrained in UX, UY UZ, and UR1
            for index in region[2]
                Boun[gdl_n*(fixed_point-1)+index] = 0.0; ## bottom left nodes are restrained at x y, z, and ROTX direction
                Boun1[gdl_n*(fixed_point-1)+index] = 1;
                #Boun[gdl_n*bottom_left[i]-3] = 2; # Restrained in Z-dir (Used to check for Reaction force)
            end
        end
    end

    #loaded_region = [[[0 10;0 200;65 70],[3,4,5]],[[100 110;0 200;65 70],[3,4,5]],[[190 200;0 200;65 70],[3,4,5]]]
    
    for region in loaded_region
        loaded_points = intersect(findall(region[1][1,1] .<= Positions[:,1] .<= region[1][1,2]), findall(region[1][2,1] .<= Positions[:,2] .<= region[1][2,2]), findall(region[1][3,1] .<= Positions[:,3] .<= region[1][3,2]));
        println(loaded_points)
        for loaded_point in loaded_points # Nodes at the Bottom left are Restrained in UX, UY UZ, and UR1
            for i in eachindex(region[2])
                Boun[gdl_n*(loaded_point-1)+region[2][i]] = region[3][i]; ## bottom left nodes are restrained at x y, z, and ROTX direction
                Boun1[gdl_n*(loaded_point-1)+region[2][i]] = 3;
                #Boun[gdl_n*bottom_left[i]-3] = 2; # Restrained in Z-dir (Used to check for Reaction force)
            end
        end
    end

    if plot_boundary == "Yes"
        text_ = []
        for i in axes(Positions,1)
            push!(text_, join(map(string, Boun[i*gdl_n-5:i*gdl_n])))
        end
        plt3d = PlotlyJS.scatter3d(; x=Positions[:, 1], y=Positions[:, 2], z=Positions[:, 3],text=text_,
        mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))
        layout = Layout(margin=attr(l=40, r=40, t=40, b=40), scene_camera = attr(
            up=attr(x=0, y=0, z=1),
            center=attr(x=0, y=0, z=0),
            eye=attr(x=1.7, y=1.2, z=1.2)))
        points_position = plot(plt3d, layout)
        display(points_position)
    end

        ######## bottom nodes, simply-support beam
    #bottom_all= findall(Positions[:,3] .== 0 ); # Bottom Nodes of the Cube
    
    # for i in bottom_all
    #     if support_start-mean_d_x .<== Positions[i,1].<== support_start+mean_d_x
    #         append!(bottom_left,i)
    #     elseif support_end-mean_d_x .<== Positions[i,1].<== support_end+mean_d_x
    #         append!(bottom_right,i)
    #     end
    # end
    # for i in bottom_all
    #     if  Positions[i,1].== support_start+mean_d_x
    #         append!(bottom_left,i)
    #     elseif  Positions[i,1].== support_end-mean_d_x
    #         append!(bottom_right,i)
    #     end
    # end   


    

    # println("Number of restrained DOFs = ", size(findall(Boun.==1),1))
    # println("Number of imposed displacement DOFs = ", size(findall(Boun.>1),1))
    #top_all = findall(Positions[:,3] .== dimen3)
# plt3d = PlotlyJS.scatter3d(; x=Positions[bottom_left, 1], y=Positions[bottom_left, 2], z=Positions[bottom_left, 3],text=collect(1:length(Positions)/3),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))
# layout = Layout(margin=attr(l=40, r=40, t=40, b=40), scene_camera = attr(
#     up=attr(x=0, y=0, z=1),
#     center=attr(x=0, y=0, z=0),
#     eye=attr(x=1.7, y=1.2, z=1.2)))
# points_position = plot(plt3d, layout)
