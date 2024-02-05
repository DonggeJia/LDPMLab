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
F_residual = internal - external
#Point_loading=Geometry_parameters[13,1];
output_direction_dis = zeros(Int64, gdl);
#output_displacement_directions
for region in output_displacement_directions
    output_points = intersect(findall(region[1][1, 1] .<= Positions[:, 1] .<= region[1][1, 2]), findall(region[1][2, 1] .<= Positions[:, 2] .<= region[1][2, 2]), findall(region[1][3, 1] .<= Positions[:, 3] .<= region[1][3, 2]))
    for output_point in output_points # Nodes at the Bottom left are Restrained in UX, UY UZ, and UR1
        for index in region[2]
            output_direction_dis[gdl_n*(output_point-1)+index] = 1 ## bottom left nodes are restrained at x y, z, and ROTX direction
        end
    end
end
output_direction_load = zeros(Int64, gdl);
#output_load_directions
for region in output_load_directions
    output_points = intersect(findall(region[1][1, 1] .<= Positions[:, 1] .<= region[1][1, 2]), findall(region[1][2, 1] .<= Positions[:, 2] .<= region[1][2, 2]), findall(region[1][3, 1] .<= Positions[:, 3] .<= region[1][3, 2]))
    for output_point in output_points # Nodes at the Bottom left are Restrained in UX, UY UZ, and UR1
        for index in region[2]
            output_direction_load[gdl_n*(output_point-1)+index] = 2 ## bottom left nodes are restrained at x y, z, and ROTX direction
        end
    end
end

displace_12 = displace[output_direction_dis.==1, 1:step_interval:end]
F_residual1 = F_residual[output_direction_load.==2, 1:step_interval:end]

Force_disp_collect_load = [sum(F_residual1[:, i]) for i in axes(F_residual1, 2)]
Force_disp_collect_dis = [mean(displace_12[:, i]) for i in axes(displace_12, 2)]
Force_disp_collect = hcat(Force_disp_collect_dis, Force_disp_collect_load)
CSV.write(string(load_dis_out_name, "load_dis.csv"), Tables.table(Force_disp_collect), writeheader=false)




loading_dis_mark = output_direction_load + output_direction_dis

if plot_dis_load_region == "Yes"
    text_ = []
    for i in axes(Positions, 1)
        if loading_dis_mark[i*gdl_n-5:i*gdl_n] == [0, 0, 0, 0, 0, 0]
            push!(text_, " ")
        else
            push!(text_, join(map(string, loading_dis_mark[i*gdl_n-5:i*gdl_n]))) # 1 indicates displacement extraction; 2 indicates load, 3 indicates both 
        end
    end
    plt3d = PlotlyJS.scatter3d(; x=Positions[:, 1], y=Positions[:, 2], z=Positions[:, 3], text=text_,
        mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))
    layout = Layout(margin=attr(l=40, r=40, t=40, b=40), scene_camera=attr(
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
