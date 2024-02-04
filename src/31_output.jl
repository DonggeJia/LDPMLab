


function PlotForceDisp(bottom_left,bottom_right,top_middle,delta )
    ## total force
    bot_node=[bottom_left;bottom_right]
    #bot_node= top_middle
    S3_bottom=zeros(size(bot_node,1),size(external,2))   # node force at z
    for j in 1:size(bot_node,1)
        S3_bottom[j,:]=external[6*bot_node[j]-3,:]
    end
    S3_bottom_sum=zeros(size(external,2));
    for i in 1:size(external,2)
        S3_bottom_sum[i]=sum(S3_bottom[:,i])  
    end
    ## disp
    disp_z_top=zeros(size(top_middle,1),size(external,2))
    for j in 1:size(top_middle,1)
        disp_z_top[j,:]=displace[6*top_middle[j]-3,:]
    end

    ## detailed curve
    Force_disp_detail=zeros(size(external,2),2)
    Force_disp_detail[:,1]=(-disp_z_top[1,:])'
    # Force_disp_detail[:,2]=S3_bottom_sum[:]/Geometry_parameters[2,1]
    Force_disp_detail[:,2]=S3_bottom_sum[:]/1000
    # ## simply curve
    inc = size(displace,2) # Number of Converged Steps to Plot
    steps = collect(1:delta:inc)

    Force_disp_collect=zeros(size(steps,1),2)
    for i in eachindex(steps)
        step = steps[i]
        Force_disp_collect[i,1] = Force_disp_detail[step,1]
        Force_disp_collect[i,2] = Force_disp_detail[step,2]
    end



   

    return Force_disp_detail,Force_disp_collect
end



Force_disp_detail,Force_disp_collect = PlotForceDisp(bottom_left,bottom_right,top_middle,100 )

# ## CSV.write("force_disp_detail.csv",  Tables.table(Force_disp_detail), writeheader=false)
CSV.write(string(Modelname,".csv"),  Tables.table(Force_disp_collect), writeheader=false)





# for i in [bottom_left;bottom_right;top_middle] 
#     println(Positions[i,:])
# end








# #### abstract from top nodes
# # disp
# restrain_top_nodes = findall(Boun.==3)   ## # Restrained top nodes in Z-dir 
# disp_z_top=displace[restrain_top_nodes,:]   # grab disp for all top nodes

# #
# ## top nodes
# #
# Reaction_top=external[restrain_top_nodes,:] # grab force for all top nodes
# ReactionForce_top=zeros(size(Reaction_top,2)) # define total force of top nodes 
# for i in 1:size(Reaction_top,2)
#     ReactionForce_top[i]=sum(Reaction_top[:,i])
# end
# TopNode_Force_disp=[-disp_z_top[1,:] ReactionForce_top/Geometry_parameters[2,1]]
# CSV.write("TopNode_Force_disp.csv",  Tables.table(TopNode_Force_disp), writeheader=false)



# #bot_node= top_middle
# S1_bottom=zeros(size(bot_node,1),size(external,2))
# S2_bottom=zeros(size(bot_node,1),size(external,2))
# S3_bottom=zeros(size(bot_node,1),size(external,2))
# S4_bottom=zeros(size(bot_node,1),size(external,2))
# S5_bottom=zeros(size(bot_node,1),size(external,2))
# S6_bottom=zeros(size(bot_node,1),size(external,2))

# for j in 1:size(bot_node,1)
#     S1_bottom[j,:]=external[6*bot_node[j]-5,:]
#     S2_bottom[j,:]=external[6*bot_node[j]-4,:]
#     S3_bottom[j,:]=external[6*bot_node[j]-3,:]
#     S4_bottom[j,:]=external[6*bot_node[j]-2,:]
#     S5_bottom[j,:]=external[6*bot_node[j]-1,:]
#     S6_bottom[j,:]=external[6*bot_node[j]-0,:]
# end


# S1_bottom_sum=zeros(size(external,2));
# S2_bottom_sum=zeros(size(external,2));
# S3_bottom_sum=zeros(size(external,2));
# S4_bottom_sum=zeros(size(external,2));
# S5_bottom_sum=zeros(size(external,2));
# S6_bottom_sum=zeros(size(external,2));

# for i in 1:size(external,2)
#     S1_bottom_sum[i]=sum(S1_bottom[:,i])
#     S2_bottom_sum[i]=sum(S2_bottom[:,i])
#     S3_bottom_sum[i]=sum(S3_bottom[:,i])
# end
