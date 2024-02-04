using DelimitedFiles
#rm("tetrahedra_mechanical_100.dat")
tetrahedron = readdlm("randompoints_Kupher.dat")

using PlotlyJS
barpointpair = Vector{Int64}[]
for rawite = 1:Int(length(tetrahedron) / 4)
    push!(barpointpair, [tetrahedron[rawite, 1], tetrahedron[rawite, 2]])
    push!(barpointpair, [tetrahedron[rawite, 1], tetrahedron[rawite, 3]])
    push!(barpointpair, [tetrahedron[rawite, 1], tetrahedron[rawite, 4]])
    push!(barpointpair, [tetrahedron[rawite, 2], tetrahedron[rawite, 3]])
    push!(barpointpair, [tetrahedron[rawite, 2], tetrahedron[rawite, 4]])
    push!(barpointpair, [tetrahedron[rawite, 3], tetrahedron[rawite, 4]])
end


barpointpair2 = deepcopy(barpointpair)
for i=1:Int(length(barpointpair))
  if barpointpair2[i][1]<barpointpair2[i][2]
    inter2 = barpointpair2[i][1]
    barpointpair2[i][1]=barpointpair2[i][2]
    barpointpair2[i][2]=inter2  
  end
end
barpointpair2 = unique(barpointpair2)

# using SymRCM
# conn = mapreduce(permutedims, vcat, barpointpair2);
# nfens = maximum(conn);
# adjgr = SymRCM.adjgraph(conn, nfens)
# degrees = SymRCM.nodedegrees(adjgr)
# numbering1 = symrcm(adjgr, degrees)



# for deletsame1 = 1:Int(length(barpointpair))
#     ggg = 0
#     nnn = Int(length(barpointpair2))
#     deletsame2 = 1
#     while deletsame2 <= nnn
#         if barpointpair[deletsame1] == barpointpair2[deletsame2]||(barpointpair[deletsame1][1] == barpointpair2[deletsame2][2]&&barpointpair[deletsame1][2] == barpointpair2[deletsame2][1])
#             ggg = ggg + 1
#             if ggg != 1
#                 splice!(barpointpair2, deletsame2)
#                 nnn = Int(length(barpointpair2))
#                 deletsame2 -= 1
#             end
#         end
#         deletsame2 += 1
#     end
# end
# function f1(x)
#     d = countmap(x)
#     return [key for (key, val) in d]
# end
# barpointpair2 = f1(barpointpair)

# using PlotlyJS
# trace = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(barpointpair2))+1)
# # ser = [[1 3 2]; [2 4 4]]
# for jjj = 1:Int(length(barpointpair2))
#     twopoints = [pointsfinal[barpointpair2[jjj][1]], pointsfinal[barpointpair2[jjj][2]]]
#         trace[jjj] = PlotlyJS.scatter3d(; x=[twopoints[1][1], twopoints[2][1]], y=[twopoints[1][2], twopoints[2][2]], z=[twopoints[1][3], twopoints[2][3]], showlegend=false, mode="lines",
#         line=attr(size=12, width=4,color ="rgb(200,165,132)"))
# end
# twopoints = [pointsfinal[barpointpair2[56][1]], pointsfinal[barpointpair2[56][2]]]
# trace[Int(length(barpointpair2))+1] = PlotlyJS.scatter3d(; x=[twopoints[1][1], twopoints[2][1]], y=[twopoints[1][2], twopoints[2][2]], z=[twopoints[1][3], twopoints[2][3]], showlegend=false, mode="lines",
# line=attr(size=10, width=10,color ="red"))
# LOPL = PlotlyJS.plot(trace)
# open("./#56 element.html", "w") do io
#     PlotlyBase.to_html(io, LOPL.plot)
# end
# trace[2] = PlotlyJS.scatter3d(; x=ser[:, 1], y=ser[:, 2], z=ser[:, 3], mode="lines",
#     marker=attr(color="#1f77b4", size=12, symbol="circle",
#         line=attr(color="rgb(0,0,0)", width=0)),
#     line=attr(color="#1f77b4", width=1))
# trace[3] = PlotlyJS.scatter3d(; x=ser[:, 1], y=ser[:, 2], z=ser[:, 3], mode="lines",
#     marker=attr(color="#1f77b4", size=12, symbol="circle",
#         line=attr(color="rgb(0,0,0)", width=0)),
#     line=attr(color="#1f77b4", width=1))


tttt = convert(Matrix{Int}, tetrahedron)
verticeoftetrahedron = pointsfinal[tttt]
using LinearAlgebra
voronoifaces = Array{Any}(undef, Int(length(tetrahedron) / 4))
for tetranumber = 1:Int(length(tetrahedron) / 4)
       pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 2])
       distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 2]) /pointsdistance
    edgecenter1 =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 2]])/2+pointsradiusfinal[tttt[tetranumber, 2]])*distancedirection+verticeoftetrahedron[tetranumber, 2]
    facecenter1pre = (verticeoftetrahedron[tetranumber, 3] + edgecenter1) / 2
       pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 3])
       distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 3]) /pointsdistance
    edgecenter2 =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 3]])/2+pointsradiusfinal[tttt[tetranumber, 3]])*distancedirection+verticeoftetrahedron[tetranumber, 3]
    facecenter2pre = (verticeoftetrahedron[tetranumber, 2] + edgecenter2) / 2
       pointsdistance = norm(verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 3])
       distancedirection = (verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 3]) /pointsdistance
    edgecenter3 =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 2]] - pointsradiusfinal[tttt[tetranumber, 3]])/2+pointsradiusfinal[tttt[tetranumber, 3]])*distancedirection+verticeoftetrahedron[tetranumber, 3]
    facecenter3pre = (verticeoftetrahedron[tetranumber, 1] + edgecenter3) / 2
    facecenter1_1 = (facecenter3pre + facecenter2pre + facecenter1pre) / 3
    centroid1pre = (facecenter1_1 + verticeoftetrahedron[tetranumber, 4]) / 2

      pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 2])
      distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 2]) /pointsdistance
    edgecenter1_ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 2]])/2+pointsradiusfinal[tttt[tetranumber, 2]])*distancedirection+verticeoftetrahedron[tetranumber, 2]
    facecenter1pre_ = (verticeoftetrahedron[tetranumber, 4] + edgecenter1_) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 4])
      distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 4]) /pointsdistance
    edgecenter2_ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 4]])/2+pointsradiusfinal[tttt[tetranumber, 4]])*distancedirection+verticeoftetrahedron[tetranumber, 4]
    facecenter2pre_ = (verticeoftetrahedron[tetranumber, 2] + edgecenter2_) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 4])
      distancedirection = (verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 4]) /pointsdistance
    edgecenter3_ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 2]] - pointsradiusfinal[tttt[tetranumber, 4]])/2+pointsradiusfinal[tttt[tetranumber, 4]])*distancedirection+verticeoftetrahedron[tetranumber, 4]
    facecenter3pre_ = (verticeoftetrahedron[tetranumber, 1] + edgecenter3_) / 2
    facecenter2_2 = (facecenter3pre_ + facecenter2pre_ + facecenter1pre_) / 3
    centroid2pre = (facecenter2_2 + verticeoftetrahedron[tetranumber, 3]) / 2

      pointsdistance = norm(verticeoftetrahedron[tetranumber, 3] - verticeoftetrahedron[tetranumber, 2])
      distancedirection = (verticeoftetrahedron[tetranumber, 3] - verticeoftetrahedron[tetranumber, 2]) /pointsdistance
    edgecenter1__ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 3]] - pointsradiusfinal[tttt[tetranumber, 2]])/2+pointsradiusfinal[tttt[tetranumber, 2]])*distancedirection+verticeoftetrahedron[tetranumber, 2]
    facecenter1pre__ = (verticeoftetrahedron[tetranumber, 4] + edgecenter1__) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 3] - verticeoftetrahedron[tetranumber, 4])
      distancedirection = (verticeoftetrahedron[tetranumber, 3] - verticeoftetrahedron[tetranumber, 4]) /pointsdistance
    edgecenter2__ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 3]] - pointsradiusfinal[tttt[tetranumber, 4]])/2+pointsradiusfinal[tttt[tetranumber, 4]])*distancedirection+verticeoftetrahedron[tetranumber, 4]
    facecenter2pre__ = (verticeoftetrahedron[tetranumber, 2] + edgecenter2__) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 4])
      distancedirection = (verticeoftetrahedron[tetranumber, 2] - verticeoftetrahedron[tetranumber, 4]) /pointsdistance
    edgecenter3__ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 2]] - pointsradiusfinal[tttt[tetranumber, 4]])/2+pointsradiusfinal[tttt[tetranumber, 4]])*distancedirection+verticeoftetrahedron[tetranumber, 4]
    facecenter3pre__ = (verticeoftetrahedron[tetranumber, 3] + edgecenter3__) / 2
    facecenter3_3 = (facecenter3pre__ + facecenter2pre__ + facecenter1pre__) / 3
    centroid3pre = (facecenter3_3 + verticeoftetrahedron[tetranumber, 1]) / 2

      pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 3])
      distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 3]) /pointsdistance
    edgecenter1___ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 3]])/2+pointsradiusfinal[tttt[tetranumber, 3]])*distancedirection+verticeoftetrahedron[tetranumber, 3]
    facecenter1pre___ = (verticeoftetrahedron[tetranumber, 4] + edgecenter1___) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 4] - verticeoftetrahedron[tetranumber, 3])
      distancedirection = (verticeoftetrahedron[tetranumber, 4] - verticeoftetrahedron[tetranumber, 3]) /pointsdistance
    edgecenter2___ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 4]] - pointsradiusfinal[tttt[tetranumber, 3]])/2+pointsradiusfinal[tttt[tetranumber, 3]])*distancedirection+verticeoftetrahedron[tetranumber, 3]
    facecenter2pre___ = (verticeoftetrahedron[tetranumber, 1] + edgecenter2___) / 2
      pointsdistance = norm(verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 4])
      distancedirection = (verticeoftetrahedron[tetranumber, 1] - verticeoftetrahedron[tetranumber, 4]) /pointsdistance
    edgecenter3___ =  ((pointsdistance - pointsradiusfinal[tttt[tetranumber, 1]] - pointsradiusfinal[tttt[tetranumber, 4]])/2+pointsradiusfinal[tttt[tetranumber, 4]])*distancedirection+verticeoftetrahedron[tetranumber, 4]
    facecenter3pre___ = (verticeoftetrahedron[tetranumber, 3] + edgecenter3___) / 2
    facecenter4_4 = (facecenter3pre___ + facecenter2pre___ + facecenter1pre___) / 3
    centroid4pre = (facecenter4_4 + verticeoftetrahedron[tetranumber, 2]) / 2

    centroid = (centroid1pre + centroid2pre + centroid3pre + centroid4pre) / 4
#the order of facets:12edgemidpoint-124facecenter-centroid;
#12edgemidpoint-123facecenter-centroid; 13edgemidpoint-123facecenter-centroid;
#13edgemidpoint-134facecenter-centroid; 23edgemidpoint-123facecenter-centroid;
#23edgemidpoint-234facecenter-centroid; 14edgemidpoint-134facecenter-centroid;
#14edgemidpoint-124facecenter-centroid; 24edgemidpoint-124facecenter-centroid;
#24edgemidpoint-234facecenter-centroid; 34edgemidpoint-134facecenter-centroid;
#34edgemidpoint-234facecenter-centroid;
    voronoifaces[tetranumber] = [[centroid, edgecenter1, facecenter2_2], 
    [centroid, edgecenter1, facecenter1_1], [centroid, edgecenter2, facecenter1_1],
     [centroid, edgecenter2, facecenter4_4], [centroid, edgecenter3, facecenter1_1], 
     [centroid, edgecenter3, facecenter3_3], [centroid, edgecenter2_, facecenter4_4],
      [centroid, edgecenter2_, facecenter2_2], [centroid, edgecenter3_, facecenter2_2], 
      [centroid, edgecenter3_, facecenter3_3], [centroid, edgecenter2___, facecenter4_4], 
      [centroid, edgecenter2___, facecenter3_3]]
end
# tracevoroface = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(tetrahedron) / 4), 12)

# for ppp = 1:Int(length(tetrahedron) / 4)

#     for mmm = 1:12
#         tracevoroface[ppp, mmm] = mesh3d(
#             x=[voronoifaces[ppp][mmm][1][1], voronoifaces[ppp][mmm][2][1], voronoifaces[ppp][mmm][3][1], voronoifaces[ppp][mmm][1][1]],
#             y=[voronoifaces[ppp][mmm][1][2], voronoifaces[ppp][mmm][2][2], voronoifaces[ppp][mmm][3][2], voronoifaces[ppp][mmm][1][2]],
#             z=[voronoifaces[ppp][mmm][1][3], voronoifaces[ppp][mmm][2][3], voronoifaces[ppp][mmm][3][3], voronoifaces[ppp][mmm][1][3]],
#              color="rgb(153, 217, 234)",
#             # Intensity of each vertex, which will be interpolated and color-coded
#             # i, j and k give the vertices of triangles
#             # here we represent the 4 triangles of the tetrahedron surface
#             i=[0, 0, 0, 1],
#             j=[1, 2, 3, 2],
#             k=[2, 3, 1, 3],
#             # name="y",
#             showscale=false)
#     end
#     # for mmm = 1:12
#     #     tracevoroface[ppp, mmm] = isosurface(
#     #         x=[voronoifaces[ppp][mmm][1][1], voronoifaces[ppp][mmm][2][1], voronoifaces[ppp][mmm][3][1], voronoifaces[ppp][mmm][1][1]],
#     #         y=[voronoifaces[ppp][mmm][1][2], voronoifaces[ppp][mmm][2][2], voronoifaces[ppp][mmm][3][2], voronoifaces[ppp][mmm][1][2]],
#     #         z=[voronoifaces[ppp][mmm][1][3], voronoifaces[ppp][mmm][2][3], voronoifaces[ppp][mmm][3][3], voronoifaces[ppp][mmm][1][3]],
#     #         value=[1.4,1.4,1.4,1.4])
#     # end
# end
# tracevoroface1 = reshape(tracevoroface, Int(length(tetrahedron) / 4) * 12)
# layout = Layout(autosize=false, width=500, height=500,
#     margin=attr(l=0, r=0, b=0, t=65), scene_camera = attr(
#         up=attr(x=0, y=0, z=1),
#         center=attr(x=0, y=0, z=0),
#         eye=attr(x=1.7, y=1.2, z=1.2)))
# tesselation_mechanical = PlotlyJS.plot(tracevoroface1, layout)
# #savefig(tesselation_mechanical, "tesselation_mechanical.svg")
# open("./50 facets&elements.html", "w") do io
#     PlotlyBase.to_html(io, tesselation_mechanical.plot)
# end

# pointsofhytra = Array{Any}(missing, length(tetrahedron)/4, 4, 3)
# isosurface(
#     x=[0,0,0,0],
#     y=[1,0,1,1],
#     z=[1,1,0,1],
#     value=[1.4,1.4,1.4,1.4]
#     # isomin=2,
#     # isomax=8,)

# using DelimitedFiles
# writedlm("vertice1oftetrahedron.xls", pointsresult[tttt[:,1],:])
# writedlm("vertice2oftetrahedron.xls", pointsresult[tttt[:,2],:])
# writedlm("vertice3oftetrahedron.xls", pointsresult[tttt[:,3],:])
# writedlm("vertice4oftetrahedron.xls", pointsresult[tttt[:,4],:])
# trace11 = isosurface(
#             x=[voronoifaces[1][1][1][1], voronoifaces[1][1][2][1], voronoifaces[1][1][3][1], voronoifaces[1][1][1][1]],
#             y=[voronoifaces[1][1][1][2], voronoifaces[1][1][2][2], voronoifaces[1][1][3][2], voronoifaces[1][1][1][2]],
#             z=[voronoifaces[1][1][1][3], voronoifaces[1][1][2][3], voronoifaces[1][1][3][3], voronoifaces[1][1][1][3]],
#             value=[1.4,1.4,1.4,1.4])