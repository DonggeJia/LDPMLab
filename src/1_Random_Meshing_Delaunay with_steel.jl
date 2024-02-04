cd(@__DIR__)
using LinearAlgebra
using DelimitedFiles
using Statistics
using CSV
using Tables
using TetGen
using PlotlyJS
mutable struct RawTetGenIO{T}
    pointlist::Matrix
end
# read values of the Geometrical parameters

height_S = 20 #mm
traverse_distribution = [dimen2 / 4, dimen2 / 4 * 2, dimen2 / 4 * 3]
diameter_S = 4 #mm

Positions, Size = pointsresult, pointsdiameterfinal
#
## Delaunay
#
input = TetGen.RawTetGenIO{Cdouble}(pointlist=Positions')
triangulation = TetGen.tetrahedralize(input, "Q")
TRI = triangulation.tetrahedronlist'
#Positions[TRI[2,:],2:3]
JJKK = [[1, 2], [2, 3], [3, 4], [1, 3], [1, 4], [2, 4]]
eliminate_TRI = []
for i in axes(TRI, 1) 
    for LLKKJ in JJKK
        poin3 = Positions[TRI[i, LLKKJ[1]], 2:3]
        poin2 = Positions[TRI[i, LLKKJ[2]], 2:3]
        for lolp in traverse_distribution
            poin1 = [lolp, height_S]
            A = (poin2[2] - poin3[2])
            B = (poin3[1] - poin2[1])
            C = poin2[1] * poin3[2] - poin3[1] * poin2[2]
            t = dot(poin1-poin2,poin3-poin2)/dot(poin3-poin2,poin3-poin2)
            if 0<t<1 && (abs(A * poin1[1] + B * poin1[2] + C) / sqrt(A^2 + B^2)) < diameter_S / 2
                append!(eliminate_TRI, i)
                @goto next_tetrahe
            end
        end
    end
    @label next_tetrahe
end
#Positions[setdiff(1:end, unique(reshape(TRI[eliminate_TRI,:],:))),:]
TRI = TRI[setdiff(1:end, eliminate_TRI), :]
Positions = Positions[sort(unique(reshape(TRI,:))),:]
Size = Size[sort(unique(reshape(TRI,:))),:]
## 


n_nodes = size(Positions, 1)                                                                # 6 degrees for each nodes 
gdl_n = 6;                                                                               # 6 degrees for each nodes 
gdl = gdl_n * n_nodes;                                                                     # Total Number of DOFs


#
## check volume is 0
#

n_tet = size(TRI, 1);
V = ones(n_tet)
for i in 1:n_tet
    A = Positions[TRI[i, 1], :]                                                                                                      # 1st Node of the Connection
    B = Positions[TRI[i, 2], :]                                                                                                      # 2nd Node of the Connection
    C = Positions[TRI[i, 3], :]                                                                                                      # 2nd Node of the Connection
    D = Positions[TRI[i, 4], :]
    b_1 = A - B # Coordinates of P1 in the Reference System for the Current Point
    c_1 = A - C # Coordinates of P2 in the Reference System for the Current Point
    d_1 = A - D # Coordinates of P3 in the Reference System for the Current Point
    V[i] = abs(dot(b_1, cross(c_1, d_1))) / 6 # mass of the Tetrahedron
end

V = sort(V)
V_total = sum(V)


#get positions of steel nodes, steel_element, and steel_concrete_bond_element
#the position of steel
# height_S = 20 #mm
# traverse_distribution = [dimen2 / 4, dimen2 / 4 * 2, dimen2 / 4 * 3]
# diameter_S = 4 #mm
#search particle nodes around steel 1
nodes_surround_S_total = []
nodes_S_total = []
for j in eachindex(traverse_distribution)
    nodes_surround_steel = []
    nodes_steel = []
    for i = 1 : Int(length(Positions)/3)
        if ((Positions[i, 2] - traverse_distribution[j])^2 + (Positions[i, 3] - height_S)^2)^0.5 < da / 2 * 1.5
            push!(nodes_surround_steel, i)
            push!(nodes_steel, [Positions[i, 1], traverse_distribution[j], height_S])
        end
    end
    push!(nodes_surround_S_total, nodes_surround_steel)
    push!(nodes_S_total, nodes_steel)
end

for j = 1:length(traverse_distribution)
    unique!(x -> Positions[x][1], nodes_surround_S_total[j])
    unique!(nodes_S_total[j])
end

nodes_surround_steel_sequence= [[l for l in nodes_surround_S_total[1]]; [l for l in nodes_surround_S_total[2]]; [l for l in nodes_surround_S_total[3]]]


nodes_surround_steel_coor= [[Positions[l,:] for l in nodes_surround_S_total[1]]; [Positions[l,:] for l in nodes_surround_S_total[2]]; [Positions[l,:] for l in nodes_surround_S_total[3]]]
nodes_of_steel_coor = [nodes_S_total[1]; nodes_S_total[2]; nodes_S_total[3]]
nodes_surround_steel_coor_ma = mapreduce(permutedims, vcat, nodes_surround_steel_coor)
nodes_of_steel_coor_ma = mapreduce(permutedims, vcat, nodes_of_steel_coor)

all_nodes_with_steel = cat(Positions, nodes_of_steel_coor_ma; dims =1)

# nodes_on_steel_surface = []
# for i in eachindex(nodes_S_total[1])
#     push!(nodes_on_steel_surface, nodes_of_steel_coor[i]+(nodes_surround_steel_coor[i]-nodes_of_steel_coor[i])/norm(nodes_surround_steel_coor[i]-nodes_of_steel_coor[i])*diameter_S*2)
# end
# nodes_on_steel_surface = mapreduce(permutedims, vcat, nodes_on_steel_surface)

# # Specify the file path
# file_path = "D:\\voro++\\voro++-0.4.6\\examples\\walls\\pack_cylinder"  # Replace with your desired file path

# # Data to write
# data = Any[collect(1:length(nodes_S_total[1])) nodes_on_steel_surface 4*ones(length(nodes_S_total[1]))]

# # Open the file in write mode
# writedlm(file_path, data, ' ')

# plt3d = PlotlyJS.scatter3d(; x=nodes_on_steel_surface[:, 1], y=nodes_on_steel_surface[:, 2], z=nodes_on_steel_surface[:, 3],text=collect(1:length(nodes_on_steel_surface)),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))

# layout = Layout(margin=attr(l=40, r=40, t=40, b=40),
#                 scene=attr(aspectmode="data",
#                            camera=attr(up=attr(x=0, y=0, z=1),
#                                        center=attr(x=0, y=0, z=0),
#                                        eye=attr(x=1.7, y=1.2, z=1.2))))

# ll = plot(plt3d, layout)

# PlotlyJS.add_trace!(ll, scatter3d(; x=nodes_of_steel_coor_ma[:, 1], y=nodes_of_steel_coor_ma[:, 2], z=nodes_of_steel_coor_ma[:, 3],text=collect(1:length(nodes_of_steel_coor)),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(0, 127, 127)")))
# ll


