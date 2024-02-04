
mutable struct RawTetGenIO{T}
    pointlist::Matrix
end
# read values of the Geometrical parameters

Positions, Size = pointsresult, pointsdiameterfinal
#
## Delaunay
#
input = TetGen.RawTetGenIO{Cdouble}(pointlist=Positions')
triangulation = TetGen.tetrahedralize(input, "Q")
TRI = triangulation.tetrahedronlist'


#
## 
#
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
height_S = 40 #mm
traverse_distribution = [dimen2 / 4, dimen2 / 4 * 2, dimen2 / 4 * 3]
diameter_S = 16 #mm
#search particle nodes around steel 1
nodes_surround_S_total = []
nodes_S_total = []
for j in eachindex(traverse_distribution)
    nodes_surround_steel = []
    nodes_steel = []
    for i = 1 : Int(length(Positions)/3)
        if diameter_S/2+7 < ((Positions[i, 2] - traverse_distribution[j])^2 + (Positions[i, 3] - height_S)^2)^0.5 < diameter_S*1+7 
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

nodes_surround_steel_radius = pointsdiameterfinal[nodes_surround_steel_sequence]./2
nodes_surround_steel_coor= [[Positions[l,:] for l in nodes_surround_S_total[1]]; [Positions[l,:] for l in nodes_surround_S_total[2]]; [Positions[l,:] for l in nodes_surround_S_total[3]]]
nodes_of_steel_coor = [nodes_S_total[1]; nodes_S_total[2]; nodes_S_total[3]]
nodes_surround_steel_coor_ma = mapreduce(permutedims, vcat, nodes_surround_steel_coor)
nodes_of_steel_coor_ma = mapreduce(permutedims, vcat, nodes_of_steel_coor)

all_nodes_with_steel = cat(Positions, nodes_of_steel_coor_ma; dims =1)

nodes_on_steel_surface = []
for i in eachindex(nodes_S_total[1])
    push!(nodes_on_steel_surface, nodes_of_steel_coor[i]+(nodes_surround_steel_coor[i]-nodes_of_steel_coor[i])/norm(nodes_surround_steel_coor[i]-nodes_of_steel_coor[i])*diameter_S/2)
end
nodes_on_steel_surface = mapreduce(permutedims, vcat, nodes_on_steel_surface)

steel_uniques_ini = [[[i, i+1] for i = 1:length(nodes_S_total[1])-1]; [[i, i+1] for i = length(nodes_S_total[1])+1:length(nodes_S_total[1])+length(nodes_S_total[2])-1]; [[i, i+1] for i = length(nodes_S_total[1])+length(nodes_S_total[2])+1:length(nodes_S_total[1])+length(nodes_S_total[2])+ length(nodes_S_total[3])-1]]
steel_uniques_ini = mapreduce(permutedims, vcat, steel_uniques_ini)
coor_of_steel_elements = nodes_of_steel_coor[steel_uniques_ini]
coor_of_steel_bond_elements = [[nodes_of_steel_coor_ma[i,:], nodes_surround_steel_coor_ma[i,:]] for i in eachindex(nodes_surround_steel_coor)]

# Specify the file path
# file_path = "D:\\voro++\\voro++-0.4.6\\examples\\walls\\pack_cylinder.txt"  # Replace with your desired file path

# # Data to write
# data = Any[collect(1:length(nodes_S_total[1])) nodes_on_steel_surface]

# # Open the file in write mode
# writedlm(file_path, data, ' ')

using DelimitedFiles
writedlm("points on steel for sectioning steel surface.xls", nodes_on_steel_surface)

plt3d = PlotlyJS.scatter3d(; x=nodes_on_steel_surface[:, 1], y=nodes_on_steel_surface[:, 2], z=nodes_on_steel_surface[:, 3],text=collect(1:length(nodes_on_steel_surface)),
mode="markers+text", opacity=0.9, marker=attr(color="rgb(127, 127, 127)"))

layout = Layout(margin=attr(l=40, r=40, t=40, b=40),
                scene=attr(aspectmode="data",
                           camera=attr(up=attr(x=0, y=0, z=1),
                                       center=attr(x=0, y=0, z=0),
                                       eye=attr(x=1.7, y=1.2, z=1.2))))

ll = plot(plt3d, layout)
#PlotlyJS.add_trace!(ll, (scatter(x=d1, y=F.(d1), mode="lines"))) 
# PlotlyJS.add_trace!(ll, scatter3d(; x=nodes_of_steel_coor_ma[:, 1], y=nodes_of_steel_coor_ma[:, 2], z=nodes_of_steel_coor_ma[:, 3],text=collect(1:length(nodes_of_steel_coor)),
# mode="markers+text", opacity=0.9, marker=attr(color="rgb(0, 127, 127)")))
# ll


function plot_cylinder(p1, p2, radius, cylinder_color)
    n_points=50
    # Length of the cylinder
    h = norm(p2 - p1)
    
    # Direction vector from p1 to p2
    v = (p2 - p1) / h

    # Create a rotation matrix to align the cylinder along the direction vector
    z_axis = [0, 0, 1]
    axis = cross(z_axis, v)
    angle = acos(dot(z_axis, v))
    R = RodriguesVec(axis, angle) # Using Rodrigues' rotation formula

    # Create mesh grid for the parameterized surface
    u = range(0, 2pi, length=n_points)
    v = range(0, h, length=n_points)
    us = ones(length(v)) * u'
    vs = v * ones(length(u))'

    # Surface parameterization
    X = radius * cos.(us)
    Y = radius * sin.(us)
    Z = vs

    # Apply rotation and translation
    XYZ = [X[:] Y[:] Z[:]]
    XYZ_rotated = (R * XYZ')'
    XYZ_translated = XYZ_rotated .+ p1'

    # Reshape back to matrices
    X_rotated = reshape(XYZ_translated[:, 1], size(X))
    Y_rotated = reshape(XYZ_translated[:, 2], size(Y))
    Z_rotated = reshape(XYZ_translated[:, 3], size(Z))

    # Create and plot the surface
    return surface(x=X_rotated, y=Y_rotated, z=Z_rotated, colorscale=[[0, cylinder_color], [1, cylinder_color]], opacity=0.5)
    
end

# Function to create a Rodrigues' rotation matrix
function RodriguesVec(k, θ)
    k = k / norm(k)
    cosθ = cos(θ)
    sinθ = sin(θ)
    return [cosθ + k[1]^2 * (1 - cosθ)        k[1]*k[2]*(1 - cosθ) - k[3]*sinθ  k[1]*k[3]*(1 - cosθ) + k[2]*sinθ;
            k[2]*k[1]*(1 - cosθ) + k[3]*sinθ  cosθ + k[2]^2*(1 - cosθ)         k[2]*k[3]*(1 - cosθ) - k[1]*sinθ;
            k[3]*k[1]*(1 - cosθ) - k[2]*sinθ  k[3]*k[2]*(1 - cosθ) + k[1]*sinθ  cosθ + k[3]^2*(1 - cosθ)]
end

# Example usage
p1 = [1, 2, 3]
p2 = [4, 5, 6]
radius = 1
fk = "rgb($(rand(0:255)),$(rand(0:255)),$(rand(0:255)))"
plot(plot_cylinder(p1, p2, diameter_S/2, fk))



num_points = 50
# Function to generate points for a rounded plan
function generate_rounded_plan(center, radius, normal, num_points=50)
    # Create a rotation matrix to align with the normal vector
    z_axis = [0, 0, 1]
    rotation_axis = cross(z_axis, normal)
    rotation_angle = acos(dot(normal, z_axis) / (norm(normal) * norm(z_axis)))
    rotation_matrix = I + sin(rotation_angle) * skew_symmetric(rotation_axis) +
                      (1 - cos(rotation_angle)) * skew_symmetric(rotation_axis)^2

    # Generate points in the plane
    θ = range(0, 2π, length=num_points)
    r = range(0, radius, length=num_points)
    x, y, z = [], [], []
    for ri in r
        for ti in θ
            p = center + rotation_matrix * [ri * cos(ti), ri * sin(ti), 0]
            push!(x, p[1])
            push!(y, p[2])
            push!(z, p[3])
        end
    end
    return surface(x=reshape(x, num_points, num_points), 
    y=reshape(y, num_points, num_points), 
    z=reshape(z, num_points, num_points),
    colorscale=[[0, "rgb(117,153,189)"], [1, "blue"]],showscale=false)
end

# Skew-symmetric matrix for a vector
function skew_symmetric(v)
    return [  0   -v[3]  v[2];
             v[3]   0   -v[1];
            -v[2]  v[1]   0  ]
end

# Define the parameters for the rounded plan
center = [0, 0, 0]  # Center of the disk
radius = 1          # Radius
normal = [0, 0, 1]  # Normal vector (pointing upwards)

# Generate points for the rounded plan

plot(generate_rounded_plan(center, radius, normal))



function sphere1(r, C)   # r: radius; C: center [cx,cy,cz]
    global n = 40
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = C[1] .+ r*cos.(u) * sin.(v)'
    y = C[2] .+ r*sin.(u) * sin.(v)'
    z = C[3] .+ r*ones(n) * cos.(v)'
    return x, y, z
end



# Function to generate vertices of a regular polygon in 3D space
function polygon_vertices(c, r, h, num_sides)
    # Normalize the normal vector
    h /= norm(h)

    # Generate a basis for the plane
    v1 = cross(h, [1.0, 0.0, 0.0])
    if norm(v1) < 1e-10
        v1 = cross(h, [0.0, 1.0, 0.0])
    end
    v1 /= norm(v1)
    v2 = cross(h, v1)

    # Calculate vertices
    θ = range(0, stop=2π, length=num_sides+1)[1:end-1]
    return [c + r * cos(t) * v1 + r * sin(t) * v2 for t in θ]
end

# Function to plot a regular polygonal face using PlotlyJS
function plot_polygonal_face(c, r, h, pentagen_color)
    vertices = polygon_vertices(c, r, h, 5)

    # Extracting x, y, z coordinates
    xs = [v[1] for v in vertices]
    ys = [v[2] for v in vertices]
    zs = [v[3] for v in vertices]

    # Add the center point
    push!(xs, c[1])
    push!(ys, c[2])
    push!(zs, c[3])

    # Create and plot the mesh
    return mesh3d(x=xs, y=ys, z=zs, color=pentagen_color, opacity=0.5)
end

# Example usage
c = [190.18020497880843, 149.49115824794146, 47.98380110419605]
r = 8
h = [0.0, -0.8349031835171843, 13.099752391573688]
plot_polygonal_face(c, r, h, "blue")


trace = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(coor_of_steel_elements)/2))
trace3 = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(coor_of_steel_elements)/2))

# ser = [[1 3 2]; [2 4 4]]
for jjj = 1:Int(length(coor_of_steel_elements)/2)
    color_3 = "rgb($(rand(0:255)),$(rand(0:255)),$(rand(0:255)))"
    twopoints = [coor_of_steel_elements[jjj,1], coor_of_steel_elements[jjj,2]]
        trace[jjj] = PlotlyJS.scatter3d(; x=[twopoints[1][1], twopoints[2][1]], y=[twopoints[1][2], twopoints[2][2]], z=[twopoints[1][3], twopoints[2][3]], showlegend=false, mode="lines",
        line=attr(size=120, width=4,color =color_3))
        trace3[jjj] = generate_rounded_plan((twopoints[1]+twopoints[2])/2, diameter_S/2, [1,0,0])
end

trace1 = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(coor_of_steel_bond_elements)))
trace4 = Array{GenericTrace{Dict{Symbol,Any}}}(undef, Int(length(coor_of_steel_bond_elements)))

# ser = [[1 3 2]; [2 4 4]]
for jjj = 1:Int(length(coor_of_steel_bond_elements))
    color_1 = "rgb($(rand(0:255)),$(rand(0:255)),$(rand(0:255)))"
    twopoints = [coor_of_steel_bond_elements[jjj][1], coor_of_steel_bond_elements[jjj][2]]
        trace1[jjj] = PlotlyJS.scatter3d(; x=[twopoints[1][1], twopoints[2][1]], y=[twopoints[1][2], twopoints[2][2]], z=[twopoints[1][3], twopoints[2][3]], showlegend=false, mode="lines",
        line=attr(size=12, width=4,color = color_1))
        println((twopoints[2]-twopoints[1])/norm(twopoints[2]-twopoints[1])*diameter_S/2+twopoints[1], 8, twopoints[2]-twopoints[1])
        trace4[jjj] = plot_polygonal_face((twopoints[2]-twopoints[1])/norm(twopoints[2]-twopoints[1])*diameter_S/2+twopoints[1], 8, twopoints[2]-twopoints[1], color_1)
end
traceball = GenericTrace{Dict{Symbol,Any}}[]
for defg = 1:Int(length(nodes_surround_steel_radius))
    if nodes_surround_steel_radius[defg] != 0
    kok = sphere1(nodes_surround_steel_radius[defg], nodes_surround_steel_coor[defg])#pointsradiusfinal[defg], pointsfinal[defg])
    push!(traceball,PlotlyJS.surface(x=kok[1],y=kok[2],z=kok[3],showscale=false,colorscale=[[0, "rgb(189,189,189)"], [1, "rgb(189,189,189)"]]))
    end
end

#plot_cylinder(p1, p2, diameter_S/2)
layout = Layout(
    margin=attr(l=0, r=0, b=0, t=65), scene=attr(aspectmode="data"), scene_camera = attr(
        up=attr(x=0, y=0, z=1),
        center=attr(x=0, y=0, z=0),
        eye=attr(x=1.7, y=1.2, z=1.2)))

demenstrate_steel_bond = PlotlyJS.plot([trace1;trace;trace3;traceball;trace4], layout)

#PlotlyJS.add_trace!(demenstrate_steel_bond, trace3[20]) 

# #savefig(tesselation_mechanical, "tesselation_mechanical.svg")
open("./demenstrate_steel_bond.html", "w") do io
     PlotlyBase.to_html(io, demenstrate_steel_bond.plot)
end

