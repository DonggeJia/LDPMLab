# using DelimitedFiles, BenchmarkTools
point_original_inf = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\particle_orginal_information.custom1")
#(1) ID of particle; coordinates of particle; radius of particle; (4)coordinates of vertices surronding a particle; 
#(5) A list of bracketed sequence of vertices that make up each face; (6)areas of each face; (7) volume; (8) neighboring particle
point_ID = point_original_inf[:, 1]
partical_amount = length(point_ID)
point_coordinates = point_original_inf[:, 2:4]
particle_radius = point_original_inf[:, 5]
surrounding_vertices = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\surrounding vertices coordinates.custom1")
surrounding_vertices = [filter!(x -> x != "", surrounding_vertices[i, :]) for i = 1:partical_amount]
surrounding_vertices = [[parse.(Float64, match(r"\((.*),(.*),(.*)\)", surrounding_vertices[i][j]).captures) for j = 1:length(surrounding_vertices[i])] for i = 1:partical_amount]
sequence_vertices_offace = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\sequence of vertices of faces.custom1")
sequence_vertices_offace = [filter!(x -> x != "", sequence_vertices_offace[i, :]) for i = 1:partical_amount]
sequence_vertices_offace = [[parse.(Int64, [m.match for m in eachmatch(r"\d+", str)]) .+ 1 for str in sequence_vertices_offace[i]] for i = 1:partical_amount]
areas_offace = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\areas of each face.custom1")
areas_offace = [filter!(x -> x != "", areas_offace[i, :]) for i = 1:partical_amount]
volume7 = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\volume.custom1")
neighboring_particles = readdlm("D:\\voro-master-ready\\voro-master\\examples\\custom\\neighboring particle.custom1")
neighboring_particles = [filter!(x -> x != "", neighboring_particles[i, :]) for i = 1:partical_amount]
#plot bars
# for a particle i, its jth face
using PlotlyJS
trace2 = []
# ser = [[1 3 2]; [2 4 4]]
for jjj = 1:partical_amount
    for kkk = 1:length(sequence_vertices_offace[jjj])
        severalpoints = surrounding_vertices[jjj][sequence_vertices_offace[jjj][kkk]]
        for lll = 1:length(severalpoints)-1
            push!(trace2, [[severalpoints[lll][1], severalpoints[lll+1][1]], [severalpoints[lll][2], severalpoints[lll+1][2]], [severalpoints[lll][3], severalpoints[lll+1][3]]])
        end
    end
end
unique!(trace2)
unique!(x -> x[1][1] * x[1][2] + x[2][1] * x[2][2] + x[3][1] * x[3][2], trace2)

trace = GenericTrace{Dict{Symbol,Any}}[]
# ser = [[1 3 2]; [2 4 4]]
for i = 1:length(trace2)
    push!(trace, PlotlyJS.scatter3d(; x=trace2[i][1], y=trace2[i][2], z=trace2[i][3], mode="lines"))

end
plot(trace)
#calculate the intersection between bars and cylinder
#function of the cylinder
using SymPy

steel_radius = steel_dia / 2
for i = 1:length(trace2)
    if (trace2[i][1][1] - steel_center_line_x)^2 + (trace2[i][3][1] - steel_center_line_z)^2 < steel_radius^2 && (trace2[i][1][2] - steel_center_line_x)^2 + (trace2[i][3][2] - steel_center_line_z)^2 < steel_radius^2
        trace2[i] = 0.0
    elseif (trace2[i][1][1] - steel_center_line_x)^2 + (trace2[i][3][1] - steel_center_line_z)^2 < steel_radius^2
        println(i, "_1")


        x, y, z = symbols("x y z")
        linsolve((x + y + z - 1, x + y + 2*z - 3, x - y), (x, y, z))
        
        result = (nlsolve(find_intersection, vec(zeros(3, 1)))).zero

        if abs(result[1][1] - trace2[i][1][2]) < abs(trace2[i][1][1] - trace2[i][1][2])
            trace2[i][1][1] = result[1][1]
            trace2[i][2][1] = result[1][2]
            trace2[i][3][1] = result[1][3]
        else
            trace2[i][1][1] = result[2][1]
            trace2[i][2][1] = result[2][2]
            trace2[i][3][1] = result[2][3]
        end

    elseif (trace2[i][1][2] - steel_center_line_x)^2 + (trace2[i][3][2] - steel_center_line_z)^2 < steel_radius^2
        println(i, "_2")

        function find_intersection(F, x)
            F[1] = (x[1] - steel_center_line_x)^2 + (x[3] - steel_center_line_z)^2 - steel_radius^2
            F[2] = (x[1] - trace2[i][1][2]) / (trace2[i][1][1] - trace2[i][1][2]) - (x[2] - trace2[i][2][2]) / (trace2[i][2][1] - trace2[i][2][2])
            F[3] = (x[1] - trace2[i][1][2]) / (trace2[i][1][1] - trace2[i][1][2]) - (x[3] - trace2[i][3][2]) / (trace2[i][3][1] - trace2[i][3][2])
            F[3] = (x[3] - trace2[i][3][2]) / (trace2[i][3][1] - trace2[i][3][2]) - (x[2] - trace2[i][2][2]) / (trace2[i][2][1] - trace2[i][2][2])
        end
        result = (nlsolve(find_intersection, vec(zeros(3, 1))))#.zero
        if abs(result[1][1] - trace2[i][1][1]) < abs(trace2[i][1][2] - trace2[i][1][1])
            trace2[i][1][2] = result[1][1]
            trace2[i][2][2] = result[1][2]
            trace2[i][3][2] = result[1][3]
        else
            trace2[i][1][2] = result[2][1]
            trace2[i][2][2] = result[2][2]
            trace2[i][3][2] = result[2][3]
        end
    end
end
using NLsolve

function f!(F, x)
    F[1] = x[1]^2-18
end
res = (nlsolve(f!, [4.0])).zero
using SymPy
x, y, z = symbols("x y z")
result= linsolve((x + y + z - 1, x + y + 2*z - 3, x - y), (x, y, z))
lol

lol = Symbol(result)
[i for i in Symbol]
lil = String(lol)
a, b, c, d = symbols("a, b, c, d", real=True)
nonlinsolve([b + a, a - 2*b], [a, b])
int_s = parse.(, match(r".*\((.*),(.*),(.*)\).*", lil).captures)