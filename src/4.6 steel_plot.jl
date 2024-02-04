
file = open("show steels.vtk", "w")

println(
    file,
    "# vtk DataFile Version 3.0
3D spheres example
ASCII
DATASET POLYDATA
POINTS $length(traverse_distribution)*2 float"
)
#for i=1:length(pointsfinal)
#println(file, facets_points[1][1], " ", facets_points[1][2]+3, " ", facets_points[1][3]+3)
for j in eachindex(traverse_distribution)
    println(file, 0.0, " ", traverse_distribution[j], " ", height_S)
    println(file, dimen1, " ", traverse_distribution[j], " ", height_S)
end
#end

println(file, "LINES $length(traverse_distribution) $length(traverse_distribution)*3")
for i in eachindex(traverse_distribution)
    println(file, 2, " ", 2 * i - 2, " ", 2 * i - 1)
end
#println(file, "\n") # This will add another line feed

close(file)

#run(`paraview /vtk_files/encastre_deck_$(steps[1]).vtk`)



# @save string("Result_crack", Modelname,".jld2")  w Unique_Connections


# maximum(w[:,1830])


# w_select=zeros(size(w,1),4)
# w_select[:,2]=w[:,1830]
# w_select[:,3]=w[:,3700]
# w_select[:,4]=w[:,5964]

# @save string("Result_crack", Modelname,"select_4step.jld2")  w_select Unique_Connections
# CSV.write("w4000_theory.csv",  Tables.table(w[:,4000]), writeheader=false)
# CSV.write("len.csv",  Tables.table(len), writeheader=false)
# CSV.write("sigma_N.csv",  Tables.table(sigma_N[:,3000]), writeheader=false)
# CSV.write("eps_n3000.csv",  Tables.table(eps_n[:,3000]), writeheader=false)
# CSV.write("eps_t3000.csv",  Tables.table(eps_t[:,3000]), writeheader=false)
# CSV.write("sigma_t3000.csv",  Tables.table(sigma_T[:,3000]), writeheader=false)
# Positions[Unique_Connections[127,:],:]
# Positions[Unique_Connections[150,:],:]
