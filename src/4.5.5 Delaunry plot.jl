function mesh_plot_Delauney(filename)
file0 = string(filename, "Delauney.vtk")
file = open(file0, "w")

println(file, "# vtk DataFile Version 3.0
3D lines example
ASCII
DATASET POLYDATA
POINTS $(Int(gdl/6)) float")

for i=1:Int(gdl/6)
println(file, Positions[i,1], " ", Positions[i,2], " ", Positions[i,3])
end
println(file, "LINES $(size(Unique_Connections, 1) - size(steel_uniques, 1) - size(steel_bond_uniques, 1)) $(3*(size(Unique_Connections, 1) - size(steel_uniques, 1) - size(steel_bond_uniques, 1)))")
for i=1:size(Unique_Connections, 1) - size(steel_uniques, 1) - size(steel_bond_uniques, 1)
  println(file, 2, " ", Connect[i].ID1-1, " ", Connect[i].ID2-1)
end
#println(file, "\n") # This will add another line feed

close(file)
end

mesh_plot_Delauney(filename)