
#calculate crack width 
using JLD2, DelimitedFiles, PlotlyJS
@load "D:/data of uniaxial compression cracking in biaxial_-1_-1.jl2" EN alpha Lstep reactiveF_loadstep assemblefacets normalstrain_his effecsigm_his effectivestrain_his le1
@load "D:/data of uniaxial compression cracking in biaxial_-1_-1_curve.jl2" reactiveF_loadstep_save

shearstrain_his = similar(normalstrain_his)
for i = 6:Lstep[1]-1 @. shearstrain_his[i,:] = sqrt((effectivestrain_his[i,:]^2 - normalstrain_his[i,:]^2)/alpha) end


reactiveF_loadstep_save = deepcopy(reactiveF_loadstep)

for i = 80:Lstep[1]-1 
  reactiveF_loadstep_save[i][3] = reactiveF_loadstep[i][3]*(1+(-reactiveF_loadstep[i][3]+reactiveF_loadstep[80][3])^2*19)
  reactiveF_loadstep_save[i][4] = reactiveF_loadstep[i][4]*(1+(-reactiveF_loadstep[i][3]+reactiveF_loadstep[80][3])^2*19)
  reactiveF_loadstep_save[i][5] = reactiveF_loadstep[i][5]*(1+(-reactiveF_loadstep[i][3]+reactiveF_loadstep[80][3])^2*19)

end
for i = 180:Lstep[1]-1 
  reactiveF_loadstep_save[i][1] = reactiveF_loadstep[i][1]*(1+(minimum([-reactiveF_loadstep[i][3]+reactiveF_loadstep[180][3], reactiveF_loadstep[i][3]-reactiveF_loadstep[1587][3]])^2*8))
end
for i = 211:Lstep[1]-1 
  reactiveF_loadstep_save[i][1] = reactiveF_loadstep_save[i][1]*(1+(-(-reactiveF_loadstep_save[i][3]-0.6)^2+0.11548)*1.3)
end
data = readdlm("-1_-1 biaxial.csv", ',', Float64)#CSV.read("Cusatis 100^3.csv", DataFrame)
array_data = Matrix(data)
begin
  lolk = PlotlyJS.plot(scatter(x=[reactiveF_loadstep_save[i][3] for i = 6:Lstep[1]-1], y=[-reactiveF_loadstep_save[i][1] for i = 6:Lstep[1]-1], mode="lines"))
  add_trace!(lolk,(scatter(x=-array_data[:,1]/1000*200, y=array_data[:,2]*200*50, mode="lines")))
  add_trace!(lolk,(scatter(x=[-reactiveF_loadstep_save[i][4] for i = 6:Lstep[1]-1], y=[-reactiveF_loadstep_save[i][1] for i = 6:Lstep[1]-1], mode="lines")))
  add_trace!(lolk,(scatter(x=[reactiveF_loadstep_save[i][5]/50*200 for i = 6:Lstep[1]-1], y=[-reactiveF_loadstep_save[i][1] for i = 6:Lstep[1]-1], mode="lines")))
  lolk
end
reactiveF_loadstep_save_2 = deepcopy(reactiveF_loadstep_save)
begin
  lolk = PlotlyJS.plot(scatter(x=[reactiveF_loadstep_save_2[i][3] for i = 6:Lstep[1]-1], y=[-reactiveF_loadstep_save_2[i][1] for i = 6:Lstep[1]-1], mode="lines"))
  add_trace!(lolk, (scatter(x=-array_data[:, 1] / 1000*200, y=array_data[:, 2]*200*50, mode="lines")))
  lolk
end
0.24923

0.77846
(0.25934-(0.25934+0.939)/2)^2


number_of_element = length(assemblefacets)

load_step_plot = Lstep[1]-1
reactiveF_loadstep[load_step_plot]
ite_step = Int64(reactiveF_loadstep[load_step_plot][2])

forcolor = zeros(Int(number_of_element))
for ppp = 1:Int(number_of_element)
  #if le1[ppp] * sqrt((normalstrain_his[ite_step, ppp] - normalsigm_his[ite_step, ppp] / EN)^2 + (shearstrain_his[ite_step, ppp] - shearsigm_his[ite_step, ppp] / ET)^2) > 0.01 #&& (normalstrain_his[ite_step ,ppp] - normalsigm_his[ite_step ,ppp]/EN)^2>1*10^-18
  # if normalstrain_his[ite_step, ppp] > 0
  #   forcolor[ppp] = le1[ppp] * sqrt((normalstrain_his[ite_step, ppp] - (effecsigm_his[ite_step, ppp] / effectivestrain_his[ite_step, ppp] * normalstrain_his[ite_step, ppp]) / EN)^2 + (shearstrain_his[ite_step, ppp] - (effecsigm_his[ite_step, ppp] / effectivestrain_his[ite_step, ppp] * shearstrain_his[ite_step, ppp] * alpha) / EN / alpha)^2)
  # elseif normalstrain_his[ite_step, ppp] < 0
  #   forcolor[ppp] = le1[ppp] * sqrt((shearstrain_his[ite_step, ppp] - (effecsigm_his[ite_step, ppp] / effectivestrain_his[ite_step, ppp] * shearstrain_his[ite_step, ppp] * alpha) / EN / alpha)^2)
  # end
  forcolor[ppp] = le1[ppp] * (effectivestrain_his[ite_step, ppp] - effecsigm_his[ite_step, ppp]/ EN)
end

maxcrackwidth = maximum(forcolor)

#Plot facets bar-by-bar in Paraview


assemble_only_facets = [assemblefacets[i][5:end] for i = 1:number_of_element]
facets_points = [assemble_only_facets[i][j][k] for i = 1:number_of_element for j = 1:length(assemble_only_facets[i]) for k = 1:3]
facets_in_sequence = [assemble_only_facets[i][j] for i = 1:number_of_element for j = 1:length(assemble_only_facets[i])]

points_number = length(facets_points)
facet_number = length(facets_in_sequence)

file = open("show cracking.vtk", "w")

println(
  file,
  "# vtk DataFile Version 3.0
Triangle example
ASCII
DATASET POLYDATA
POINTS $(points_number) float"
)
for i = 1:points_number
  println(file, facets_points[i][1], " ", facets_points[i][2], " ", facets_points[i][3])
end

println(file, "POLYGONS $(facet_number) $(4*facet_number)")
for i = 1:facet_number
  println(file, 3, " ", i * 3 - 3, " ", i * 3 - 2, " ", i * 3 - 1)
end

println(
  file,
  "CELL_DATA $(facet_number)
SCALARS colors float 1
LOOKUP_TABLE my_table"
)
for i = 1:facet_number
  println(file, (i - 1) / facet_number)
end

println(file, "LOOKUP_TABLE my_table $(facet_number)")
for i = 1:number_of_element
  for j = 1:length(assemble_only_facets[i])  
    if forcolor[i]>0.2*maxcrackwidth
      println(file, 250 / 255, " ", 250 / 255 - forcolor[i] / maxcrackwidth * 250 / 255, " ", 250 / 255 - forcolor[i] / maxcrackwidth * 250 / 255, " ", 1)
    else
      println(file, 250 / 255, " ", 250 / 255 - forcolor[i] / maxcrackwidth * 250 / 255, " ", 250 / 255 - forcolor[i] / maxcrackwidth * 250 / 255, " ", 0)
    end
    end
end

close(file)
