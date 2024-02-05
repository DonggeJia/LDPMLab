
using Format
function vtk_cracks(sigma_n, sigma_t, epsilon_n, epsilon_t, E, α, Unique_Connections_Triangles, len, displace, steps, filename)
    
    w = zeros(size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques)), size(displace,2))


    for ii = 1:size(displace,2)
        for jj = 1:size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques))
            if sigma_n[jj,ii]>0

                #w_n = len[jj]*(epsilon_n[jj,ii]-sigma_n[jj,ii]/E[jj])
                #w_t = len[jj]*(epsilon_t[jj,ii]-sigma_t[jj,ii]/E[jj]/α)
                 w_n = len[jj]*(epsilon_n[jj,ii])
                 w_t = len[jj]*(epsilon_t[jj,ii])
                ## w_t = 0
                w[jj,ii] = sqrt(w_n^2 + w_t^2)
            end
        end
    end


    num_cel = zeros(size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques)),1);
    for i=1:size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques))
        num_cel[i,1] = Int64(size(Unique_Connections_Triangles[i].Non_Projected,1));
    end

    num_triangles = Int64(sum(num_cel));
    num_points = 3*num_triangles;
    num_entries = num_triangles*4;

    cell_types_decimal = ones(num_triangles,1)*5;
    cell_types = [floor(Int,x) for x in cell_types_decimal]


    connectivity = Array{Int64,2}(undef, num_triangles,4)


    ## connectivity = zeros(num_triangles,4);
    for i=1:num_triangles
        connectivity[i,:] = [3 i*3-3 i*3-2 i*3-1];
    end


    if isdir("vtk_files/") == false # Check if directory to contain post-processing files exists
        mkdir("vtk_files") # If not, create directory
    end

    for mm = 1:size(steps,1)
        ii = steps[mm]

        file = string(filename,"_", ii, ".vtk")
        open(file, "w") do f
            write(f, "# vtk DataFile Version 4.0 \n")
            write(f, "Unstructured grid legacy vtk file with point scalar data \n")
            write(f, "ASCII \n")
            write(f, "DATASET UNSTRUCTURED_GRID \n")
            write(f, "POINTS $num_points double \n")




            #### coordinator  ???
            for jj = 1:size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques))
               for kk=1:Int64(num_cel[jj,1])
                write(f,   "$(Unique_Connections_Triangles[jj].Non_Projected[kk,1]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,2]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,3]) \n")
                write(f,   "$(Unique_Connections_Triangles[jj].Non_Projected[kk,4]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,5]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,6]) \n")
                write(f,   "$(Unique_Connections_Triangles[jj].Non_Projected[kk,7]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,8]) $(Unique_Connections_Triangles[jj].Non_Projected[kk,9]) \n")
                end
            end



            write(f, "CELLS $num_triangles $num_entries \n")

            for jj = 1:size(connectivity,1)
                write(f, "$(connectivity[jj,1]) $(connectivity[jj,2]) $(connectivity[jj,3])  $(connectivity[jj,4])\n")
            end

            write(f, "CELL_TYPES $num_triangles\n")

            for ii =1: size(cell_types,1)
            write(f, " $(cell_types[ii]) \n")                  
            end


            write(f, "CELL_DATA  $num_triangles\n")            
            write(f, "SCALARS crack_openings double\n")                
            write(f, "LOOKUP_TABLE default\n") 
            for jj = 1:size(Unique_Connections_Triangles,1)-(length(steel_uniques) + length(steel_bond_uniques))
               for k=1:Int64(num_cel[jj,1])
                write(f, "$(w[jj,ii]) \n")                
               end
            end


        end

    end
   return cell_types,w

end

Unique_Connections_Triangles=Unique_Connections_Elements

len=[Connect[i].Length for i in axes(Connect,1)]


n_inc = size(displace,2)-1
steps = [Int(round(relative_time_of_cracking[1]*n_inc)),Int(round(relative_time_of_cracking[2]*n_inc)),Int(round(relative_time_of_cracking[3]*n_inc))]
cell_types,w=vtk_cracks(sigma_N, sigma_T, eps_n, eps_t, E, alpha, Unique_Connections_Triangles,len, displace, steps, crack_plot_dirc_and_name)

