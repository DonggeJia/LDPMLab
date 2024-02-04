# Evaluation of the Connectivity Matrix
struct UniqueElementsStructure
    Projected :: Array{Float64,2}                                                 # Projected Elements composing the connection
    Non_Projected :: Array{Float64,2}                                             # Non-projected Elements composing the connection
    ID1 :: Int64                                                                  # ID of first node of the strut
    ID2 :: Int64                                                                  # ID of second node of the strut
    Tet :: Vector{Int64}                                                          # ID of tetrahedron
    Area :: Vector{Float64}                                                       # Area of each triangle (collated in a vector)
end
struct ConnectivityVector
    ID1 :: Int64                                                                  # ID of first node of the strut
    ID2 :: Int64                                                                  # ID of second node of the strut
    A :: Array{Float64,1}                                                         # Position of first node
    B :: Array{Float64,1}                                                         # Position of second node
    Area :: Float64                                                               # Total (projected) area of the connection
    Centroid :: Array{Float64,1}                                                  # Effective centroid of the connection
    Local_Axes :: Array{Float64,1}                                                # n, l, m local axes (collated in a vector)
    Length :: Float64                                                             # total length
    Length_Agg :: Float64                                                         # Length of the total aggregate portion
    Length_Matrix :: Float64                                                      # Length of the matrix portion
    DOFs :: Array{Int64,1}                                                        # Degrees of Freedom Associated with Nodes I and J
end
struct CompatibilityMatrix
    BN1 :: Array{Float64,1}
    BN2 :: Array{Float64,1}
    BL1 :: Array{Float64,1}
    BL2 :: Array{Float64,1}
    BM1 :: Array{Float64,1}
    BM2 :: Array{Float64,1}
end

function Build_Connect(Elements, Positions, Size)

    Non_Unique_Connections = [Elements.ID1 Elements.ID2];                         # Retrieve all the connections from the Elements structures
    Unique_Connections = unique(Non_Unique_Connections, dims=1);                  # Unique Connections (no repetitions)
    nel = size(Unique_Connections,1);

    Unique_Connections_Elements = Array{UniqueElementsStructure}(undef,nel);      # Initialize the Unique Elements structure
    Connect = Array{ConnectivityVector}(undef,nel);                               # Initialize the Connectivity structure

    
    @inbounds for i=1:nel
        
        Elements_ID = vec(all(Unique_Connections[i,:]' .== Non_Unique_Connections, dims=2));                                               # Indices for all the appearances of current Unique Connection
        Unique_Connections_Elements[i] = UniqueElementsStructure(Elements.Projected[Elements_ID,:], Elements.Non_Projected[Elements_ID,:], Unique_Connections[i,1], Unique_Connections[i,2], Elements.Tet[Elements_ID], Array{Float64,1}(undef,count(Elements_ID)));
        Connection_Area = 0.;
        Connection_Centroid = [0. 0. 0.];
        @inbounds for j=1:count(Elements_ID)
            P1 = Unique_Connections_Elements[i].Projected[j,1:3];                  # Point 1 of Triangle
            P2 = Unique_Connections_Elements[i].Projected[j,4:6];                  # Point 2 of Triangle
            P3 = Unique_Connections_Elements[i].Projected[j,7:9];                  # Point 3 of Triangle
            Triangle_Centroid = (P1+P2+P3)'/3;                                     # Centroid of the Triangle
            Triangle_Area = 1/2*norm(cross(P2-P1,P3-P1));
            Connection_Area += Triangle_Area;                                      # Connected Area is Updated per Every Triangle
            Connection_Centroid = (Connection_Centroid*(Connection_Area-Triangle_Area)+Triangle_Centroid*Triangle_Area)/Connection_Area;    # Refresh Center of Mass Value
            Unique_Connections_Elements[i].Area[j] = Triangle_Area;
        end
        id_1 = Unique_Connections[i,1];
        id_2 = Unique_Connections[i,2];
        A = Positions[id_1,:];                                                      # Node A of the Strut
        B = Positions[id_2,:];                                                      # Node B of the Strut
        len = norm(B-A);                                                            # Length of the Strut
        Agg_1 = Size[id_1];                                                         # Diameter of the 1st Aggregate
        Agg_2 = Size[id_2];                                                         # Diameter of the 2nd Aggregate
        len_agg = (Agg_1 + Agg_2)/2;                                                # Total Length of Aggregate Included in the Strut
        len_m = len - len_agg;                                                      # Effective Matrix Length
        uv_1 = (B-A)/norm(B-A);                                                     # Unit Vector normal to the strut
        if uv_1[3] == 0 && uv_1[1] == uv_1[2];                                      # Choose between the two previously defined vectors
            perp = [uv_1[2]-uv_1[3], uv_1[1], uv_1[1]];
        else
            perp = [uv_1[3], uv_1[3], uv_1[1]-uv_1[2]];
        end
        uv_2 = cross(perp,uv_1);
        uv_2 /= norm(uv_2);
        uv_3 = cross(uv_1,uv_2);
        dofs = collect([gdl_n*id_1-(gdl_n-1):gdl_n*id_1;gdl_n*id_2-(gdl_n-1):gdl_n*id_2]);
        Connect[i] = ConnectivityVector(id_1, id_2, A, B, Connection_Area, vec(Connection_Centroid), [uv_1; uv_2; uv_3], len, len_agg, len_m, dofs);
    end
    return Connect, Unique_Connections, Unique_Connections_Elements, nel
end

Connect, Unique_Connections, Unique_Connections_Elements, nel= Build_Connect(Elements, Positions, Size);




#
# Compatibility Matrix
#
function Compatibility_Matrix(Connect)
B = Array{CompatibilityMatrix}(undef, nel);
@inbounds for i=1:nel
    I = Connect[i].A;
    J = Connect[i].B;
    Center = Connect[i].Centroid;
    A1 = [1. 0. 0. 0. Center[3]-I[3] I[2]-Center[2]; 0. 1. 0. I[3]-Center[3] 0. Center[1]-I[1]; 0. 0. 1. Center[2]-I[2] I[1]-Center[1] 0.];    # B_1, is 3*6
    A2 = [1. 0. 0. 0. Center[3]-J[3] J[2]-Center[2]; 0. 1. 0. J[3]-Center[3] 0. Center[1]-J[1]; 0. 0. 1. Center[2]-J[2] J[1]-Center[1] 0.];    # B_2, is 3*6
    BN1 = vec(1/Connect[i].Length*Connect[i].Local_Axes[1:3]'*A1);                    #  1*6= (1*3)*(3*6)
    BN2 = vec(1/Connect[i].Length*Connect[i].Local_Axes[1:3]'*A2);                    #  1*6= (1*3)*(3*6)
    BL1 = vec(1/Connect[i].Length*Connect[i].Local_Axes[4:6]'*A1);                    #  1*6= (1*3)*(3*6)
    BL2 = vec(1/Connect[i].Length*Connect[i].Local_Axes[4:6]'*A2);                    #  1*6= (1*3)*(3*6)
    BM1 = vec(1/Connect[i].Length*Connect[i].Local_Axes[7:9]'*A1);                    #  1*6= (1*3)*(3*6)
    BM2 = vec(1/Connect[i].Length*Connect[i].Local_Axes[7:9]'*A2);                    #  1*6= (1*3)*(3*6)
    B[i] = CompatibilityMatrix(BN1, BN2, BL1, BL2, BM1, BM2);                         # B is 3*12, B_[i,1] = {[-n*A1 n*A2; -l*A1 l*A2; -m*A1 m*A2]}; # Full Compatibility Matrix
    
end
return B
end
B = Compatibility_Matrix(Connect);



## Steel_steel element and steel_concrete_bond element
steel_uniques_ini = [[[i, i+1] for i = 1:length(nodes_S_total[1])-1]; [[i, i+1] for i = length(nodes_S_total[1])+1:length(nodes_S_total[1])+length(nodes_S_total[2])-1]; [[i, i+1] for i = length(nodes_S_total[1])+length(nodes_S_total[2])+1:length(nodes_S_total[1])+length(nodes_S_total[2])+ length(nodes_S_total[3])-1]]

steel_uniques = [steel_uniques_ini[i] + [n_nodes, n_nodes] for i in eachindex(steel_uniques_ini)]

steel_bond_uniques = [[i+n_nodes, nodes_surround_steel_sequence[i]] for i in eachindex(nodes_surround_steel_sequence)]

len_steel_bond = []
for j = 1:length(traverse_distribution)

append!(len_steel_bond, (norm(nodes_S_total[j][1]-nodes_S_total[j][2])+norm(nodes_S_total[j][1]-[nodes_S_total[j][1][]))/2)
for i = 2:length(nodes_S_total[j])-1
append!(len_steel_bond, (norm(nodes_S_total[j][i]-nodes_S_total[j][i+1])+norm(nodes_S_total[j][i]-nodes_S_total[j][i-1]))/2)
end
