# The Previous Triangulation defines 12 Elements per each tetrahedron
# whose areas are calculated and projected on the plane perprendicular to
# the connection in order to obtain the effective area of the connection
struct ElementsStructure
    Projected :: Array{Float64,2}
    Non_Projected :: Array{Float64,2}
    ID1 :: Vector{Int64}
    ID2 :: Vector{Int64}
    Tet :: Vector{Int64}
end
function Project_Areas(Positions, TRI,Gravity_final)
    n_tet = size(TRI,1);
    Elements = ElementsStructure(Array{Float64,2}(undef,12*n_tet,9), Array{Float64,2}(undef,12*n_tet,9), Array{Int64,1}(undef,0), Array{Int64,1}(undef,0), Array{Int64,1}(undef,0));

    @inbounds for i=1:n_tet
        A = Positions[TRI[i,1],:];                                                                                                      # 1st Node of the Connection
        B = Positions[TRI[i,2],:];                                                                                                      # 2nd Node of the Connection
        C = Positions[TRI[i,3],:];                                                                                                      # 2nd Node of the Connection
        D = Positions[TRI[i,4],:];                                                                                                      # 2nd Node of the Connection
        G_ABCD =Gravity_final[i,1:3];                                                                                                        # Retrieve G_ABCD fromGravity_final
        M_AB =Gravity_final[i,4:6];                                                                                                          # Retrieve M_AB fromGravity_final
        M_AC =Gravity_final[i,7:9];                                                                                                          # Retrieve M_AC fromGravity_final
        M_AD =Gravity_final[i,10:12];                                                                                                        # Retrieve M_AD fromGravity_final
        M_BC =Gravity_final[i,13:15];                                                                                                        # Retrieve M_BC fromGravity_final
        M_BD =Gravity_final[i,16:18];                                                                                                        # Retrieve M_BD fromGravity_final
        M_CD =Gravity_final[i,19:21];                                                                                                        # Retrieve M_CD fromGravity_final
        G_ABC =Gravity_final[i,22:24];                                                                                                       # Retrieve G_ABC fromGravity_final
        G_ABD =Gravity_final[i,25:27];                                                                                                       # Retrieve G_ABD fromGravity_final
        G_ACD =Gravity_final[i,28:30];                                                                                                       # Retrieve G_ACD fromGravity_final
        G_BCD =Gravity_final[i,31:33];                                                                                                       # Retrieve G_BCD fromGravity_final, Effective Center of Mass of the BCD Triangle
        #
        # AB Connection
        #
        uv_AB = (B-A)/norm(B-A);                                                                                                        # Unit Vector for the AB Direction
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_AB,uv_AB))*uv_AB;                                                                       # Projection of G_ABCD Onto a Plane (go through mid point of AB) Perpendicular to AB
        G_ABC_projected = G_ABC-(dot(G_ABC-M_AB,uv_AB))*uv_AB;                                                                          # Projection of G_ABC Onto a Plane Perpendicular to AB
        G_ABD_projected = G_ABD-(dot(G_ABD-M_AB,uv_AB))*uv_AB;                                                                          # Projection of G_ABD Onto a Plane Perpendicular to AB
        Elements.Projected[12i-11:12i-10,:] = [M_AB' G_ABCD_projected' G_ABC_projected'; M_AB' G_ABCD_projected' G_ABD_projected'];    # Points of the Triangle Projected onto the Plane
        Elements.Non_Projected[12i-11:12i-10,:] = [M_AB' G_ABCD' G_ABC'; M_AB' G_ABCD' G_ABD'];                                        # Points of the Triangle
        id1 = min(TRI[i,1], TRI[i,2]);                                                                                                  # Retrieve Lowest ID Value
        id2 = max(TRI[i,1], TRI[i,2]);                                                                                                  # Retrieve Highest ID Value
        push!(Elements.ID1, id1, id1);                                                                                                 # ID1
        push!(Elements.ID2, id2, id2);                                                                                                 # ID2
        push!(Elements.Tet, i, i);                                                                                                     # Tetrahedron
        #
        # AC Connection
        #
        uv_AC = (C-A)/norm(C-A); 
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_AC,uv_AC))*uv_AC; 
        G_ABC_projected = G_ABC-(dot(G_ABC-M_AC,uv_AC))*uv_AC; 
        G_ACD_projected = G_ACD-(dot(G_ACD-M_AC,uv_AC))*uv_AC; 
        Elements.Projected[12i-9:12i-8,:] = [M_AC' G_ABCD_projected' G_ABC_projected'; M_AC' G_ABCD_projected' G_ACD_projected']; 
        Elements.Non_Projected[12i-9:12i-8,:] = [M_AC' G_ABCD' G_ABC'; M_AC' G_ABCD' G_ACD']; 
        id1 = min(TRI[i,1], TRI[i,3]); 
        id2 = max(TRI[i,1], TRI[i,3]); 
        push!(Elements.ID1, id1, id1); 
        push!(Elements.ID2, id2, id2); 
        push!(Elements.Tet, i, i); 
        #
        # AD Connection
        #
        uv_AD = (D-A)/norm(D-A);                                                                                                    
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_AD,uv_AD))*uv_AD;                                                                    
        G_ABD_projected = G_ABD-(dot(G_ABD-M_AD,uv_AD))*uv_AD;                                                                       
        G_ACD_projected = G_ACD-(dot(G_ACD-M_AD,uv_AD))*uv_AD;                                                                       
        Elements.Projected[12i-7:12i-6,:] = [M_AD' G_ABCD_projected' G_ABD_projected'; M_AD' G_ABCD_projected' G_ACD_projected'];    
        Elements.Non_Projected[12i-7:12i-6,:] = [M_AD' G_ABCD' G_ABD'; M_AD' G_ABCD' G_ACD'];                                       
        id1 = min(TRI[i,1], TRI[i,4]);                                                                                             
        id2 = max(TRI[i,1], TRI[i,4]);                                                        
        push!(Elements.ID1, id1, id1);                              
        push!(Elements.ID2, id2, id2);                                                                                          
        push!(Elements.Tet, i, i);                                                                       
        #
        # BC Connection
        #
        uv_BC = (C-B)/norm(C-B); 
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_BC,uv_BC))*uv_BC; 
        G_ABC_projected = G_ABC-(dot(G_ABC-M_BC,uv_BC))*uv_BC; 
        G_BCD_projected = G_BCD-(dot(G_BCD-M_BC,uv_BC))*uv_BC;
        Elements.Projected[12i-5:12i-4,:] = [M_BC' G_ABCD_projected' G_ABC_projected'; M_BC' G_ABCD_projected' G_BCD_projected']; 
        Elements.Non_Projected[12i-5:12i-4,:] = [M_BC' G_ABCD' G_ABC'; M_BC' G_ABCD' G_BCD']; 
        id1 = min(TRI[i,2], TRI[i,3]); 
        id2 = max(TRI[i,2], TRI[i,3]); 
        push!(Elements.ID1, id1, id1); 
        push!(Elements.ID2, id2, id2); 
        push!(Elements.Tet, i, i); 
        #
        # BD Connection
        #
        uv_BD = (D-B)/norm(D-B)
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_BD,uv_BD))*uv_BD; 
        G_ABD_projected = G_ABD-(dot(G_ABD-M_BD,uv_BD))*uv_BD; 
        G_BCD_projected = G_BCD-(dot(G_BCD-M_BD,uv_BD))*uv_BD; 
        Elements.Projected[12i-3:12i-2,:] = [M_BD' G_ABCD_projected' G_ABD_projected'; M_BD' G_ABCD_projected' G_BCD_projected']; 
        Elements.Non_Projected[12i-3:12i-2,:] = [M_BD' G_ABCD' G_ABD'; M_BD' G_ABCD' G_BCD']; 
        id1 = min(TRI[i,2], TRI[i,4]); 
        id2 = max(TRI[i,2], TRI[i,4]); 
        push!(Elements.ID1, id1, id1);
        push!(Elements.ID2, id2, id2); 
        push!(Elements.Tet, i, i);
        #
        # CD Connection
        #
        uv_CD = (D-C)/norm(D-C)
        G_ABCD_projected = G_ABCD-(dot(G_ABCD-M_CD,uv_CD))*uv_CD;                                                                  # Projection of G_ABCD Onto a Plane Perpendicular to CD
        G_ACD_projected = G_ACD-(dot(G_ACD-M_CD,uv_CD))*uv_CD;                                                                     # Projection of G_ACD Onto a Plane Perpendicular to CD
        G_BCD_projected = G_BCD-(dot(G_BCD-M_CD,uv_CD))*uv_CD;                                                                     # Projection of G_BCD Onto a Plane Perpendicular to CD
        Elements.Projected[12i-1:12i,:] = [M_CD' G_ABCD_projected' G_ACD_projected'; M_CD' G_ABCD_projected' G_BCD_projected'];   # Points of the Triangle Projected onto the Plane
        Elements.Non_Projected[12i-1:12i,:] = [M_CD' G_ABCD' G_ACD'; M_CD' G_ABCD' G_BCD'];                                       # Points of the Triangle
        id1 = min(TRI[i,3], TRI[i,4]);                                                                                             # Retrieve Lowest ID Value
        id2 = max(TRI[i,3], TRI[i,4]);                                                                                             # Retrieve Highest ID Value
        push!(Elements.ID1, id1, id1);                                                                                            # ID1
        push!(Elements.ID2, id2, id2);                                                                                            # ID2
        push!(Elements.Tet, i, i);                                                                                                # Tetrahedron
    end
    return Elements
end
Elements = Project_Areas(Positions, TRI,Gravity_final);

#Triangles=Project_Areas(Positions, TRI,Gravity_final);
# ###check
# n_elements_check=size(Elements.ID1,1)
# strut_check=zeros(n_elements_check,2)
# strut_check[:,1]=Elements.ID1[:]
# strut_check[:,2]=Elements.ID2[:]
# strut_unique= unique(strut_check, dims=1); 
