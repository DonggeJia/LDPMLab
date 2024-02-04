# Import Material Mat_parameters from the "Mat_parameters.txt" file and evaluate the constitutive law
 #Mpa
function ConstitutiveLawMat_parameters(Connect, Mat_parameters)
    # This function calculates the Softening-Controlling Mat_parameters for Every Lattice Element in Order to Preserve the Correct Energy Dissipation
    n_conn = length(Connect);
    n_conn_concrete = size(Unique_Connections, 1)
    K_t = zeros(n_conn);
    K_s = zeros(n_conn);
    E = zeros(n_conn);
    @inbounds for i = 1:n_conn_concrete
        E[i] = Connect[i].Length/(Connect[i].Length_Agg/Mat_parameters[2] + Connect[i].Length_Matrix/Mat_parameters[1]);      # Effective Elastic Modulus
        lcr_t = 2*E[i]*Mat_parameters[6]/(Mat_parameters[3]^2);                                                               # Characteristic Length for Mode I Fracture
        lcr_s = 2*Mat_parameters[8]*E[i]*Mat_parameters[7]/Mat_parameters[5]^2;                                                   # Characteristic Length for Mode II Fracture
        if Connect[i].Length <= lcr_t
            K_t[i] = 2*E[i]/(lcr_t/Connect[i].Length-1);                                                              # Softening Exponential Coefficient, Parameter that Governs the Post-Peak Behavior in Tension
        else
            K_t[i] = 1e5;                                                                                             # If the actual length is greater than the critical value, use a fictitious (very big!) softening modulus
        end
        if Connect[i].Length <= lcr_s
            K_s[i] = 2*Mat_parameters[8]*E[i]/(lcr_s/Connect[i].Length-1);                                                # Parameter that Governs the Post-Peak Behavior in Shear
        else
            K_s[i] = 1e5;                                                                                             # If the actual length is greater than the critical value, use a fictitious (very big!) softening modulus
        end
        #n_t(i,1) = log(K_t(i,1)/(K_t(i,1)-K_s(i,1)))/log(1-2*omega_zero/pi); % Exponent for the Tension-Shear Behavior
    end
    @inbounds for i = n_conn_concrete+1:n_conn_concrete+size(steel_uniques, 1)
        #E[i] = E_steel;      # steel Elastic Modulus
        
        #n_t(i,1) = log(K_t(i,1)/(K_t(i,1)-K_s(i,1)))/log(1-2*omega_zero/pi); % Exponent for the Tension-Shear Behavior
    end
    @inbounds for i = n_conn_concrete+size(steel_uniques, 1)+1:n_conn_concrete+size(steel_uniques, 1)+size(steel_bond_uniques, 1)
        #E[i] = E_bond
        #n_t(i,1) = log(K_t(i,1)/(K_t(i,1)-K_s(i,1)))/log(1-2*omega_zero/pi); % Exponent for the Tension-Shear Behavior
    end
    return E, K_t, K_s
end

function TetrahedronStiffnessMatrix(TRI, Connect, Unique_Connections, Mat_parameters, Elements, E)
    n_tet = size(TRI,1);
    Tet_Elements = Vector{Array}(undef,n_tet);                                                                        # Elements Composing Each Tetrahedron
    Tet_Stiffness = Vector{Array}(undef,n_tet);                                                                       # Stiffness Matrix for Each Tetrahedron
    Tet_Mass = Vector{Array}(undef,n_tet);                                                                            # Mass Matrix for Each Tetrahedron
    Tet_B = Array{Array}(undef,n_tet, 12);                                                                            # Compatibility Matrix B for each Triangle in Each Tetrahedron
    Connectivity = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];                                                                    # Edge Connectivity in a Tetrahedron (general)
    Cell_ = [1 2 3 4 5 6; 1 2 7 8 9 10; 3 4 7 8 11 12; 5 6 9 10 11 12];                                               # Every Row Contains Indices to the Facets that Create one Cell
    @inbounds for i = 1:n_tet
        Tet_Elements[i] = Elements.Projected[i*12-11:i*12,:];
        Tet_Stiffness[i] = zeros(24,24);                                                                              # Initialize Stiffness Matrix for the Current Tetrahedron
        Tet_Mass[i] = zeros(24);                                                                                      # Initialize Mass Matrix for the Current Tetrahedron
        AB = sort([TRI[i,1],TRI[i,2]]);                                                                               # IDs of A and B points (sorted)
        AC = sort([TRI[i,1],TRI[i,3]]);                                                                               # IDs of A and C points (sorted)
        AD = sort([TRI[i,1],TRI[i,4]]);                                                                               # IDs of A and D points (sorted)
        BC = sort([TRI[i,2],TRI[i,3]]);                                                                               # IDs of B and C points (sorted)
        BD = sort([TRI[i,2],TRI[i,4]]);                                                                               # IDs of B and D points (sorted)
        CD = sort([TRI[i,3],TRI[i,4]]);                                                                               # IDs of C and D points (sorted)
        edges_id = [AB'; AC'; AD'; BC'; BD'; CD'];                                                                    # IDs of the 6 Edges
        @inbounds for j=1:6                                                                                           # Cycle Over number of edges
            id = findfirst(all(edges_id[j,:]' .== Unique_Connections, dims=2));
            n = Connect[id].Local_Axes[1:3];
            l = Connect[id].Local_Axes[4:6];
            m = Connect[id].Local_Axes[7:9];
            if TRI[i,Connectivity[j,1]] < TRI[i,Connectivity[j,2]]
                I_ = Connect[id].A;
                J_ = Connect[id].B;
            else
                I_ = Connect[id].B;
                J_ = Connect[id].A;
            end
            Compliance = [E[id] 0 0; 0 Mat_parameters[8]*E[id] 0; 0 0 Mat_parameters[8]*E[id]];
            @inbounds for k=1:2 # Cycle Over the Two Elements that Interest Each Lattice Element
                Center = (Tet_Elements[i][2j+k-2,1:3] + Tet_Elements[i][2j+k-2,4:6] + Tet_Elements[i][2j+k-2,7:9])/3; # Center of Gravity of the Triangle
                A1 = [1. 0. 0. 0. Center[3]-I_[3] I_[2]-Center[2]; 0. 1. 0. I_[3]-Center[3] 0. Center[1]-I_[1]; 0. 0. 1. Center[2]-I_[2] I_[1]-Center[1] 0.]; # B_1
                A2 = [1. 0. 0. 0. Center[3]-J_[3] J_[2]-Center[2]; 0. 1. 0. J_[3]-Center[3] 0. Center[1]-J_[1]; 0. 0. 1. Center[2]-J_[2] J_[1]-Center[1] 0.]; # B_2
                Tet_B[i,2j+k-2] = [-n'*A1 n'*A2; -l'*A1 l'*A2; -m'*A1 m'*A2];
            end
            K_edge = Tet_B[i,2j-1]'*Compliance*Tet_B[i,2j-1] + Tet_B[i,2j]'*Compliance*Tet_B[i,2j]; # Stiffness Matrix of the Current Lattice Element [12x12]
            id1 = Connectivity[j,1];
            id2 = Connectivity[j,2];
            Tet_Stiffness[i][6id1-5:6id1, 6id1-5:6id1] += K_edge[1:6,1:6];
            Tet_Stiffness[i][6id1-5:6id1, 6id2-5:6id2] += K_edge[1:6,7:12];
            Tet_Stiffness[i][6id2-5:6id2, 6id1-5:6id1] += K_edge[7:12,1:6];
            Tet_Stiffness[i][6id2-5:6id2, 6id2-5:6id2] += K_edge[7:12,7:12];
        end
        @inbounds for j=1:4 # Cycle Over the 4 Vertices of the Tetrahedron
            M = zeros(6); # Number of Facets per Node per Tetrahedron
            @inbounds for k=1:6 # Cycle Over Number of Facets Pertaining Each Node
                P1 = Tet_Elements[i][Cell_[j, k],1:3];
                P2 = Tet_Elements[i][Cell_[j, k],4:6];
                P3 = Tet_Elements[i][Cell_[j, k],7:9];
                P4 = Positions[TRI[i,j],:];
                b_1 = P1-P4; # Coordinates of P1 in the Reference System for the Current Point
                c_1 = P2-P4; # Coordinates of P2 in the Reference System for the Current Point
                d_1 = P3-P4; # Coordinates of P3 in the Reference System for the Current Point
                mass = Mat_parameters[11]*abs(dot(b_1,cross(c_1,d_1)))/6; # mass of the Tetrahedron
                x_1 = b_1[1]; # X Coordinate of 1st Point
                x_2 = c_1[1]; # Y Coordinate of 1st Point
                x_3 = d_1[1]; # Z Coordinate of 1st Point
                y_1 = b_1[2]; # X Coordinate of 2nd Point
                y_2 = c_1[2]; # Y Coordinate of 2nd Point
                y_3 = d_1[2]; # Z Coordinate of 2nd Point
                z_1 = b_1[3]; # X Coordinate of 3rd Point
                z_2 = c_1[3]; # Y Coordinate of 3rd Point
                z_3 = d_1[3]; # Z Coordinate of 3rd Point
                I_x = mass/10*(y_1^2+y_1*y_2+y_2^2+y_1*y_3+y_2*y_3+y_3^2+z_1^2+z_1*z_2+z_2^2+z_1*z_3+z_2*z_3+z_3^2);
                I_y = mass/10*(x_1^2+x_1*x_2+x_2^2+x_1*x_3+x_2*x_3+x_3^2+z_1^2+z_1*z_2+z_2^2+z_1*z_3+z_2*z_3+z_3^2);
                I_z = mass/10*(x_1^2+x_1*x_2+x_2^2+x_1*x_3+x_2*x_3+x_3^2+y_1^2+y_1*y_2+y_2^2+y_1*y_3+y_2*y_3+y_3^2);
                M += [mass, mass, mass, I_x, I_y, I_z];
            end
            Tet_Mass[i][6j-5:6j] = M;
        end
    end
    return Tet_Mass, Tet_Stiffness
end


function MassScaling(Tet_Mass, Tet_Stiffness, TRI)
    n_tet = size(TRI,1);
    ω_n = zeros(n_tet); # Initialize Natural Frequency Vector
        @inbounds for i = 1:n_tet
            ei = eigen(Tet_Stiffness[i], diagm(Tet_Mass[i]));
            ω_n[i] = sqrt(maximum(real(ei.values)));
        end
        println("")
        println("Maximum frequency is ", maximum(ω_n), ".")
        println("Average frequency is ", mean(ω_n), ".")
        println("Median frequency is ", median(ω_n), ". Corresponding Δt is ", 2/median(ω_n), ".")
        println("Input desired Δt to be used for mass scaling:")
        # s = readline()               ### input a number here
        # Δt = parse(Float64, chomp(s));
        Δt= round(2/median(ω_n),digits=5);
        max_ω = 2/Δt;
        scale = findall(ω_n .> max_ω)
        @inbounds for i=1:size(scale,1)
            Tet_Mass[scale[i]] *= ω_n[scale[i]]^2/max_ω^2;
        end
    return Δt, Tet_Mass, scale
end



#
function MassMatrix(gdl, TRI, Tet_Mass, damping)
     M = zeros(gdl);
     @inbounds for i=1:size(TRI,1)
         M_ = Tet_Mass[i];
         M[6*TRI[i,1]-5:6*TRI[i,1]] += Tet_Mass[i][1:6];
         M[6*TRI[i,2]-5:6*TRI[i,2]] += Tet_Mass[i][7:12];
         M[6*TRI[i,3]-5:6*TRI[i,3]] += Tet_Mass[i][13:18];
         M[6*TRI[i,4]-5:6*TRI[i,4]] += Tet_Mass[i][19:24];
     end
     C = M.*damping; # Calculate Damping Matrix
     return M, C
end


Mat_parameters = LDPM_w_steel.mechanical_parameters; # [E_m E_a σ_t σ_c σ_s G_t G_s α n_t n_c ρ kc_1 kc_2 K_c ζ]
const k1 = -Mat_parameters[3]*Mat_parameters[4]/Mat_parameters[5]^2; # Stress Space Elliptic Boundary Parameter
const k2 = -Mat_parameters[3]*Mat_parameters[4]; # Stress Space Elliptic Boundary Parameter
const sigma_t = Mat_parameters[3];# Extract Tension Limit from Mat_parameters Vector
const sigma_c = Mat_parameters[4]; # Extract Compression Limit from Mat_parameters Vector
const sigma_s = Mat_parameters[5]; # Extreact Shear Limit from Mat_parameters Vector
const alpha = Mat_parameters[8]; # Extract Coupling Parameter
const K_c = Mat_parameters[14]; # Extract Compressive Exponential Initial Slop Parameter
const kc2 = Mat_parameters[13]; # Lateral Confinement Parameter 1
const kc1 = Mat_parameters[12]; # Lateral Confinement Parameter 2
const n_c = Mat_parameters[10]; # Extract Compressive Exponential Parameter
const n_t = Mat_parameters[9]; # Extract Tensile Exponential Parameter

E, K_t, K_s = ConstitutiveLawMat_parameters(Connect, Mat_parameters); # Calculate Effective Young's Modulus and the Exponential Softening Parameter for Every Connection According to the Required Fracture Energy
Tet_Mass, Tet_Stiffness = TetrahedronStiffnessMatrix(TRI, Connect, Unique_Connections, Mat_parameters, Elements, E);


Δt, Tet_Mass,scale = MassScaling(Tet_Mass, Tet_Stiffness, TRI);

M, C = MassMatrix(gdl, TRI, Tet_Mass, Mat_parameters[15]);


