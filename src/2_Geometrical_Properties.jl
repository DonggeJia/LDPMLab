# Evaluation of the volume and centroids of the Tethraedrons and the
# Aggregate Volumes
#
struct Tet_Structure
    Vert :: Array{Float64,2}                                                                                           # Positions of the Vertices
    V :: Float64                                                                                                       # Volume of the Tetrahedron
    ind :: Array{Int64,1}                                                                                              # Indices of the Vertices in the Global Displacement Vector # degree (x y z) of each node
end

function Geometrical_Props(Positions, TRI, Size)
    n_tet = size(TRI,1);                                                                                               # Number of tetrahedron
    edge_l = zeros(n_tet,6)                                                                                            # the 6 edges in each tetrahedron
    V = zeros(n_tet)                                                                                                   # tetrahedron volume 
    G = zeros(n_tet,3)                                                                                                 # Centroid Coordinates of each tetrahedron
    M_AB = zeros(n_tet,3)                                                                                              # Coordinates of the effective mid-point of the strut
    M_AC = zeros(n_tet,3)
    M_AD = zeros(n_tet,3)
    M_BC = zeros(n_tet,3)
    M_BD = zeros(n_tet,3)
    M_CD = zeros(n_tet,3)
    Vagg = zeros(n_tet,4)                                                                                              # Volume of aggregate at point A B C D
    G_A = zeros(n_tet,3)                                                                                               # Center of the aggregate at point A
    G_B = zeros(n_tet,3)
    G_C = zeros(n_tet,3)
    G_D = zeros(n_tet,3)
    G_ABCD = zeros(n_tet,3)                                                                                             # effective Center of  of tetrahedron
    G_ABC = zeros(n_tet,3)                                                                                              # effective Center of  of triangle ABC
    G_ABD = zeros(n_tet,3)
    G_ACD = zeros(n_tet,3)
    G_BCD = zeros(n_tet,3)
    V_eff = zeros(n_tet)
    A_eff = zeros(n_tet,4)

    Tet = Array{Tet_Structure}(undef,n_tet)

    @inbounds for i=1:n_tet
        ####
        #    Evaluation of the Volume and the Centroid of the Tetrahedrons
        ####
        A = Positions[TRI[i,1],:];                                                                                       # 1st vertex of the tethraedron
        B = Positions[TRI[i,2],:];                                                                                       # 2nd vertex of the tethraedron
        C = Positions[TRI[i,3],:];                                                                                       # 3rd vertex of the tethraedron
        D = Positions[TRI[i,4],:];                                                                                       # 4th vertex of the tethraedron
        b = B - A;                                                                                                       # 1st edge
        c = C - A;                                                                                                       # 2nd edge
        d = D - A;                                                                                                       # 3rd edge
        V[i] = abs(dot(b,cross(c,d)))/6;                                                                                 # Volume of the Tethrhedron
        G[i,:] = (A + B + C + D)/4;                                                                                      # Centroid Coordinates
        ind = [collect((6*TRI[i,1]-5):(6*TRI[i,1]-3));collect((6*TRI[i,2]-5):(6*TRI[i,2]-3));collect((6*TRI[i,3]-5):(6*TRI[i,3]-3));collect((6*TRI[i,4]-5):(6*TRI[i,4]-3))];
                                                                    # Indices of the Vertices in the Global Displacement Vector,  degree (x y z) of each node
        Tet[i] = Tet_Structure([A';B';C';D'], V[i], ind);

        ####
        #    The geometry of each edge
        ####
        # AB Segment
        l_AB = norm(A-B);                                                                                               # length of the edge
        l_AB_ = l_AB - Size[TRI[i,1]]/2 - Size[TRI[i,2]]/2;                                                             # effective length of the strut
        uv_ab = b/l_AB;                                                                                                 # Unit Vector on the line between A and B
        M_AB[i,:] = A + uv_ab*(Size[TRI[i,1]]/2 + l_AB_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # AC Segment
        l_AC = norm(A-C);                                                                                               # length of the edge
        l_AC_ = l_AC - Size[TRI[i,1]]/2 - Size[TRI[i,3]]/2;                                                             # effective length of the strut
        uv_ac = c/l_AC;                                                                                                 # Unit Vector on the line between A and C
        M_AC[i,:] = A + uv_ac*(Size[TRI[i,1]]/2 + l_AC_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # AD Segment
        l_AD = norm(A-D);                                                                                               # length of the edge
        l_AD_ = l_AD - Size[TRI[i,1]]/2 - Size[TRI[i,4]]/2;                                                             # effective length of the strut
        uv_ad = d/l_AD;                                                                                                 # Unit Vector on the line between A and D
        M_AD[i,:] = A + uv_ad*(Size[TRI[i,1]]/2 + l_AD_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # BC Segment
        l_BC = norm(B-C);                                                                                               # length of the edge
        l_BC_ = l_BC - Size[TRI[i,2]]/2 - Size[TRI[i,3]]/2;                                                             # effective length of the strut
        uv_bc = (C - B)/l_BC;                                                                                           # Unit Vector on the line between B and C
        M_BC[i,:] = B + uv_bc*(Size[TRI[i,2]]/2 + l_BC_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # BD Segment
        l_BD = norm(B-D);                                                                                               # length of the edge
        l_BD_ = l_BD - Size[TRI[i,2]]/2 - Size[TRI[i,4]]/2;                                                             # effective length of the strut
        uv_bd = (D - B)/l_BD;                                                                                           # Unit Vector on the line between B and D
        M_BD[i,:] = B + uv_bd*(Size[TRI[i,2]]/2 + l_BD_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # CD Segment
        l_CD = norm(D-C);                                                                                               # length of the edge
        l_CD_ = l_CD - Size[TRI[i,3]]/2 - Size[TRI[i,4]]/2;                                                             # effective length of the strut
        uv_cd = (D - C)/l_CD;                                                                                           # Unit Vector on the line between C and D
        M_CD[i,:] = C + uv_cd*(Size[TRI[i,3]]/2 + l_CD_/2 );                                                            # Coordinates of the effective mid-point of the strut
        # Matrix Containing the Effective Length of the Struts
        edge_l[i,:] = [l_AB_ l_AC_ l_AD_ l_BC_ l_BD_ l_CD_];


        #
        # Volume of the Spherical Tetrahedrons Generated by the Intersection
        # Between Aggregates and Paste

        #
        # Point A of the Tethraedron - to cal. centriod of the 1st Aggregate
        #
        bac = acos(dot(uv_ab,uv_ac));                                                                                   # Angle Between Vectors AB & AC
        bad = acos(dot(uv_ab,uv_ad));                                                                                   # Angle Between Vectors AB & AD
        cad = acos(dot(uv_ac,uv_ad));                                                                                   # Angle Between Vectors AC & AD
        ang1_a = acos((cos(bac) - cos(bad)*cos(cad))/(sin(bad)*sin(cad)));                                              # dihedral angle "AD" Angle on the sphere obtained by the sine rule
        ang2_a = acos((cos(bad) - cos(bac)*cos(cad))/(sin(bac)*sin(cad)));                                              # dihedral angle "AC" Angle on the sphere obtained by the sine rule
        ang3_a = acos((cos(cad) - cos(bac)*cos(bad))/(sin(bac)*sin(bad)));                                              # dihedral angle "AB" Angle on the sphere obtained by the sine rule
        d_a=3*V[i] /(abs(norm(cross(C-B,D-B)))*0.5)                                                                     # Distances Between A and the Plane Containing BCD      
        r_a = min(d_a,Size[TRI[i,1]]/2);                                                                                # This is a Check to Prevent an Aggregate from Crossing one of the Tetrahedron Opposite Faces, Leading to Wrong Results
        Vagg[i,1] = ((ang1_a + ang2_a + ang3_a - pi)*r_a^3)/3;                                                          # Volume of aggregate A,  V = 1/3 * (spherical excess) * R^3
        vec1_a = (cross(uv_ac,uv_ad))/(norm(cross(uv_ac,uv_ad)))*sign(dot(cross(uv_ac,uv_ad),uv_ab));                   # (AC x AD)/|(AC x AD)| normal unit vector
        vec2_a = (cross(uv_ab,uv_ad))/(norm(cross(uv_ab,uv_ad)))*sign(dot(cross(uv_ab,uv_ad),uv_ac));                   # (AB x AD)/|(AB x AD)| normal unit vector
        vec3_a = (cross(uv_ab,uv_ac))/(norm(cross(uv_ab,uv_ac)))*sign(dot(cross(uv_ab,uv_ac),uv_ad));                   # (AB x AC)/|(AB x AC)| normal unit vector
        G_A[i,:] = A + (3/4*r_a/(2*(ang1_a + ang2_a + ang3_a - pi)))*(vec1_a*cad + vec2_a*bad + vec3_a*bac);            # Center of Mass of Aggregate "A"
        
        #
        # Point B of the Tethraedron - cal. centriod of the 2nd Aggregate
        #
        abc = acos(dot(-uv_ab,uv_bc));                                                                                   # Angle Between Vectors BA & BC
        abd = acos(dot(-uv_ab,uv_bd));                                                                                   # Angle Between Vectors BA & BD
        cbd = acos(dot(uv_bc,uv_bd));                                                                                    # Angle Between Vectors BC & BD
        ang1_b = acos((cos(abc) - cos(abd)*cos(cbd))/(sin(abd)*sin(cbd)));                                               # "D" Angle on the sphere obtained by the sine rule
        ang2_b = acos((cos(abd) - cos(abc)*cos(cbd))/(sin(abc)*sin(cbd)));                                               # "C" Angle on the sphere obtained by the sine rule
        ang3_b = acos((cos(cbd) - cos(abc)*cos(abd))/(sin(abc)*sin(abd)));                                               # "A" Angle on the sphere obtained by the sine rule
        d_b= 3* V[i] /(abs(norm(cross(C-A,D-A)))*0.5)                                                                    # Distances Between B and the Plane Containing ACD
        r_b = min(d_b,Size[TRI[i,2]]/2);                                                                                 # This is a Check to Prevent an Aggregate from Crossing one of the Tetrahedron Opposite Faces, Leading to Wrong Results
        Vagg[i,2] = ((ang1_b + ang2_b + ang3_b - pi)*r_b^3)/3;                                                           # V = 1/3 * (spherical excess) * R^3
        vec1_b = (cross(uv_bc,uv_bd))/(norm(cross(uv_bc,uv_bd)))*sign(dot(cross(uv_bc,uv_bd),-uv_ab));                   # (BC x BD)/|(BC x BD)|
        vec2_b = (cross(-uv_ab,uv_bd))/(norm(cross(-uv_ab,uv_bd)))*sign(dot(cross(-uv_ab,uv_bd),uv_bc));                 # (BA x BD)/|(BA x BD)|
        vec3_b = (cross(-uv_ab,uv_bc))/(norm(cross(-uv_ab,uv_bc)))*sign(dot(cross(-uv_ab,uv_bc),uv_bd));                 # (BA x BC)/|(BA x BC)|
        G_B[i,:] = B + (3/4*r_b/(2*(ang1_b + ang2_b + ang3_b - pi)))*(vec1_b*cbd + vec2_b*abd + vec3_b*abc);             # Center of Mass of Aggregate "B"
        
        #
        # Point C of the Tethraedron - 3rd Aggregate
        #
        acb = acos(dot(-uv_ac,-uv_bc));                                                                                   # Angle Between Vectors CA & CB
        acd = acos(dot(-uv_ac,uv_cd));                                                                                    # Angle Between Vectors CA & CD
        bcd = acos(dot(-uv_bc,uv_cd));                                                                                    # Angle Between Vectors CB & CD
        ang1_c = acos((cos(acb) - cos(acd)*cos(bcd))/(sin(acd)*sin(bcd)));                                                # "D" Angle on the sphere obtained by the sine rule
        ang2_c = acos((cos(acd) - cos(acb)*cos(bcd))/(sin(acb)*sin(bcd)));                                                # "B" Angle on the sphere obtained by the sine rule
        ang3_c = acos((cos(bcd) - cos(acb)*cos(acd))/(sin(acb)*sin(acd)));                                                # "A" Angle on the sphere obtained by the sine rule
        d_c=3* V[i] /(abs(norm(cross(B-A,D-A)))*0.5)                                                                      # Distances Between C and the Plane Containing ABD
        r_c = min(d_c,(Size[TRI[i,3]]/2))
        Vagg[i,3] = ((ang1_c + ang2_c + ang3_c - pi)*r_c^3)/3;                                                            # V = 1/3 * (spherical excess) * R^3
        vec1_c = (cross(-uv_bc,uv_cd))/(norm(cross(-uv_bc,uv_cd)))*sign(dot(cross(-uv_bc,uv_cd),-uv_ac));                 # (CB x CD)/|(CB x CD)|
        vec2_c = (cross(-uv_ac,uv_cd))/(norm(cross(-uv_ac,uv_cd)))*sign(dot(cross(-uv_ac,uv_cd),-uv_bc));                 # (CA x CD)/|(CA x CD)|
        vec3_c = (cross(-uv_ac,-uv_bc))/(norm(cross(-uv_ac,-uv_bc)))*sign(dot(cross(-uv_ac,-uv_bc),uv_cd));               # (CA x CB)/|(CA x CB)|
        G_C[i,:] = C + (3/4*r_c/(2*(ang1_c + ang2_c + ang3_c - pi)))*(vec1_c*bcd + vec2_c*acd + vec3_c*acb);              # Center of Mass of Aggregate "C"
        
        #
        # Point D of the Tethraedron - 4th Aggregate
        #
        adb = acos(dot(-uv_ad,-uv_bd));                                                                                    # Angle Between Vectors DA & DB
        adc = acos(dot(-uv_ad,-uv_cd));                                                                                    # Angle Between Vectors DA & DC
        bdc = acos(dot(-uv_bd,-uv_cd));                                                                                    # Angle Between Vectors DB & DC
        ang1_d = acos((cos(adb) - cos(adc)*cos(bdc))/(sin(adc)*sin(bdc)));                                                 # "C" Angle on the sphere obtained by the sine rule
        ang2_d = acos((cos(adc) - cos(adb)*cos(bdc))/(sin(adb)*sin(bdc)));                                                 # "B" Angle on the sphere obtained by the sine rule
        ang3_d = acos((cos(bdc) - cos(adb)*cos(adc))/(sin(adb)*sin(adc)));                                                 # "A" Angle on the sphere obtained by the sine rule
        d_d= 3* V[i] /(abs(norm(cross(B-A,C-A)))*0.5)                                                                      # Distances Between D and the Plane Containing ABC
        r_d = min(d_d,(Size[TRI[i,4]]/2))
        Vagg[i,4] = ((ang1_d + ang2_d + ang3_d - pi)*r_d^3)/3;                                                             # V = 1/3 * (spherical excess) * R^3
        vec1_d = (cross(-uv_bd,-uv_cd))/(norm(cross(-uv_bd,-uv_cd)))*sign(dot(cross(-uv_bd,-uv_cd),-uv_ad));               # (DB x DC)/|(DB x DC)|
        vec2_d = (cross(-uv_ad,-uv_cd))/(norm(cross(-uv_ad,-uv_cd)))*sign(dot(cross(-uv_ad,-uv_cd),-uv_bd));               # (DA x DC)/|(DA x DC)|
        vec3_d = (cross(-uv_ad,-uv_bd))/(norm(cross(-uv_ad,-uv_bd)))*sign(dot(cross(-uv_ad,-uv_bd),-uv_cd));               # (DA x DB)/|(DA x DB)|
        G_D[i,:] = D + (3/4*r_d/(2*(ang1_d + ang2_d + ang3_d - pi)))*(vec1_d*bdc + vec2_d*adc + vec3_d*adb);               # Center of Mass of Aggregate "D"
        #
        # Effective Center of Mass of the Volume Obtained by Subtracting From
        # the Volume of the Tetrahedron the Counterparts of Aggregate Volumes
        #
        V_eff[i] = V[i] - Vagg[i,1] - Vagg[i,2] - Vagg[i,3] - Vagg[i,4];                                                   # Effective Volume
        if V_eff[i] > 0
            G_ABCD[i,:] = (V[i]*G[i,:] - Vagg[i,1]*G_A[i,:] - Vagg[i,2]*G_B[i,:] - Vagg[i,3]*G_C[i,:] - Vagg[i,4]*G_D[i,:])/(V_eff[i]); # Effective Center of Mass
        else()
            G_ABCD[i,:] = G[i,:]
        end
        
        #
        # 2D
        #
        # Evaluation of the Properties of the 2D Faces of the Tetrahedron
        #

        #
        # ABC Triangle
        #
        A_ABC = 1/2*norm(cross(C-A,B-A));                                                                                    # Area of the Triangle
        G_ABC_ = (A + B + C)/3;                                                                                              # Centroid of the Triangle
        # Center of Mass of Circular Sector "A"
        d_a_abc = norm(cross(A-B,A-C))/norm(B-C);                                                                            # Distance Between Point A and the Line for BC
        r_a_abc = min(d_a_abc,(Size[TRI[i,1]]/2));                                                                           # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABC_A = (r_a_abc^2)*bac/2;                                                                                           # Area of the Circular Sector Defined by the Intersection of Aggregate A with ABC Triangle
        G_A_ABC = A + ((uv_ab + uv_ac)/norm(uv_ab + uv_ac))*2/3*r_a_abc*(sin(bac/2)/(bac/2));                                # Center of Mass of Circular Sector "A"
        # Center of Mass of Circular Sector "B"
        d_b_abc = norm(cross(B-A,B-C))/norm(A-C);                                                                            # Distance Between Point B and the Line for AC
        r_b_abc = min(d_b_abc,(Size[TRI[i,2]]/2));                                                                           # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABC_B = (r_b_abc^2)*abc/2;                                                                                           # Area of Circular Sector "B"
        G_B_ABC = B + ((-uv_ab + uv_bc)/norm(-uv_ab + uv_bc))*2/3*r_b_abc*(sin(abc/2)/(abc/2));                              # Center of Mass of Circular Sector "B"=4/3*r*sin(0.5angle)/angle
        # Center of Mass of Circular Sector "C"        
        d_c_abc = norm(cross(C-A,C-B))/norm(A-B);                                                                            # Distance Between Point C and the Line for AB
        r_c_abc = min(d_c_abc,(Size[TRI[i,3]]/2));                                                                           # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABC_C = (r_c_abc^2)*acb/2;                                                                                           # Area of Circular Sector "C"
        #cal. Effective center, effective Area and Effective Center of Mass of the ABC        
        G_C_ABC = C + ((-uv_ac - uv_bc)/norm((-uv_ac - uv_bc))*2/3*r_c_abc*sin(acb/2)/(acb/2));                              # Center of Mass of Circular Sector "C"
        A_eff[i,1] = A_ABC - ABC_A -ABC_B - ABC_C;                                                                           # Effective Area
        if A_eff[i,1] > 0
            G_ABC[i,:] = (A_ABC*G_ABC_ - ABC_A*G_A_ABC - ABC_B*G_B_ABC - ABC_C*G_C_ABC)/(A_eff[i,1]);                        # Effective Center of Mass of the ABC Triangle
        else()
            G_ABC[i,:] = G_ABC_
        end

        #
        # ABD Triangle
        #
        A_ABD = 1/2*norm(cross(B-A,D-A));                                                                                     # Area of the Triangle
        G_ABD_ = (A + B + D)/3;                                                                                               # Centroid of the Triangle
        # Center of Mass of Circular Sector "A"
        d_a_abd = norm(cross(A-B,A-D))/norm(B-D);                                                                             # Distance Between Point A and the Line for BD
        r_a_abd = min(d_a_abd,(Size[TRI[i,1]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABD_A = (r_a_abd^2)*bad/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate A with ABD Triangle
        G_A_ABD = A + ((uv_ab + uv_ad)/norm(uv_ab + uv_ad))*2/3*r_a_abd*(sin(bad/2)/(bad/2));                                 # Center of Mass of Circular Sector "A"
        # Center of Mass of Circular Sector "B"
        d_b_abd = norm(cross(B-A,B-D))/norm(A-D);                                                                             # Distance Between Point B and the Line for AD
        r_b_abd = min(d_b_abd,(Size[TRI[i,2]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABD_B = (r_b_abd^2)*abd/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate B with ABD Triangle
        G_B_ABD = B + ((-uv_ab + uv_bd)/norm(-uv_ab + uv_bd))*2/3*r_b_abd*(sin(abd/2)/(abd/2));                               # Center of Mass of Circular Sector "B"
        # Center of Mass of Circular Sector "D"
        d_d_abd = norm(cross(D-A,D-B))/norm(A-B);                                                                             # Distance Between Point D and the Line for AB
        r_d_abd = min(d_d_abd,(Size[TRI[i,4]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ABD_D = (r_d_abd^2)*adb/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate D with ABD Triangle
        #cal. Effective center, effective Area and Effective Center of Mass of the ABC 
        G_D_ABD = D + ((-uv_ad - uv_bd)/norm(-uv_ad - uv_bd))*2/3*r_d_abd*(sin(adb/2)/(adb/2));                               # Center of Mass of Circular Sector "D"
        A_eff[i,2] = A_ABD - ABD_A - ABD_B - ABD_D;                                                                           # Effective Area
        if A_eff[i,2] > 0
            G_ABD[i,:] = (A_ABD*G_ABD_ - ABD_A*G_A_ABD - ABD_B*G_B_ABD - ABD_D*G_D_ABD)/(A_eff[i,2]);                         # Effective Center of Mass of the ABD Triangle
        else()
            G_ABD[i,:] = G_ABD_
        end

        #
        # ACD Triangle
        #
        A_ACD = 1/2*norm(cross(C-A,D-A));                                                                                     # Area of the Triangle
        G_ACD_ = (A + C + D)/3;                                                                                               # Centroid of the Triangle
        # Center of Mass of Circular Sector "A"
        d_a_acd = norm(cross(A-C,A-D))/norm(C-D);                                                                             # Distance Between Point A and the Line for CD
        r_a_acd = min(d_a_acd,(Size[TRI[i,1]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ACD_A = (r_a_acd^2)*cad/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate A with ACD Triangle
        G_A_ACD = A + ((uv_ac + uv_ad)/norm(uv_ac + uv_ad))*2/3*r_a_acd*(sin(cad/2)/(cad/2));                                 # Center of Mass of Circular Sector "A"
        # Center of Mass of Circular Sector "C"
        d_c_acd = norm(cross(C-A,C-D))/norm(A-D);                                                                             # Distance Between Point C and the Line for AD
        r_c_acd = min(d_c_acd,(Size[TRI[i,3]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ACD_C = (r_c_acd^2)*acd/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate C with ACD Triangle
        G_C_ACD = C + ((-uv_ac + uv_cd)/norm(-uv_ac + uv_cd))*2/3*r_c_acd*(sin(acd/2)/(acd/2));                               # Center of Mass of Circular Sector "C"
        # Center of Mass of Circular Sector "D"
        d_d_acd = norm(cross(D-A,D-C))/norm(A-C);                                                                             # Distance Between Point D and the Line for AC
        r_d_acd = min(d_d_acd,(Size[TRI[i,4]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        ACD_D = (r_d_acd^2)*adc/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate D with ACD Triangle
        #cal. Effective center, effective Area and Effective Center of Mass of the ABC        
        G_D_ADC = D + ((-uv_ad - uv_cd)/norm(-uv_ad - uv_cd))*2/3*r_d_acd*(sin(adc/2)/(adc/2));                               # Center of Mass of Circular Sector "D"
        A_eff[i,3] = A_ACD - ACD_A - ACD_C - ACD_D;                                                                           # Effective Area
        if A_eff[i,3] > 0
            G_ACD[i,:] = (A_ACD*G_ACD_ - ACD_A*G_A_ACD - ACD_C*G_C_ACD - ACD_D*G_D_ADC)/(A_eff[i,3]);                         # Effective Center of Mass of the ACD Triangle
        else()
            G_ACD[i,:] = G_ACD_
        end

        #
        # BCD Triangle
        #
        A_BCD = 1/2*norm(cross(C-B,D-B));                                                                                     # Area of the Triangle
        G_BCD_ = (B + C + D)/3;                                                                                               # Centroid of the Triangle
        # Center of Mass of Circular Sector "B"
        d_b_bcd = norm(cross(B-C,B-D))/norm(C-D);                                                                             # Distance Between Point B and the Line for CD
        r_b_bcd = min(d_b_bcd,(Size[TRI[i,2]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        BCD_B = (r_b_bcd^2)*cbd/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate B with BCD Triangle
        G_B_BCD = B + ((uv_bc + uv_bd)/norm(uv_bc + uv_bd))*2/3*r_b_bcd*(sin(cbd/2)/(cbd/2));                                 # Center of Mass of Circular Sector "B"
        # Center of Mass of Circular Sector "C"
        d_c_bcd = norm(cross(C-B,C-D))/norm(B-D);                                                                             # Distance Between Point C and the Line for BD
        r_c_bcd = min(d_c_bcd,(Size[TRI[i,3]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        BCD_C = (r_c_bcd^2)*bcd/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate C with BCD Triangle
        G_C_BCD = C +((-uv_bc + uv_cd)/norm(-uv_bc + uv_cd))*2/3*r_c_bcd*(sin(bcd/2)/(bcd/2));                                # Center of Mass of Circular Sector "C"
        # Center of Mass of Circular Sector "D"
        d_d_bcd = norm(cross(D-B,D-C))/norm(B-C);                                                                             # Distance Between Point D and the Line for BC
        r_d_bcd = min(d_d_bcd,(Size[TRI[i,4]]/2));                                                                            # Check to Prevent the Aggregate to Cross the Opposite Edge, Leading to Wrong Results
        BCD_D = (r_d_bcd^2)*bdc/2;                                                                                            # Area of the Circular Sector Defined by the Intersection of Aggregate D with BCD Triangle
        #cal. Effective center, effective Area and Effective Center of Mass of the ABC         
        G_D_BCD = D + ((-uv_bd - uv_cd)/norm(-uv_bd - uv_cd))*2/3*r_d_bcd*(sin(bdc/2)/(bdc/2));                                                                                      # Center of Mass of Circular Sector "D"
        A_eff[i,4] = A_BCD - BCD_B - BCD_C - BCD_D;                                                                           # Effective Area
        if A_eff[i,4] > 0
            G_BCD[i,:] = (A_BCD*G_BCD_ - BCD_B*G_B_BCD - BCD_C*G_C_BCD - BCD_D*G_D_BCD)/(A_eff[i,4]);                         # Effective Center of Mass of the BCD Triangle
        else()
            G_BCD[i,:] = G_BCD_
        end

    end
        
    
    Gravity_final = hcat(G_ABCD, M_AB, M_AC, M_AD, M_BC, M_BD, M_CD, G_ABC, G_ABD, G_ACD, G_BCD)
    
    ## effective V and A
    if minimum(V_eff) < 0
        warning("Minimum Effective Volume is Negative!")
    end
    if min(minimum(A_eff)) < 0
        warning("Minimum Effective Area is Negative!")
    end


    return Gravity_final, n_tet, Tet

end


Gravity_final, n_tet, Tet = Geometrical_Props(Positions, TRI, Size);
