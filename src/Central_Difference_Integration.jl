# This file uses the pre-processed information from the previous scripts to solve the explicit dynamics problem through a Central Difference Integration scheme
#  Units mm/sec/N/MPa/ton
# Analysis Settings
function Central_Difference_Integration(t_final)

# Calculated Parameters
n_inc = Int64(ceil(t_final / Δt)); # Number of Total Increments
rampth = ones(n_inc); # Initialize Rampth Function used to apply the external imposed displacements gradually
rampth_steps = 500; # Desired number of steps to move from 0 load to full loading
rampth[1:rampth_steps] = collect(range(1/rampth_steps, step = 1/rampth_steps, stop = 1)); # Rampth over a certain number of steps

f_d, r_d, num_fd, num_rd, boundary_vel = BoundaryConditions(Boun, Boun1);
# Initialization of Global Vectors (Lots of them!)
# F = zeros(gdl,1); % Initialize the Internal Forces Vector
F_fd = zeros(num_fd); # Initialize the Free DOFs Internal Forces Vector
F_rd = zeros(num_rd); # Initialize the Restrained DOFs Internal Forces
R = zeros(gdl); # Initialize the External Forces Vector
R_fd = zeros(num_fd); # Initialize the Free DOFs external Forces Vector
R_rd = zeros(num_rd); # Initialize the Restrained DOFs external Forces

M_fd = M[f_d];
M_rd = M[r_d];
C_rd = C[r_d];
#
U = zeros(gdl); # Initialize Current Displacement Vector
U_fd = zeros(num_fd); # Initialize Free DOFs Displacement Vector
U_rd = zeros(num_rd); # Initialize Restrained DOFs Displacement Vector
# V_rd = zeros(num_rd,1); % Initialize Restrained DOFs Velocity Vector (Known!)
A_rd = zeros(num_rd); # Initialize Restrained DOFs Acceleration Vector (Known!)
# U_minus = zeros(gdl,1); % Initialize Previous Step Displacement Vector
U_minus_fd = zeros(num_fd); # Initialize Free DOFs Displacement Vector
U_minus_rd = zeros(num_rd); # Initialize Restrained DOFs Displacement Vector
U_plus = zeros(gdl); # Initialize Next Step Displacement Vector
U_plus_fd = zeros(num_fd); # Initialize Free DOFs Displacement Vector
U_plus_rd = zeros(num_rd); # Initialize Restrained DOFs Displacement Vector
V_plus_rd = zeros(num_rd); # Initialize Restrained DOFs Velocity Vector
#
displace = zeros(gdl,n_inc+1); # Initialize Matrix to Store Displacement Values
eps_ = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Deformation
eps_n = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Axial Deformation
eps_t = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Tangential Deformation
eps_vol = zeros(nel); # Initialize Vector Containing Volumetric Deformation
sigma_eff = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Effective Stress
sigma_N = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Axial Stress
sigma_T = zeros(nel,n_inc+1); # Initialize Matrix to Store Element Tangential Stress
# omega = zeros(nel,n_inc+1); % Initialize Matrix to Store Coupling Strain
qi = [Array{Float64}(undef,12) for i in 1:nel]; # Initialize Cell to Store Element's Internal Forces # qi = zeros(12,nel)
eps_pos_max = zeros(nel); # Initialize Vector to Store Maximum Strain Achieved by Element (Used to Detect Unloading)
eps_neg_max = zeros(nel); # Initialize Vector to Store Minimim Strain Achieved by Element (Used to Detect Unloading)
eps_n_max = zeros(nel); # Initialize Vector to Store Maximum Positive Axial Strain
eps_t_max = zeros(nel); # Initialize Vector to Store Maximum Shear Strain
internal = zeros(gdl,n_inc+1); # Initialize Matrix to Store the Internal Forces Vector
external = zeros(gdl,n_inc+1); # Initialize Matrix to Store the External Forces Vector
#
n_i = 1; # 1st Step of the Explicit Integration Requires some Extra Calculations for Trajectory Initialization
@inbounds for j = 1:num_fd; # Cycle over Free DOFs
   initial_vel_fd = 0;
   initial_acc_fd = 0;
   U_minus_fd[j] = U_fd[j] + initial_acc_fd/2*Δt^2-Δt*initial_vel_fd;
   R_tilda = R_fd[j]-F_fd[j] + U_fd[j]*2*M_fd[j]/Δt^2-U_minus_fd[j]*M_fd[j]/Δt^2;
   U_plus_fd[j] = (R_tilda*Δt^2)/M_fd[j];
end
@inbounds for j=1:num_rd; # Cycle over Restrained/ImposedDisp DOFs
    initial_vel = boundary_vel[j]*rampth[n_i];
    initial_acc = 0;
    U_minus_rd[j] = U_rd[j]+initial_acc/2*Δt^2-Δt*initial_vel;
    U_plus_rd[j] = initial_acc/2*Δt^2+initial_vel*Δt;
    R_rd[j] = F_rd[j]+ U_plus_rd[j]*(M_rd[j]/Δt^2+C_rd[j]/(2*Δt)) - U_rd[j]*M_rd[j]*2/Δt^2 + U_minus_rd[j]*(M_rd[j]/Δt^2-C_rd[j]/(2*Δt));
end
R[f_d] = R_fd; # Update External Forces for Free DOFs
R[r_d] = R_rd; # Update External Forces for Restrained DOFs
U_plus[f_d] = U_plus_fd; # Update U_plus Vector
U_plus[r_d] = U_plus_rd; # Update U_plus Vector
U_minus = copy(U); # Advance to Next Step (Current Displacement Becomes Previous Step Displacement)
U = copy(U_plus); # Advance to Next Step (Next-Step Displacement Becomes Current Displacement)
V_rd = boundary_vel .* rampth[n_i]; # Initialize Velocity Vector Equal to the Value Calculated for the Restrained/Imposed DOFs
F = zeros(gdl);
@inbounds for iel = 1:nel
    #id_1 = Connect[iel].ID1; # ID of 1st Node
    #id_2 = Connect[iel].ID2; # ID of 2nd Node
    area = Connect[iel].Area; # extract cross-sectional area
    length_ = Connect[iel].Length; # Length of the Current Connection
    el_disp = U[Connect[iel].DOFs]; # Element Displacement Vector [U_I ; U_J]
    qi[iel], eps_[iel,n_i:n_i+1], eps_n[iel,n_i:n_i+1], eps_t[iel,n_i:n_i+1], sigma_eff[iel,n_i:n_i+1], sigma_N[iel,n_i:n_i+1], sigma_T[iel,n_i:n_i+1] = Element_Response(
    eps_n[iel,n_i:n_i+1], eps_t[iel,n_i:n_i+1], eps_[iel,n_i:n_i+1], E[iel], sigma_N[iel,n_i:n_i+1], sigma_T[iel,n_i:n_i+1], sigma_eff[iel,n_i:n_i+1], eps_pos_max[iel], eps_neg_max[iel], eps_n_max[iel], eps_t_max[iel], K_t[iel], K_s[iel], eps_vol[iel], area, length_, el_disp, B[iel]);
    F[Connect[iel].DOFs] += qi[iel]; # Update Internal Forces Vector with the Value of Internal Forces Computed
end
internal[:,n_i+1] = F; # Store Force Values
external[:,n_i+1] = R; # Store Reactions at the Boundaries
displace[:,n_i+1] = U; # Store Displacement Values
@inbounds for j=1:nel
    eps_pos_max[j], eps_n_max[j], eps_neg_max[j], eps_t_max[j] = eps_update!(eps_[j,n_i+1],eps_n[j,n_i+1],eps_t[j,n_i+1],eps_pos_max[j],eps_neg_max[j],eps_n_max[j],eps_t_max[j]);
end
eps_vol = zeros(nel)
eps_vol_tet = zeros(n_tet)
for i=1:n_tet
    eps_vol_tet[i] = Tet_Vol_Def(Tet[i], reshape(U[Tet[i].ind],3,4)');
end
@inbounds for i=1:nel
    eps_vol[i] = Con_Vol_Def(eps_vol_tet, Unique_Connections_Elements[i], Connect[i].Area)
end
#fprintf('Step %d of %d Computed. \n', n_i, n_inc);
for n_i=2:n_inc
    F_fd = F[f_d];
    R_fd = R[f_d];
    U_fd = U[f_d];
    U_minus_fd = U_minus[f_d]; # Refresh Free Dof Displacement Vector
    @inbounds for j=1:num_fd; # Cycle over Free DOFs (Parallelized)
        R_tilda = R_fd[j]-F_fd[j]+U_fd[j]*2*M_fd[j]/Δt^2-U_minus_fd[j]*M_fd[j]/Δt^2;
        U_plus_fd[j] = (R_tilda*Δt^2)/M_fd[j];
    end
    U_rd = U[r_d]; # Refresh Restrained/Imposed DOFs Displacement Vector
    F_rd = F[r_d]; # Refresh Restrained/Imposed DOFs Internal Forces Vector
    U_minus_rd = U_minus[r_d];
    @inbounds for j=1:num_rd; # Cycle over Restrained/ImposedDisp DOFs (Parallelized)
        V_plus_rd[j] = boundary_vel[j]*rampth[n_i];
        A_rd[j] = (V_plus_rd[j]-V_rd[j])/Δt;
        U_plus_rd[j] = U_rd[j]+V_plus_rd[j]*Δt+A_rd[j]/2*Δt^2;
        R_rd[j] = F_rd[j]+ U_plus_rd[j]*(M_rd[j]/Δt^2+C_rd[j]/(2*Δt)) - U_rd[j]*M_rd[j]*2/Δt^2 + U_minus_rd[j]*(M_rd[j]/Δt^2-C_rd[j]/(2*Δt));
    end
R[f_d] = R_fd;
R[r_d] = R_rd;
U_plus[f_d] = U_plus_fd; # Rebuild U_plus Vector
U_plus[r_d] = U_plus_rd; # Rebuild U_plus Vector
U_minus = copy(U); # Advance to Next Step
U = copy(U_plus); # Advance to Next Step
V_rd = copy(V_plus_rd); # Advance to Next Step (Restrained DOFs Velocity)
F = zeros(gdl);
@inbounds for iel = 1:nel
    #id_1 = Connect[iel].ID1; # ID of 1st Node
    #id_2 = Connect[iel].ID2; # ID of 2nd Node
    area = Connect[iel].Area; # extract cross-sectional area
    length_ = Connect[iel].Length; # Length of the Current Connection
    el_disp = U[Connect[iel].DOFs]; # Element Displacement Vector [U_I ; U_J]
    qi[iel], eps_[iel,n_i:n_i+1], eps_n[iel,n_i:n_i+1], eps_t[iel,n_i:n_i+1], sigma_eff[iel,n_i:n_i+1], sigma_N[iel,n_i:n_i+1], sigma_T[iel,n_i:n_i+1] = Element_Response(
    eps_n[iel,n_i:n_i+1], eps_t[iel,n_i:n_i+1], eps_[iel,n_i:n_i+1], E[iel], sigma_N[iel,n_i:n_i+1], sigma_T[iel,n_i:n_i+1], sigma_eff[iel,n_i:n_i+1], eps_pos_max[iel], eps_neg_max[iel], eps_n_max[iel], eps_t_max[iel], K_t[iel], K_s[iel], eps_vol[iel], area, length_, el_disp, B[iel]);
    F[Connect[iel].DOFs] += qi[iel]; # Update Internal Forces Vector with the Value of Internal Forces Computed
end
internal[:,n_i+1] = F; # Store Force Values
external[:,n_i+1] = R; # Store Reactions at the Boundaries
displace[:,n_i+1] = U; # Store Displacement Values
@inbounds for j=1:nel
    eps_pos_max[j], eps_n_max[j], eps_neg_max[j], eps_t_max[j] = eps_update!(eps_[j,n_i+1],eps_n[j,n_i+1],eps_t[j,n_i+1],eps_pos_max[j],eps_neg_max[j],eps_n_max[j],eps_t_max[j]);
end
# [ eps_pos_max,eps_neg_max,eps_n_max,eps_t_max ] = eps_update_2( nel,eps(:,n),eps_n(:,n),eps_t(:,n),eps_pos_max,eps_neg_max,eps_n_max,eps_t_max ); # Refresh Maximum and Minimum Epsilon Reached Up to Current Step
eps_vol = zeros(nel)
eps_vol_tet = zeros(n_tet)
@inbounds for i=1:n_tet
    eps_vol_tet[i] = Tet_Vol_Def(Tet[i], reshape(U[Tet[i].ind],3,4)');
end
@inbounds for i=1:nel
    eps_vol[i] = Con_Vol_Def(eps_vol_tet, Unique_Connections_Elements[i], Connect[i].Area)
end
if n_i % 100 == 0
println("Computed Step #", n_i," of ", n_inc,  ".")
end
end
return displace, eps_, eps_n, eps_t, sigma_eff, sigma_N, sigma_T, internal, external
#
end


## t_final=0.2
#const n_inc = Int64(ceil(t_final / Δt)); # Number of Total Increments
displace, eps_, eps_n, eps_t, sigma_eff, sigma_N, sigma_T, internal, external = Central_Difference_Integration(t_final)
