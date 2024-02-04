function BoundaryConditions(Boun, Boun1)
    # BOUN_VEL Construct Vector Containing Velocity of Restrained DOFs (Boundary Condition)
    #  The function builds the velocity vector for the Restrained/Imposed DOFs
    f_d = Boun1 .== 0; # Logical Indexing of Free DOFs
    r_d = f_d .== 0; # Logical Indexing of Restrained/Imposed DOFS
    num_fd = count(f_d); # Number of Free DOFs rd
    num_rd = count(r_d); # Number of Restrained/Imposed DOFs
    # boundary_vel = boun_vel( imp_vel,Boun,rd ); % Build Velocity Vector for Restrained/Imposed DOFs
    boundary_vel = Boun[r_d]
    return f_d, r_d, num_fd, num_rd, boundary_vel
end


