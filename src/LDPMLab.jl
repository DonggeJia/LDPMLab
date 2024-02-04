
module LDPMLab

# import all packages needed
using JLD2
using DelimitedFiles
using Statistics
using CSV
using Tables
using TetGen
#using Plots
using LinearAlgebra
using DataFrames
using Printf
using PlotlyJS
using Format

mutable struct ldpm
    geometry_parameters::Array{Float64,1}                                                 # Projected Elements composing the connection
    #Non_Projected::Array{Float64,2}                                             # Non-projected Elements composing the connection
    mechanical_parameters::Array{Float64,1}                                                                  # ID of first node of the strut
    #delta_t::Float64                                                                  # ID of second node of the strut
    #Tet::Vector{Int64}                                                          # ID of tetrahedron
    #Area::Vector{Float64}                                                       # Area of each triangle (collated in a vector)
end
LDPM = ldpm([0, 0.0], [0, 0.0])

mutable struct ldpmwsteel
    geometry_parameters::Array{Float64,1}                                                 # Projected Elements composing the connection
    #Non_Projected::Array{Float64,2}                                             # Non-projected Elements composing the connection
    steel_layout::Array{Float64,1}
    mechanical_parameters::Array{Float64,1}                                                                  # ID of first node of the strut
    #delta_t::Float64                                                                  # ID of second node of the strut
    #Tet::Vector{Int64}                                                          # ID of tetrahedron
    #Area::Vector{Float64}                                                       # Area of each triangle (collated in a vector)
end
LDPM_w_steel = ldpmwsteel([0, 0.0], [0, 0.0], [0, 0.0])


function particle_distribution(model_name, particle_save="Yes", particle_dirc_and_name="LDPM_particle_distribution")
    if typeof(model_name) == ldpm
        include("0.0 generate particles.jl")
        include("0.1 particle position.jl")
    elseif typeof(model_name) == ldpmwsteel
        include("0.0 generate particles with steel.jl")
        include("0.1 particle position with steel.jl")
    end
    if particle_save == "Yes"
        @save "$particle_dirc_and_name.jld2" pointsresult pointsdiameterfinal
    end
    d1 = collect(d0:0.05:da)
    cdf(d) = ((1 - (d0 / d)^q)) / ((1 - (d0 / da)^q))
    ll = plot(scatter(x=d1, y=cdf.(d1), mode="lines"))
    d1 = collect(0:0.05:da)
    F(d) = (d / da)^nF
    PlotlyJS.add_trace!(ll, (scatter(x=d1, y=F.(d1), mode="lines")))
    PlotlyJS.add_trace!(ll, scatter(x=diameter, y=FF, mode="markers"))#u-----
    display(ll)
end

function load_particle_distribution(particle_dirc_and_name="LDPM_particle_distribution")
    @load "$particle_dirc_and_name.jld2" pointsresult pointsdiameterfinal
end
function Meshing(model_name, mesh_plot="Yes", mesh_dirc_and_name="LDPM_mesh_facets")
    if typeof(model_name) == ldpm
        include("1_Random_Meshing_Delaunay.jl")
        include("2_Geometrical_Properties.jl") # Evaluate Geometrical Properties of Tetrahedra
        include("3_Projected_Areas.jl") # Compute Areas of Triangles Composing Every Connection
        include("4_Connectivity_Matrix.jl") # Evaluate Connections Between Particles
        if mesh_plot == "Yes"
            include("4.5 Mesh plot.jl")
        end
    elseif typeof(model_name) == ldpmwsteel
        include("1_Random_Meshing_Delaunay with_steel.jl")
        include("2_Geometrical_Properties.jl") # Evaluate Geometrical Properties of Tetrahedra
        include("3_Projected_Areas.jl") # Compute Areas of Triangles Composing Every Connection
        include("4_Connectivity_Matrix_steel.jl") # Evaluate Connections Between Particles
        if mesh_plot == "Yes"
            include("4.5 Mesh plot.jl")
            include("4.6 steel_plot.jl")
        end
    end
end



function Boundary_setting(fixed_region=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], loaded_region=[[[0 10; 0 200; 65 70], [3, 4, 5], [-0.2, -0.1, 0.1]], [[100 110; 0 200; 65 70], [3, 4, 5], [-0.2, -0.1, 0.1]], [[190 200; 0 200; 65 70], [3, 4, 5], [-0.2, -0.1, 0.1]]], plot_boundary="Yes") #unit of velocity mm/sec
    include("5_Boundary.jl")
end

function Solutions(model_name, scale_delata_time=1.0, t_final=0.8) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
    Δt = Δt * scale_delata_time
    if typeof(model_name) == ldpm
        include("6_Material_Properties.jl") # Calculate Stable Time Step
        include("11_Boundary_Conditions.jl") # Function to Calculate Boundary Velocities
        include("12_Volumetric_Deformation.jl") # Function to Update Volumetric Deformation of Lattice Elements
        include("13_Epsilon_Update.jl") # Function to Update Maximum Deformation (Used in Constitutive Law)
        include("14_Element_Response.jl") # Function Used to Obtain Element Forces
        include("15_Constitutive_Law.jl") # Constitutive Law
        #include("15_Constitutive_Law steel and bond.jl") # Constitutive Law
        include("Central_Difference_Integration.jl")
    elseif typeof(model_name) == ldpmwsteel
        include("6_Material_Properties with steel.jl") # Calculate Stable Time Step
        include("11_Boundary_Conditions.jl") # Function to Calculate Boundary Velocities
        include("12_Volumetric_Deformation.jl") # Function to Update Volumetric Deformation of Lattice Elements
        include("13_Epsilon_Update.jl") # Function to Update Maximum Deformation (Used in Constitutive Law)
        include("14_Element_Response.jl") # Function Used to Obtain Element Forces
        include("15_Constitutive_Law.jl") # Constitutive Law
        include("15_Constitutive_Law steel and bond.jl") # Constitutive Law
        include("Central_Difference_Integration_with steel.jl")
    end
end

function post_process(model_name, relative_time_of_cracking=[0.2, 0.4, 0.5], crack_plot_dirc_and_name="cracking pattern", output_displacement_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], output_load_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], step_interval=300, load_dis_out_name="200*200*70 deck", plot_dis_load_region="Yes")
    if typeof(model_name) == ldpm
        include("vtk_cracks_ non_projected.jl") # Calculate Stable Time Step

    elseif typeof(model_name) == ldpmwsteel
        include("vtk_cracks_ non_projected steel.jl") # Calculate Stable Time Step
    end
    include("out_load_displacement.jl")
end

export load_particle_distribution, particle_distribution, Meshing, Boundary_setting, Solutions, post_process,
    LDPM, LDPM_w_steel
end

# #!use of package, units: mm, N, Mpa
# #dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp
# LDPM.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1]
# #dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp, height_S, diameter_S
# LDPM_w_steel.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1, 20.0, 16.0]
# #traverse_distribution 
# LDPM_w_steel.steel_layout = [50, 100, 150]

# particle_distribution(; model_name, Geometry_save="Yes", Geometry_dirc_and_name="LDPM_geometry")
# load_particle_distribution(particle_dirc_and_name="LDPM_particle_distribution")
# Meshing(model_name, mesh_plot="Yes", mesh_dirc_and_name="LDPM_mesh_facets")
# Boundary_setting(fixed_region=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], loaded_region=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], plot_boundary="Yes")

# 45000.0 # E_m -> Initial Elastic Modulus for the Matrix [MPa]
# 45000.0 # E_a -> Initial Elastic Modulus for the Aggregates [MPa]
# 3.0 # sigma_t -> Tension Limit [MPa]
# -50.0 # sigma_c -> Compression Limit [MPa]
# 10.0 # sigma_s -> Shear Limit [MPa]
# 0.07 # G_t -> Mode I Fracture Energy [N/mm]
# 0.35 # G_s -> Mode II Fracture Energy [N/mm]
# 0.25 # alpha -> Tangential-to-Normal Stiffness Ratio (Controls Poisson's Effect) [-]
# 2.0 # n_t -> Exponent that Controls Transition of Softening Parameter [-]
# 0.8 # n_c -> Exponent that Controls Transition of Hardening Parameter [-]
# 2.5e-6 # rho -> Material's Mass Density [ton/mm^3]
# 1.0 # kc1 -> Volumetric Compression Parameter [-]
# 5.0 # kc2 -> Volumetric Compression Parameter [-]
# 11250.0 # K_c -> Parameter that Governs the Post-Peak Behavior in Compression [-] 
# 0.0 # zeta -> damping coefficient
# LDPM.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0]
# #add steel mechanical parameters
# #E_steel = 1.96*10.0^5
# #fyi = 500
# #fu = 650 #Mpa
# #epsi_sh = 0.02
# #Esh = 833.33 #Mpa
# LDPM_w_steel.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]

# Solutions(model_name, scale_delata_time=1.0, t_final=0.05) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
# post_process(model_name, relative_time_of_cracking=[0.2, 0.4, 0.5], crack_plot_dirc_and_name="cracking pattern", output_displacement_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], output_load_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], step_interval=300, load_dis_out_name="200*200*70 deck", plot_dis_load_region="Yes")
