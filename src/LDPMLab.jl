
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
using Distances
using Colors
using GeometryBasics

mutable struct ldpm
    geometry_parameters::Array{Float64,1}                                                 # Projected Elements composing the connection
    #Non_Projected::Array{Float64,2}                                             # Non-projected Elements composing the connection
    mechanical_parameters::Array{Float64,1}                                                                  # ID of first node of the strut
    #delta_t::Float64                                                                  # ID of second node of the strut
    #Tet::Vector{Int64}                                                          # ID of tetrahedron
    #Area::Vector{Float64}                                                       # Area of each triangle (collated in a vector)
end

LDPM = ldpm([0, 0.0], [0, 0.0])

mutable struct ldpmbarr
    geometry_parameters::Array{Float64,1}                                                 # Projected Elements composing the connection
    #Non_Projected::Array{Float64,2}                                             # Non-projected Elements composing the connection
    steel_layout::Array{Float64,1}
    mechanical_parameters::Array{Float64,1}                                                                  # ID of first node of the strut
    #delta_t::Float64                                                                  # ID of second node of the strut
    #Tet::Vector{Int64}                                                          # ID of tetrahedron
    #Area::Vector{Float64}                                                       # Area of each triangle (collated in a vector)
end
LDPM_bar_reforced = ldpmbarr([0, 0.0], [0, 0.0], [0, 0.0])
# 4,5,9,10,11 are deleted
#then va>4 6>5 7>6 8>7 12>8 13>9 14>10 
LDPM.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1]
#dimen1, dimen2, dimen3, da, d0, nF, magnifyp, height_S, diameter_S
LDPM_bar_reforced.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1, 20.0, 16.0]
#traverse_distribution 
LDPM_bar_reforced.steel_layout = [50, 100, 150]

LDPM.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0]

LDPM_bar_reforced.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 0.02, 833.33]


const Mat_parameters = LDPM_bar_reforced.mechanical_parameters # [E_m E_a σ_t σ_c σ_s G_t G_s α n_t n_c ρ kc_1 kc_2 K_c ζ]
const k1 = -Mat_parameters[3] * Mat_parameters[4] / Mat_parameters[5]^2 # Stress Space Elliptic Boundary Parameter
const k2 = -Mat_parameters[3] * Mat_parameters[4] # Stress Space Elliptic Boundary Parameter
const sigma_t = Mat_parameters[3]# Extract Tension Limit from Mat_parameters Vector
const sigma_c = Mat_parameters[4] # Extract Compression Limit from Mat_parameters Vector
const sigma_s = Mat_parameters[5] # Extreact Shear Limit from Mat_parameters Vector
const alpha = Mat_parameters[8] # Extract Coupling Parameter
const K_c = Mat_parameters[14] # Extract Compressive Exponential Initial Slop Parameter
const kc2 = Mat_parameters[13] # Lateral Confinement Parameter 1
const kc1 = Mat_parameters[12] # Lateral Confinement Parameter 2
const n_c = Mat_parameters[10] # Extract Compressive Exponential Parameter
const n_t = Mat_parameters[9] # Extract Tensile Exponential Parameter

const E_steel = LDPM_bar_reforced.mechanical_parameters[end-3]
const fyi = LDPM_bar_reforced.mechanical_parameters[end-2]
const epsi_sh = LDPM_bar_reforced.mechanical_parameters[end-1]
const Esh = LDPM_bar_reforced.mechanical_parameters[end] #Mpa


"""
    particle_distribution(model_name, particle_save="Yes", particle_dirc_and_name="LDPM_particle_distribution")

Excute particle distribution in material volume and save the data of particle distribution.

### Input

- `model_name` -- symbol describing the model type used
   - `:LDPM` -- uses LDPM
   - `:LDPM_bar_reforced` -- uses rank1 update
- `particle_save` -- symbol describing wether saving particle distribution
    - `:Yes` -- save a JLD2 file with particle coordinates and their corresponding diameters. By default, `:Yes` is used. Other strings except `Yes` means not storing the current particle distribution.
    Saving particle distribution is encouraged so that the following solutions can use the same model mesh, excluding the influence of meshing on solutions.
- `particle_dirc_and_name` -- a string indicating the filefolder and filename saving particle distribution. By default, `LDPM_particle_distribution`  is used. This argument can be `D:/Juliafiles/LDPM_particle_distribution` as well.

"""
function Particle_distribution(model_name, particle_save="Yes", particle_dirc_and_name="LDPM_particle_distribution")
    cd(@__DIR__)
    if typeof(model_name) == ldpm
        include("0.0 generate particles.jl")
        include("0.1 particle position.jl")
    elseif typeof(model_name) == ldpmbarr
        include("0.0 generate particles with steel.jl")
        include("0.1 particle position with steel.jl")
    end
    if particle_save == "Yes"
        @save "$particle_dirc_and_name.jld2"
    end
    d1 = collect(d0:0.05:da)
    cdf(d) = ((1 - (d0 / d)^q)) / ((1 - (d0 / da)^q))
    ll = plot(scatter(x=d1, y=cdf.(d1), mode="lines"))
    d1 = collect(0:0.05:da)
    F(d) = (d / da)^nF
    PlotlyJS.add_trace!(ll, (scatter(x=d1, y=F.(d1), mode="lines")))
    PlotlyJS.add_trace!(ll, scatter(x=diameter, y=FF, mode="markers"))#u-----
    display(ll)
    return pointsresult, pointsdiameterfinal
end

global filename = "0"
global steel_uniques = []
global steel_bond_uniques = []

function Meshing(model_name, mesh_plot="Yes", mesh_dirc_and_name="LDPM_mesh_facets")

    if typeof(model_name) == ldpm
        include("1_Random_Meshing_Delaunay.jl")
        include("2_Geometrical_Properties.jl") # Evaluate Geometrical Properties of Tetrahedra
        include("3_Projected_Areas.jl") # Compute Areas of Triangles Composing Every Connection
        include("4_Connectivity_Matrix.jl") # Evaluate Connections Between Particles
        if mesh_plot == "Yes"
            include("4.5 Mesh plot.jl")
        end
    elseif typeof(model_name) == ldpmbarr
        include("1_Random_Meshing_Delaunay with_steel.jl")
        include("2_Geometrical_Properties.jl") # Evaluate Geometrical Properties of Tetrahedra
        include("3_Projected_Areas.jl") # Compute Areas of Triangles Composing Every Connection
        include("4_Connectivity_Matrix_steel.jl") # Evaluate Connections Between Particles
        if mesh_plot == "Yes"
            global filename = mesh_dirc_and_name
            include("4.5 Mesh plot.jl")
            include("4.6 Steel plot.jl")
        end
    end
    return gdl, Positions, gdl_n, Connect, TRI, Unique_Connections, Elements, Tet, Unique_Connections_Elements, B
end


"""
eigenbox(A[, method=Rohn()])

Returns an enclosure of all the eigenvalues of `A`. If `A` is symmetric, then the
output is a real interval, otherwise it is a complex interval.

### Input

- `A` -- square interval matrix
- `method` -- method used to solve the symmetric interval eigenvalue problem (bounding
    eigenvalues of general matrices is also reduced to the symmetric case).
    Possible values are

      - `Rohn` -- (default) fast method to compute an enclosure of the eigenvalues of
            a symmetric interval matrix
      - `Hertz` -- finds the exact hull of the eigenvalues of a symmetric interval
            matrix, but has exponential complexity.

### Algorithm


"""
global loaded_region = []
global plot_boundary = "0"
function Boundary_setting(loaded__region=[[[0 10; 0 200; 0 10], [1, 2, 3, 4, 5, 6], [0, 0, 0, 0, 0, 0]], [[190 200; 0 200; 0 10], [1, 2, 3], [0, 0, 0]], [[95 105; 0 200; 60 70], [3], [-0.2]]], plot__boundary="Yes")
    global loaded_region = loaded__region
    global plot_boundary = plot__boundary
    include("5_Boundary.jl")
    return Boun, Boun1
end
global Δt = 0.000001
global t_final = 0.2
function Solutions(model_name, scale_delata_time=1.0, t_final_=0.8) #Δt= round(2/median(ω_n),digits=5)
    global Δt = Δt * scale_delata_time
    global t_final = t_final_
    if typeof(model_name) == ldpm
        include("6_Material_Properties.jl") # Calculate Stable Time Step
        include("11_Boundary_Conditions.jl") # Function to Calculate Boundary Velocities
        include("12_Volumetric_Deformation.jl") # Function to Update Volumetric Deformation of Lattice Elements
        include("13_Epsilon_Update.jl") # Function to Update Maximum Deformation (Used in Constitutive Law)
        include("14_Element_Response.jl") # Function Used to Obtain Element Forces
        include("15_Constitutive_Law.jl") # Constitutive Law
        #include("15_Constitutive_Law steel and bond.jl") # Constitutive Law
        include("Central_Difference_Integration.jl")
    elseif typeof(model_name) == ldpmbarr
        include("6_Material_Properties with steel.jl") # Calculate Stable Time Step
        include("11_Boundary_Conditions.jl") # Function to Calculate Boundary Velocities
        include("12_Volumetric_Deformation.jl") # Function to Update Volumetric Deformation of Lattice Elements
        include("13_Epsilon_Update.jl") # Function to Update Maximum Deformation (Used in Constitutive Law)
        include("14_Element_Response.jl") # Function Used to Obtain Element Forces
        include("15_Constitutive_Law.jl") # Constitutive Law
        include("15_Constitutive_Law steel and bond.jl") # Constitutive Law
        include("Central_Difference_Integration_with steel.jl")
    end
    return displace, eps_, eps_n, eps_t, sigma_eff, sigma_N, sigma_T, internal, external
end
global relative_time_of_cracking = []
global crack_plot_dirc_and_name = "0"
global output_displacement_directions = []
global output_load_directions = []
global step_interval = 1
global load_dis_out_name = "0"
global plot_dis_load_region = "Yes"
function post_process(model_name, relative_time_of_cracking_=[0.4, 0.8, 1.0], crack_plot_dirc_and_name_="D:/cracking pattern", output_displacement_directions_=[[[90 110; 0 200; 0 10], [3]]], output_load_directions_=[[[90 110; 0 200; 60 70], [3]]], step_interval_=300, load_dis_out_name_="200*200*70 deck", plot_dis_load_region_="Yes")
    global relative_time_of_cracking = relative_time_of_cracking_
    global crack_plot_dirc_and_name = crack_plot_dirc_and_name_
    global output_displacement_directions = output_displacement_directions_
    global output_load_directions = output_load_directions_
    global step_interval = step_interval_
    global load_dis_out_name = load_dis_out_name_
    global plot_dis_load_region = plot_dis_load_region_
    if typeof(model_name) == ldpm
        include("vtk_cracks_ non_projected.jl") # Calculate Stable Time Step

    elseif typeof(model_name) == ldpmbarr
        gdl = gdlS
        include("vtk_cracks_ non_projected steel.jl") # Calculate Stable Time Step
    end
    include("out_load_displacement.jl")
end

export Particle_distribution, Meshing, Boundary_setting, Solutions, post_process,
    LDPM, LDPM_bar_reforced
end

#!use of package, units: mm, N, Mpa
# imen1, dimen2, dimen3, va, da, d0, nF, magnifyp
# LDPM.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1]
# dimen1, dimen2, dimen3, va, da, d0, nF, magnifyp, height_S, diameter_S
# LDPM_bar_reforced.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1, 20.0, 16.0]
# #traverse_distribution 
# LDPM_bar_reforced.steel_layout = [50, 100, 150]

# Particle_distribution(LDPM, "Yes", "D:/LDPM_geometry")
# @load "D:/LDPM_geometry.jld2"
# Meshing(LDPM, "Yes", "D:/LDPM_mesh_facets")
# Boundary_setting([[[0 10; 0 200; 0 10], [1, 2, 3, 4, 5, 6], [0, 0, 0, 0, 0, 0]], [[190 200; 0 200; 0 10], [1, 2, 3], [0, 0, 0]], [[95 105; 0 200; 60 70], [3], [-0.2]]], "Yes")

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
# #epsi_sh = 0.02
# #Esh = 833.33 #Mpa
# LDPM_bar_reforced.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 0.02, 833.33]

# Solutions(LDPM, 0.2, 0.8) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
# post_process(LDPM, [0.4, 0.8, 1.0], "D:/cracking pattern", [[[90 110; 0 200; 0 10], [3]]], [[[90 110; 0 200; 60 70], [3]]], 300, "200_200_70 deck", "Yes")
