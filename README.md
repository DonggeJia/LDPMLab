# LDPMLab

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/dev/)
[![Coverage](https://codecov.io/gh/DonggeJia/LDPMLab.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DonggeJia/LDPMLab.jl)

LDPMLab, is a integral, open-source and high-performance Julia package for the state of the art of Lattice discrete particle model (LDPM), first inovated by Prof Cusatis of Northwestern in 2011. This software tracks the cutting-edge advance of LDPM, incorparating its variants in simulating mechanical and frature behaviors of particle-reinforced materials with or without bar reinforcements, and mass transportation in heterogenious materials. The application can also be extended to homogeneous material when particle is no longer realistically existed, becoming a computationally more efficient sustitute for FEM.

This package is developed under the supervision of [Prof. Fascetti]() and Prof. Brigham in University of Pittsburgh. 

# Universal applications
- Specially for:
    - Realistic meso-scale simulation of mechanical failure of particle-reinforced materials, such as concrete, shale, masonry, cementicious composites, granular rocks, polymers, etc.
    - Meso-scale mass transport simulation in damaged or non-damaged particle-reinforced materials
    - Interective multiphysics simulation of particle-reinforced materials with bar reinforcement
    - Both static and dynamic loading conditions and transport boundaries are applicable
    - Advanced customizable nonlinear meso-scale constitutive laws for element-wise mechanical response and mass transportation
- More efficient sustitute for FEM in simulating non-soft materials: 
    - Setting particle zero radius for homogeneous material simulation
    - Customizing local meshing for heterogeous materials other than the particle-reinforced 

# Features
- Easy to use.
- Both static and dynamic solution strategies are incorporated. Generally, the dynamic solution is recommended for its guaranteed convergence. Static solution is more accurate for quasi-satitic performance, but may have difficulty in convergence for models with complicated constitutive laws or meshing 
- Strict implimentation of random particle distribution and distance check between adjacent particles
- Output files can be viewed in [Paraview](https://www.paraview.org/) and CSV readers

# Get started
## Installation
```
import Pkg
Pkg.add(url="git@github.com:DonggeJia/LDPMLab.git")
```

this will download the package and all the necessary dependencies for you. Next you can import the package with

```
using LDPMLab
```
and you are ready to go.

## Workflow for application
1. Speficiy geometry parameters for a LDPM model
- for particle reinforced material:
    - Enter the values of material size in x, y, z dimensions (mm), cement content (kg/mm^3), aggregate volume fraction(-), maximum particle size (mm), minimum particle size (mm), material parameter (-) for particle distribution that 0.5 corresponds to the classical Fuller curve, and scale factor for minimum distance check in particle distribution that 1.0 means particle centroids must have a distance larger than 1.0*(radius of the first particle + radius of the second particle)
    ```
    
    ```
- for 
    - Enter the values of material size in x, y, z dimensions (mm), cement content (kg/mm^3), aggregate volume fraction(-), maximum particle size (mm), minimum particle size (mm), material parameter (-) for particle distribution that 0.5 corresponds to the classical Fuller curve, and scale factor for minimum distance check in particle distribution that 1.0 means particle centroids must have a distance larger than 1.0*(radius of the first particle + radius of the second particle)


```
# imen1, dimen2, dimen3, va, da, d0, nF, magnifyp
# LDPM.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1]
# dimen1, dimen2, dimen3, va, da, d0, nF, magnifyp, height_S, diameter_S
# LDPM_bar_reforced.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1, 20.0, 16.0]
# #traverse_distribution 
# LDPM_bar_reforced.steel_layout = [50, 100, 150]

# particle_distribution(LDPM_bar_reforced, "Yes", "D:/LDPM_geometry")
# @load "D:/LDPM_geometry.jld2"
# Meshing(LDPM_bar_reforced, "Yes", "D:/LDPM_mesh_facets")
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
# #fu = 650 #Mpa
# #epsi_sh = 0.02
# #Esh = 833.33 #Mpa
# LDPM_bar_reforced.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]

# Solutions(LDPM_bar_reforced, 0.2, 0.2) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
# post_process(LDPM_bar_reforced, [0.4, 0.8, 1.0], "D:/cracking pattern", [[[90 110; 0 200; 0 10], [3]]], [[[90 110; 0 200; 60 70], [3]]], 300, "200_200_70 deck", "Yes")
```
