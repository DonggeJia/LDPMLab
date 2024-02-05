# LDPMLab

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/dev/)
[![Coverage](https://codecov.io/gh/DonggeJia/LDPMLab.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DonggeJia/LDPMLab.jl)

```
#!use of package, units: mm, N, Mpa
#dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp
LDPM.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1]
#dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp, height_S, diameter_S
LDPM_bar_reforced.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1, 20.0, 16.0]
#traverse_distribution 
LDPM_bar_reforced.steel_layout = [50, 100, 150]

particle_distribution(LDPM_bar_reforced, "Yes", "D:/LDPM_geometry")
@load "D:/LDPM_geometry.jld2"
Meshing(LDPM_bar_reforced, "Yes", "D:/LDPM_mesh_facets")
Boundary_setting([[[0 10; 0 200; 0 10], [1, 2, 3, 4, 5, 6], [0, 0, 0, 0, 0, 0]], [[190 200; 0 200; 0 10], [1, 2, 3], [0, 0, 0]], [[95 105; 0 200; 60 70], [3], [-0.2]]], "Yes")

45000.0 # E_m -> Initial Elastic Modulus for the Matrix [MPa]
45000.0 # E_a -> Initial Elastic Modulus for the Aggregates [MPa]
3.0 # sigma_t -> Tension Limit [MPa]
-50.0 # sigma_c -> Compression Limit [MPa]
10.0 # sigma_s -> Shear Limit [MPa]
0.07 # G_t -> Mode I Fracture Energy [N/mm]
0.35 # G_s -> Mode II Fracture Energy [N/mm]
0.25 # alpha -> Tangential-to-Normal Stiffness Ratio (Controls Poisson's Effect) [-]
2.0 # n_t -> Exponent that Controls Transition of Softening Parameter [-]
0.8 # n_c -> Exponent that Controls Transition of Hardening Parameter [-]
2.5e-6 # rho -> Material's Mass Density [ton/mm^3]
1.0 # kc1 -> Volumetric Compression Parameter [-]
5.0 # kc2 -> Volumetric Compression Parameter [-]
11250.0 # K_c -> Parameter that Governs the Post-Peak Behavior in Compression [-] 
0.0 # zeta -> damping coefficient
LDPM.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0]
#add steel mechanical parameters
#E_steel = 1.96*10.0^5
#fyi = 500
#fu = 650 #Mpa
#epsi_sh = 0.02
#Esh = 833.33 #Mpa
LDPM_bar_reforced.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]

Solutions(LDPM_bar_reforced, 0.2, 0.2) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
post_process(LDPM_bar_reforced, [0.4, 0.8, 1.0], "D:/cracking pattern", [[[90 110; 0 200; 0 10], [3]]], [[[90 110; 0 200; 60 70], [3]]], 300, "200_200_70 deck", "Yes")
