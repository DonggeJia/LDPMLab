# LDPMLab

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/dev/)
[![Coverage](https://codecov.io/gh/DonggeJia/LDPMLab.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DonggeJia/LDPMLab.jl)

```
#!use of package, units: mm, N, Mpa
#dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp
LDPM.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1]
#dimen1, dimen2, dimen3, c, w_over_c, da, d0, nF, Rhoc, Rhow, vair, magnifyp, height_S, diameter_S
LDPM_w_steel.geometry_parameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1, 20.0, 16.0]
#traverse_distribution 
LDPM_w_steel.steel_layout = [50, 100, 150]

particle_distribution(; model_name, Geometry_save="Yes", Geometry_dirc_and_name="LDPM_geometry")
load_particle_distribution(particle_dirc_and_name="LDPM_particle_distribution")
Meshing(model_name, mesh_plot="Yes", mesh_dirc_and_name="LDPM_mesh_facets")
Boundary_setting(fixed_region=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], loaded_region=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], plot_boundary="Yes")

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
LDPM_w_steel.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]

Solutions(model_name, scale_delata_time=1.0, t_final=0.05) #loading [velocity, direction], Δt= round(2/median(ω_n),digits=5)
post_process(model_name, relative_time_of_cracking=[0.2, 0.4, 0.5], crack_plot_dirc_and_name="cracking pattern", output_displacement_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], output_load_directions=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], step_interval=300, load_dis_out_name="200*200*70 deck", plot_dis_load_region="Yes")
```