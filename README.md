# LDPMLab

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DonggeJia.github.io/LDPMLab.jl/dev/)
[![Coverage](https://codecov.io/gh/DonggeJia/LDPMLab.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DonggeJia/LDPMLab.jl)

- LDPMLab, is a integral, open-source and high-performance Julia package for the state of the art of Lattice discrete particle model (LDPM), first inovated by Prof. Cusatis of Northwestern in 2011. 

- This software tracks the cutting-edge advance of LDPM, incorparating its variants in simulating mechanical and frature behaviors of particle-reinforced materials with or without bar reinforcements, and mass transportation in heterogenious materials. 

- The application can also be extended to homogeneous material when particle is no longer realistically existed, becoming a computationally more efficient sustitute for FEM.

- This package is developed under the supervision of Prof. Fascetti and Prof. Brigham in University of Pittsburgh. 

# Universal applications
- Specially for:
    - Realistic meso-scale simulation of mechanical failure of particle-reinforced materials, such as concrete, shale, masonry, cementicious composites, granular rocks, polymers, etc.
    - Meso-scale mass transport simulation in damaged or non-damaged particle-reinforced materials
    - Interective multiphysics simulation of particle-reinforced materials with reinforcing bars
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
### 1. Speficiy geometrical parameters for a LDPM model
- for particle reinforced material:
    - Enter the values of material size in x, y, z dimensions (mm), cement content (kg/mm<sup>3</sup>), aggregate volume fraction(-), maximum particle size (mm), minimum particle size (mm), material parameter (-) for particle distribution that 0.5 corresponds to the classical Fuller curve, and scaling factor for minimum distance check in particle distribution that 1.0 means particle centroids must have a distance larger than 1.0*(radius of the first particle + radius of the second particle)
    ```
    LDPM.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1]
    ```
- for particle reinforced material with reforcing bars
    - Enter the values of just mentioned geometrical parameters for particle reinforced material, and then the height (in z direction) of reinforcing bars (mm) and the diameter of reinforcing bars (mm) 
    ```
    LDPM_bar_reforced.geometry_parameters = [200, 200, 70, 0.734, 15, 10, 0.45, 1.1, 20.0, 16.0]
    ```
    - Enter the transverse (in y direction) distribution of reinforcing bars (mm)
    ```
    LDPM_bar_reforced.steel_layout = [50, 100, 150]
    ```
    By default, the orientation of reinforcing bars is along the x-coordinate
### 2. Distribute particles in the material volume
- for particle reinforced material
    ```
    particle_distribution(LDPM, "Yes", "D:/LDPM_geometry")
    ```
    where "Yes" means saving a JLD2 file with particle coordinates and their corresponding diameters. By default, `Yes` is used. Other strings except `Yes` means not storing the current particle distribution. Saving particle distribution is encouraged so that it can be used again whenever Julia session is restarted, excluding the influence of random particle distribution on solutions.

    "D:/LDPM_geometry" indicates the filefolder and filename for particle distribution storage.
- for particle reinforced material with reforcing bars
    ```
    particle_distribution(LDPM_bar_reforced, "Yes", "D:/LDPM_geometry")
    ```
If particle distribution is stored, you can use
```
@load "D:/LDPM_geometry.jld2"
```
to reload all the procedures you've done in steps 1 and 2.
### 3. Meshing
This step implements the meshing process using Delauney tetrahedralization and modified Voronoi tesselation.
A demonstration of the procedure from particel distribution to Delauney tetrahedralization and Voronoi meshing is
<p align="center">
    <img src="docs/src/from particel to mesh.png" width="450"/>
</p>
- for particle reinforced material
    ```
    Meshing(LDPM, "Yes", "D:/LDPM_mesh_facets")
    ```
    where "Yes" means save a vtk file for Paraview to plot the contact facets in the material volume. 
    
    "D:/LDPM_mesh_facets" indicates the filefolder and filename for storage.
- for particle reinforced material with reforcing bars
    ```
    Meshing(LDPM_bar_reforced, "Yes", "D:/LDPM_mesh_facets")
    ```
The meshing plot for particle reinforced material with reforcing bars is like
<p align="center">
    <img src="docs/src/particel reinforced material with reforcing bar.png" width="450"/>
</p>

### 4. Set a boundary condition
    
Setting boundaries requires runing a function:
`Boundary_setting(loaded__region, plot__boundary)

`loaded_region` is a varible receiving the information of displacement-controlled boundaries. The expected input is a three-layer nested vector. The first layer is `[the first boundary, the second boundary, ..., the last boundary]`. Each boundary is a vector `[Boundary condition region, freedom degrees that boundary applied on boundary nodes, displacement velocities on these boundary freedom degrees]`. `Boundary condition region` is a $3 \times 2$ matrix where the first row indicates the boundary region in x-coordinate that `Boundary condition region[1,1]` is the starting place (mm) and 

`Boundary condition region[1,2]` is the ending place (mm), the second row indicates the boundary region in y-coordinate, and the third row indicates the boundary region in z-coordinate. The unit of `displacement velocities on these boundary freedom degrees` is mm/second.

For example,
```
Boundary_setting([[[0 10; 0 200; 0 10], [1, 2, 3, 4, 5, 6], [0, 0, 0, 0, 0, 0]], [[190 200; 0 200; 0 10], [1, 2, 3], [0, 0, 0]], [[95 105; 0 200; 60 70], [3], [-0.2]]], "Yes")
```
"Yes" means the boundary will be ploted, other values for this parameter means no plot. the plot is like
<p align="center">
    <img src="docs/src/boundary setting.png" width="450"/>
</p>
The boundary condictions are marked on each point with a format `local freedom degree_velocity`. Each point will have at most 6 marks for freedom degrees in x-direction movement, y-direction movement, z-direction movement, rotation about the x-axis, rotation about the y-axis, rotation about the z-axis.

### 5. Specify mechanical parameters for a LDPM model
- for particle reinforced material
    - Enter the values of initial elastic modulus for the matrix which wraps around the particles (MPa), initial elastic modulus for particles (MPa), tension stress limit (MPa), compression stress limit (MPa), shear stress limit (MPa), Mode I fracture energy (N/mm), Mode II fracture energy (N/mm), tangential-to-normal stiffness ratio (controls poission's efffect) (-), exponent that controls transition of softening parameter (-), exponent that controls transition of hardening parameter (-), material's mass density (ton/mm<sup>3</sup>), frist volumetric compression parameter (-), second volumetric compression parameter (-), parameter that governs the post-Peak behavior in compression (-), and damping coefficient in dynamic solution (-).
    ```
    LDPM.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0]
    ```
- for particle reinforced material with reforcing bars
    - Enter the values of the just mentioned mechanical parameters for particle reinforced material, and then the paramenters for reinforcing bar's constitutive law: initial elastic modulus (MPa), yield plateau strength (MPa), hardening strain (-), and the slope of hardening (MPa) in sequence.
    ```
    LDPM_bar_reforced.mechanical_parameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 0.02, 833.33]
    ```
    Here the trilinear simplified constitutive model for steel materials is adopted.

### 6. Solution
solving the model is realized by
```
Solutions(LDPM, 0.2, 0.8)
```
for particle reinforced materials and
```
Solutions(LDPM_bar_reforced, 0.2, 0.8)
```
for particle reinforced materials with reforcing bars.

`0.2` here gives the value of scaling factor for time step. The time step is automatically determined by computing the natural frequencies of LDPM cells. `0.2` means one fifth of the calculated time step is used for better stability (-).

`1.1` indicates the value of `total time` the model will run (second). 

After this step, you already get the solution. A further step helps to output results and plot beautiful cracking pattern.
### 7. Post process
    
Post process uses a function `post_process(model_name, relative_time_of_cracking_=[0.4, 0.8, 1.0], crack_plot_dirc_and_name_="D:/cracking pattern", output_displacement_directions_=[[[90 110; 0 200; 0 10], [3]]], output_load_directions_=[[[90 110; 0 200; 60 70], [3]]], step_interval_=300, load_dis_out_name_="D:/200_200_70_deck", plot_dis_load_region_="Yes")`.

`model_name` should be `LDPM` for particle reinforced materials and `LDPM_bar_reforced` for particle reinforced materials with reforcing bars.

This function enables you output cracking data as vtk files at different time steps during the solution by `relative_time_of_cracking_`. By default, the cracking pattern plots will be generated at the time steps around `0.4` $\times$ `total time` , `0.8` $\times$ `total time`, and `1.0` $\times$ `total time`. 

`crack_plot_dirc_and_name_` indicates the filefolder and filename for the cacking plot storage.

`output_displacement_directions_` is a three-layer nested vector. The first layer is `[the first region for displacement output, the second region for displacement output, ..., the last region for displacement output]`. Each region is a vector `[regional range, freedom degrees that output]`. `regional range` is a $3 \times 2$ matrix where the first row indicates the range of the region in x-coordinate that `regional range[1,1]` is the starting place (mm) and `regional range[1,2]` is the ending place (mm), the second row indicates the range of the region in y-coordinate, and the third row indicates the range of the region in z-coordinate. 

`output_load_directions_` has the same data structure as `output_displacement_directions_`, but the output is residual load values.
If `plot_dis_load_region` is given a "Yes", the regions for displacement and load output will be plotted like
<p align="center">
    <img src="docs/src/data ouput regions.png" width="450"/>
</p>
The points that outputed data comes from have a 6-number string representing 6 degrees of freedom in sequence. `1` means the displacement data is extracted, `2` means the residual load data is extracted, and `3` means both displacement and residual load data are drawn.

`step_interval_` claims the value of step interval that the output data will be thined out so that the displacement and residual load will be extracted every `step_interval_` time step. 
Notice that for the stability of dynamical solution, the first 500 time steps in the package are shrinked, while the following time steps use same step length.

`load_dis_out_name_` indicates the filefolder and filename for the displacement and load data storage.

The cracking pattern for particle reinforced materials in Paraview is like
<p align="center">
    <img src="docs/src/data ouput regions.png" width="450"/>
</p>
## References
Excellent introductions of LDPM are:

Cusatis, G., Pelessone, D., & Mencarelli, A. (2011). Lattice discrete particle model (LDPM) for failure behavior of concrete. I: Theory. Cement and Concrete Composites, 33(9), 881-890.

Fascetti, A., Bolander, J. E., & Nisticó, N. (2018). Lattice discrete particle modeling of concrete under compressive loading: Multiscale experimental approach for parameter determination. Journal of Engineering Mechanics, 144(8), 04018058.
