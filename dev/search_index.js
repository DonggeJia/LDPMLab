var documenterSearchIndex = {"docs":
[{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Pages = [\"Tutorials.md\"]","category":"page"},{"location":"Tutorials/#LDPMLab","page":"Tutorials","title":"LDPMLab","text":"","category":"section"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Documentation for LDPMLab.","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"#!use of package, units: mm, N, Mpa #dimen1, dimen2, dimen3, c, woverc, da, d0, nF, Rhoc, Rhow, vair, magnifyp LDPM.geometryparameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1] #dimen1, dimen2, dimen3, c, woverc, da, d0, nF, Rhoc, Rhow, vair, magnifyp, heightS, diameterS LDPMwsteel.geometryparameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1, 20.0, 16.0] #traversedistribution  LDPMwsteel.steellayout = [50, 100, 150]","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"particledistribution(; modelname, Geometrysave=\"Yes\", Geometrydircandname=\"LDPMgeometry\") loadparticledistribution(particledircandname=\"LDPMparticledistribution\") Meshing(modelname, meshplot=\"Yes\", meshdircandname=\"LDPMmeshfacets\") Boundarysetting(fixedregion=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], loadedregion=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], plot_boundary=\"Yes\")","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"45000.0 # Em -> Initial Elastic Modulus for the Matrix [MPa] 45000.0 # Ea -> Initial Elastic Modulus for the Aggregates [MPa] 3.0 # sigmat -> Tension Limit [MPa] -50.0 # sigmac -> Compression Limit [MPa] 10.0 # sigmas -> Shear Limit [MPa] 0.07 # Gt -> Mode I Fracture Energy [N/mm] 0.35 # Gs -> Mode II Fracture Energy [N/mm] 0.25 # alpha -> Tangential-to-Normal Stiffness Ratio (Controls Poisson's Effect) [-] 2.0 # nt -> Exponent that Controls Transition of Softening Parameter [-] 0.8 # nc -> Exponent that Controls Transition of Hardening Parameter [-] 2.5e-6 # rho -> Material's Mass Density [ton/mm^3] 1.0 # kc1 -> Volumetric Compression Parameter [-] 5.0 # kc2 -> Volumetric Compression Parameter [-] 11250.0 # Kc -> Parameter that Governs the Post-Peak Behavior in Compression [-]  0.0 # zeta -> damping coefficient LDPM.mechanicalparameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0] #add steel mechanical parameters #Esteel = 1.96*10.0^5 #fyi = 500 #fu = 650 #Mpa #epsish = 0.02 #Esh = 833.33 #Mpa LDPMwsteel.mechanicalparameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Solutions(modelname, scaledelatatime=1.0, tfinal=0.05) #loading [velocity, direction], Δt= round(2/median(ωn),digits=5) postprocess(modelname, relativetimeofcracking=[0.2, 0.4, 0.5], crackplotdircandname=\"cracking pattern\", outputdisplacementdirections=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], outputloaddirections=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], stepinterval=300, loaddisoutname=\"20020070 deck\", plotdisload_region=\"Yes\")","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Hello","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = LDPMLab","category":"page"},{"location":"#LDPMLab","page":"Home","title":"LDPMLab","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LDPMLab.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Hello","category":"page"},{"location":"","page":"Home","title":"Home","text":"1+1 = A","category":"page"},{"location":"","page":"Home","title":"Home","text":"[1,1]+[0.2,0.4]","category":"page"},{"location":"","page":"Home","title":"Home","text":"Hello","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LDPMLab]","category":"page"},{"location":"#LDPMLab.Boundary_setting","page":"Home","title":"LDPMLab.Boundary_setting","text":"eigenbox(A[, method=Rohn()])\n\nReturns an enclosure of all the eigenvalues of A. If A is symmetric, then the output is a real interval, otherwise it is a complex interval.\n\nInput\n\nA – square interval matrix\nmethod – method used to solve the symmetric interval eigenvalue problem (bounding   eigenvalues of general matrices is also reduced to the symmetric case).   Possible values are\n- `Rohn` -- (default) fast method to compute an enclosure of the eigenvalues of\n      a symmetric interval matrix\n- `Hertz` -- finds the exact hull of the eigenvalues of a symmetric interval\n      matrix, but has exponential complexity.\n\nAlgorithm\n\n\n\n\n\n","category":"function"},{"location":"#LDPMLab.Particle_distribution","page":"Home","title":"LDPMLab.Particle_distribution","text":"particle_distribution(model_name, particle_save=\"Yes\", particle_dirc_and_name=\"LDPM_particle_distribution\")\n\nExcute particle distribution in material volume and save the data of particle distribution.\n\nInput\n\nmodel_name – symbol describing the model type used\n:LDPM – uses LDPM\n:LDPM_bar_reforced – uses rank1 update\nparticle_save – symbol describing wether saving particle distribution\n:Yes – save a JLD2 file with particle coordinates and their corresponding diameters. By default, :Yes is used. Other strings except Yes means not storing the current particle distribution.\nSaving particle distribution is encouraged so that the following solutions can use the same model mesh, excluding the influence of meshing on solutions.\nparticle_dirc_and_name – a string indicating the filefolder and filename saving particle distribution. By default, LDPM_particle_distribution  is used. This argument can be ../output examples/Juliafiles/LDPM_particle_distribution as well.\n\n\n\n\n\n","category":"function"},{"location":"references/","page":"References","title":"References","text":"CurrentModule = LDPMLab","category":"page"},{"location":"references/#LDPMLab","page":"References","title":"LDPMLab","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Documentation for LDPMLab.","category":"page"},{"location":"references/","page":"References","title":"References","text":"","category":"page"},{"location":"references/","page":"References","title":"References","text":"Hello","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"CurrentModule = LDPMLab","category":"page"},{"location":"mechanical_response/#LDPMLab","page":"Mechanical response","title":"LDPMLab","text":"","category":"section"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"Documentation for LDPMLab.","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"Hello","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"CurrentModule = LDPMLab","category":"page"},{"location":"contributing_guide/#LDPMLab","page":"Contributing","title":"LDPMLab","text":"","category":"section"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"Documentation for LDPMLab.","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"<!– @autodocs Modules = [LDPMLab] –> Hi","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"CurrentModule = LDPMLab","category":"page"},{"location":"mass_transport/#LDPMLab","page":"Mass transport","title":"LDPMLab","text":"","category":"section"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"Documentation for LDPMLab.","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"Hello","category":"page"},{"location":"Embedded_bar_Reinforcement/","page":"Bar-Reinforced LDPM","title":"Bar-Reinforced LDPM","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Embedded_bar_Reinforcement/#LDPMLab","page":"Bar-Reinforced LDPM","title":"LDPMLab","text":"","category":"section"},{"location":"Embedded_bar_Reinforcement/","page":"Bar-Reinforced LDPM","title":"Bar-Reinforced LDPM","text":"Documentation for LDPMLab.","category":"page"},{"location":"Embedded_bar_Reinforcement/","page":"Bar-Reinforced LDPM","title":"Bar-Reinforced LDPM","text":"","category":"page"},{"location":"Embedded_bar_Reinforcement/","page":"Bar-Reinforced LDPM","title":"Bar-Reinforced LDPM","text":"Hello","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Real applications/#LDPMLab","page":"Real applications","title":"LDPMLab","text":"","category":"section"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"Documentation for LDPMLab.","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"Hello","category":"page"},{"location":"API/","page":"API","title":"API","text":"CurrentModule = LDPMLab","category":"page"},{"location":"API/#LDPMLab","page":"API","title":"LDPMLab","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Documentation for LDPMLab.","category":"page"},{"location":"API/","page":"API","title":"API","text":"","category":"page"},{"location":"API/","page":"API","title":"API","text":"Hello","category":"page"}]
}