var documenterSearchIndex = {"docs":
[{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Pages = [\"Tutorials.md\"]","category":"page"},{"location":"Tutorials/#LDPMLab","page":"Tutorials","title":"LDPMLab","text":"","category":"section"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Documentation for LDPMLab.","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"#!use of package, units: mm, N, Mpa #dimen1, dimen2, dimen3, c, woverc, da, d0, nF, Rhoc, Rhow, vair, magnifyp LDPM.geometryparameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1] #dimen1, dimen2, dimen3, c, woverc, da, d0, nF, Rhoc, Rhow, vair, magnifyp, heightS, diameterS LDPMwsteel.geometryparameters = [200, 200, 70, 190 * 10^-9, 0.9, 15, 10, 0.45, 3150 * 10^-9, 1000 * 10^-9, 3.5 / 100, 1.1, 20.0, 16.0] #traversedistribution  LDPMwsteel.steellayout = [50, 100, 150]","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"particledistribution(; modelname, Geometrysave=\"Yes\", Geometrydircandname=\"LDPMgeometry\") loadparticledistribution(particledircandname=\"LDPMparticledistribution\") Meshing(modelname, meshplot=\"Yes\", meshdircandname=\"LDPMmeshfacets\") Boundarysetting(fixedregion=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], loadedregion=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], plot_boundary=\"Yes\")","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"45000.0 # Em -> Initial Elastic Modulus for the Matrix [MPa] 45000.0 # Ea -> Initial Elastic Modulus for the Aggregates [MPa] 3.0 # sigmat -> Tension Limit [MPa] -50.0 # sigmac -> Compression Limit [MPa] 10.0 # sigmas -> Shear Limit [MPa] 0.07 # Gt -> Mode I Fracture Energy [N/mm] 0.35 # Gs -> Mode II Fracture Energy [N/mm] 0.25 # alpha -> Tangential-to-Normal Stiffness Ratio (Controls Poisson's Effect) [-] 2.0 # nt -> Exponent that Controls Transition of Softening Parameter [-] 0.8 # nc -> Exponent that Controls Transition of Hardening Parameter [-] 2.5e-6 # rho -> Material's Mass Density [ton/mm^3] 1.0 # kc1 -> Volumetric Compression Parameter [-] 5.0 # kc2 -> Volumetric Compression Parameter [-] 11250.0 # Kc -> Parameter that Governs the Post-Peak Behavior in Compression [-]  0.0 # zeta -> damping coefficient LDPM.mechanicalparameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0] #add steel mechanical parameters #Esteel = 1.96*10.0^5 #fyi = 500 #fu = 650 #Mpa #epsish = 0.02 #Esh = 833.33 #Mpa LDPMwsteel.mechanicalparameters = [45000.0, 45000.0, 3.0, -50.0, 10.0, 0.07, 0.35, 0.25, 2.0, 0.8, 2.5e-6, 1.0, 5.0, 11250.0, 0.0, 1.96 * 10^5, 500, 650, 0.02, 833.33]","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Solutions(modelname, scaledelatatime=1.0, tfinal=0.05) #loading [velocity, direction], Δt= round(2/median(ωn),digits=5) postprocess(modelname, relativetimeofcracking=[0.2, 0.4, 0.5], crackplotdircandname=\"cracking pattern\", outputdisplacementdirections=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], outputloaddirections=[[[0 10; 0 200; 65 70], [3, 4, 5]], [[100 110; 0 200; 65 70], [3, 4, 5]], [[190 200; 0 200; 65 70], [3, 4, 5]]], stepinterval=300, loaddisoutname=\"20020070 deck\", plotdisload_region=\"Yes\")","category":"page"},{"location":"Tutorials/","page":"Tutorials","title":"Tutorials","text":"Modules = [LDPMLab]","category":"page"},{"location":"references/","page":"References","title":"References","text":"CurrentModule = LDPMLab","category":"page"},{"location":"references/#LDPMLab","page":"References","title":"LDPMLab","text":"","category":"section"},{"location":"references/","page":"References","title":"References","text":"Documentation for LDPMLab.","category":"page"},{"location":"references/","page":"References","title":"References","text":"","category":"page"},{"location":"references/","page":"References","title":"References","text":"Modules = [LDPMLab]","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"CurrentModule = LDPMLab","category":"page"},{"location":"mechanical_response/#LDPMLab","page":"Mechanical response","title":"LDPMLab","text":"","category":"section"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"Documentation for LDPMLab.","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"","category":"page"},{"location":"mechanical_response/","page":"Mechanical response","title":"Mechanical response","text":"Modules = [LDPMLab]","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"CurrentModule = LDPMLab","category":"page"},{"location":"contributing_guide/#LDPMLab","page":"Contributing","title":"LDPMLab","text":"","category":"section"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"Documentation for LDPMLab.","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"","category":"page"},{"location":"contributing_guide/","page":"Contributing","title":"Contributing","text":"Modules = [LDPMLab]","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"CurrentModule = LDPMLab","category":"page"},{"location":"mass_transport/#LDPMLab","page":"Mass transport","title":"LDPMLab","text":"","category":"section"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"Documentation for LDPMLab.","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"","category":"page"},{"location":"mass_transport/","page":"Mass transport","title":"Mass transport","text":"Modules = [LDPMLab]","category":"page"},{"location":"Embeding_linear_Reinforcement/","page":"Embeding linear Reinforcement","title":"Embeding linear Reinforcement","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Embeding_linear_Reinforcement/#LDPMLab","page":"Embeding linear Reinforcement","title":"LDPMLab","text":"","category":"section"},{"location":"Embeding_linear_Reinforcement/","page":"Embeding linear Reinforcement","title":"Embeding linear Reinforcement","text":"Documentation for LDPMLab.","category":"page"},{"location":"Embeding_linear_Reinforcement/","page":"Embeding linear Reinforcement","title":"Embeding linear Reinforcement","text":"","category":"page"},{"location":"Embeding_linear_Reinforcement/","page":"Embeding linear Reinforcement","title":"Embeding linear Reinforcement","text":"Modules = [LDPMLab]","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = LDPMLab","category":"page"},{"location":"#LDPMLab","page":"Home","title":"LDPMLab","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LDPMLab.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LDPMLab]","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"CurrentModule = LDPMLab","category":"page"},{"location":"Real applications/#LDPMLab","page":"Real applications","title":"LDPMLab","text":"","category":"section"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"Documentation for LDPMLab.","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"","category":"page"},{"location":"Real applications/","page":"Real applications","title":"Real applications","text":"Modules = [LDPMLab]","category":"page"},{"location":"API/","page":"API","title":"API","text":"CurrentModule = LDPMLab","category":"page"},{"location":"API/#LDPMLab","page":"API","title":"LDPMLab","text":"","category":"section"},{"location":"API/","page":"API","title":"API","text":"Documentation for LDPMLab.","category":"page"},{"location":"API/","page":"API","title":"API","text":"","category":"page"},{"location":"API/","page":"API","title":"API","text":"Modules = [LDPMLab]","category":"page"}]
}
