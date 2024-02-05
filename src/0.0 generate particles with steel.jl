
height_S = LDPM_bar_reforced.geometry_parameters[end-1] #mm
traverse_distribution = LDPM_bar_reforced.steel_layout
diameter_S = LDPM_bar_reforced.geometry_parameters[end] #mm
dimen1 = LDPM_bar_reforced.geometry_parameters[1] #mm
dimen2 = LDPM_bar_reforced.geometry_parameters[2] #mm
dimen3 = LDPM_bar_reforced.geometry_parameters[3] #mm
V = dimen1 * dimen2 * dimen3
#a_over_c = 10.7 # aggregate-to-cement ratio
da = LDPM_bar_reforced.geometry_parameters[5] #mm maximum aggregate size
d0 = LDPM_bar_reforced.geometry_parameters[6] #mm mminimum aggregate size
#q = 2.5 #material parameter
nF = LDPM_bar_reforced.geometry_parameters[7]
q = 3 - nF
va = LDPM_bar_reforced.geometry_parameters[4] #aggregate volume fraction
function particle_generation()

    va0 = (1 - (d0 / da)^nF) * va#simulated aggregate volume fraction
    Va0 = va0 * V
    V_a0 = 0.0
    particledia = Float64[]
    di = d0 * (1 - rand() * (1 - d0^q / da^q))^(-1 / q)
    V_a0 = V_a0 + pi * di^3 / 6
    while V_a0 <= Va0
        push!(particledia, di)
        di = d0 * (1 - rand() * (1 - d0^q / da^q))^(-1 / q)
        V_a0 = V_a0 + pi * di^3 / 6
    end

    diameter = d0:0.1:da
    FF = Array{Float64}(undef, length(diameter))
    for j = 1:length(diameter)
        Vaggr = 0
        for i = 1:length(particledia)
            if particledia[i] < diameter[j]
                Vaggr = (Vaggr + pi * particledia[i]^3 / 6)
            end
        end
        FF[j] = (Vaggr + (d0 / da)^nF * va * V) / va / V #we know the aggregates less than the simulated minimum aggregate accounts for (d0/da)^nF, the corresponding volome is (d0/da)^nF*va*V
    end
    # PlotlyJS.add_trace!(ll, scatter(x=diameter, y=FF, mode="markers"))#u-----
    # display(ll) 
    #u----savefig(ll, "particle generation.svg")#u------
    return particledia, FF
end
diameter = d0:0.1:da
#FF = Array{Float64}(undef, length(diameter))
particledia, FF = particle_generation()