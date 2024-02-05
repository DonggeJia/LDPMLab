
dimen1 = LDPM.geometry_parameters[1] #mm
dimen2 = LDPM.geometry_parameters[2] #mm
dimen3 = LDPM.geometry_parameters[3] #mm
V = dimen1 * dimen2 * dimen3
#a_over_c = 10.7 # aggregate-to-cement ratio
da = LDPM.geometry_parameters[5] #mm maximum aggregate size
d0 = LDPM.geometry_parameters[6] #mm mminimum aggregate size
#q = 2.5 #material parameter
nF = LDPM.geometry_parameters[7]
q = 3 - nF
# d1=collect(d0:0.05:da)
# psd(d)= q*d0^q/(1-(d0/da)^q)/d^q+1
# psd.(d1)
# plot(scatter(x=d1, y=psd.(d1), mode="lines"))
# cdf(d) = ((1-(d0/d)^q))/((1-(d0/da)^q))
# cdf.(d1)
# ll= plot(scatter(x=d1, y=cdf.(d1), mode="lines"))
# d1=collect(0:0.05:da)
# F(d) = (d/da)^nF     
# F.(d1) 
# PlotlyJS.add_trace!(ll, (scatter(x=d1, y=F.(d1), mode="lines")))   

va = LDPM.geometry_parameters[4] #aggregate volume fraction
function particle_generation()

    # Rhoc = LDPM.geometry_parameters[9] #kg/mm3
    # Rhow = LDPM.geometry_parameters[10] #kg/mm3
    # vair = LDPM.geometry_parameters[11] #typically 3%-4%
    # w = w_over_c * c
    # # d1=collect(d0:0.05:da)
    # # psd(d)= q*d0^q/(1-(d0/da)^q)/d^q+1
    # # psd.(d1)
    # # plot(scatter(x=d1, y=psd.(d1), mode="lines"))
    # # cdf(d) = ((1-(d0/d)^q))/((1-(d0/da)^q))
    # # cdf.(d1)
    # # ll= plot(scatter(x=d1, y=cdf.(d1), mode="lines"))
    # # d1=collect(0:0.05:da)
    # # F(d) = (d/da)^nF     
    # # F.(d1) 
    # # PlotlyJS.add_trace!(ll, (scatter(x=d1, y=F.(d1), mode="lines")))   

    # va = 1 - c / Rhoc - w / Rhow - vair #aggregate volume fraction
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