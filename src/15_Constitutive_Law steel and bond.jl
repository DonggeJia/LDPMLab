const E_steel = LDPM_w_steel.mechanical_parameters[end-3]
const fyi = LDPM_w_steel.mechanical_parameters[end-2]
#fu = 650 #Mpa
const epsi_sh = LDPM_w_steel.mechanical_parameters[end-1]
const Esh = LDPM_w_steel.mechanical_parameters[end] #Mpa
function Constitutive_Law_steel(eps__n)
    if eps__n*E_steel < fyi
        return eps__n*E_steel/20
    else
        if eps__n<epsi_sh
            return fyi/20
        else
            return (fyi+(eps__n-epsi_sh)*Esh)/20
        end
    end
end
using PlotlyJS
function Constitutive_Law_bond(slip__t)
    sigma_T = 54.499/(1+exp(-4.6978*(slip__t+0.0891)))-2.6149*slip__t-32.86585-0.00478
    return sigma_T/20
end
# z=collect(-2:0.0001:3)
# plot(scatter(x=z, y=Constitutive_Law_bond.(z), mode="lines"))