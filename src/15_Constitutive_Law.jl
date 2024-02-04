function Constitutive_Law( eps_n,eps_t,eps_,E,sigma_eff,sigma_N,sigma_T,eps_pos_max,eps_neg_max,eps_n_max,eps_t_max,K_t,K_s,eps_vol )
#FASCETTI_CONSTITUTIVE_LAW: takes the effective and coupling strain to compute the effective stress
# parameters = [E_m E_a sigma_t sigma_c sigma_s alpha K_c k1 k2 kc2 kc1 n_c n_t]; # Store All Constant Parameters into One Vector
#
eps_[2] = sqrt(eps_n[2]^2+alpha*eps_t[2]^2); # Compute Effective Strain
omega = atan(eps_n[2]/(sqrt(alpha)*eps_t[2])); # Coupling Strain Omega
s = sin(omega); # Sine of Omega
c = cos(omega); # Cosine of Omega
# Unstressed Behavior #
if eps_[2] == 0;
    sigma_eff[2] = 0;
    sigma_N[2] = 0;
    sigma_T[2] = 0;
    return sigma_eff,sigma_N,sigma_T,eps_
end
# Tensile (Cohesive) Behavior #
 if eps_n[2] >= 0;
     if eps_[2] >= eps_pos_max # Virgin Loading Case
        sigma_zero = (s*(sigma_c+sigma_t)+sqrt(s^2*(sigma_c+sigma_t)^2+4*k2*(s^2+k1*alpha*c^2)))/(2*(s^2+c^2*alpha*k1)); # Effective Stress Boundary
        eps_0 = sigma_zero/E; # Effective Strain at Peak Stress
        eps_1 = sqrt(eps_n_max^2+alpha*eps_t_max^2); # Maximum Effective Strain (Irreversible Damage)
        H = K_s+(K_t-K_s)*(2*omega/pi)^n_t; # Softening Exponential Coefficient
        sigma_b = sigma_zero*exp(-H/sigma_zero*max(eps_1-eps_0,0)); # Stress Boundary for Current Omega
        sigma_eff[2] = min(E*eps_[2],sigma_b); # If sigma exceeds the boundary, put it equal to sigma_b
        sigma_N[2] = sigma_eff[2]*s; # Normal Stress
        sigma_T[2] = sigma_eff[2]*c*sqrt(alpha); # Tangential Stress
    elseif eps_[2] < eps_pos_max # Unloading Case
             sigma_eff[2] = max(sigma_eff[1]+E*(eps_[2]-eps_[1]),0); # Effective Stress
             sigma_zero = (s*(sigma_c+sigma_t)+sqrt(s^2*(sigma_c+sigma_t)^2+4*k2*(s^2+k1*alpha*c^2)))/(2*(s^2+c^2*alpha*k1)); # Effective Stress Boundary
             eps_0 = sigma_zero/E; # Effective Strain at Peak Stress
             eps_1 = sqrt(eps_n_max^2+alpha*eps_t_max^2); # Maximum Effective Strain (Irreversible Damage)
             H = K_s+(K_t-K_s)*(2*omega/pi)^n_t; # Softening Exponential Coefficient
             sigma_b = sigma_zero*exp(-H/sigma_zero*max(eps_1-eps_0,0)); # Stress Boundary for Current Omega
             sigma_eff[2] = min(sigma_eff[2],sigma_b); # If the calculated sigma exceeds the stress limit for the given strain history, set it equal to sigma_b
             sigma_N[2] = sigma_eff[2]*s; # Normal Stress
             sigma_T[2] = sigma_eff[2]*c*sqrt(alpha); # Tangential Stress
    end
# Compressive (Frictional) Behavior
elseif eps_n[2] < 0
     if eps_n[2] <= eps_neg_max # Virgin Loading Case
         eps_d = eps_n[2]-eps_vol;
         rdv = eps_d/eps_vol;
         H = K_c/(1+kc2*max(rdv-kc1,0)); # Exponential Parameter for the Compressive Curve
         if eps_vol >= 0;
             sigma_b = sigma_c; # For Positive Volume Increments The Compression Boundary is Set Equal to the MesoScale Compressive Resistance of the Material
         else
             sigma_b = sigma_c*exp(H/sigma_c*(eps_n[2]-sigma_c/E)); # Exponential Softening Behavior Takes Into Account Post-Peak Rehardening due to Pore Collapse
         end
         sigma_N[2] = max(E*eps_n[2],sigma_b);
     else
         sigma_N[2] = min(sigma_N[1]+E*(eps_n[2]-eps_n[1]),0);
         eps_d = eps_n[2]-eps_vol;
         rdv = eps_d/eps_vol;
         H = K_c/(1+kc2*max(rdv-kc1,0)); # Exponential Parameter for the Compressive Curve
         if eps_vol >= 0;
             sigma_b = sigma_c; # For Positive Volume Increments The Compression Boundary is Set Equal to the MesoScale Compressive Resistance of the Material
         else
             sigma_b = sigma_c*exp(H/sigma_c*(eps_n[2]-sigma_c/E)); # Exponential Softening Behavior Takes Into Account Post-Peak Rehardening due to Pore Collapse
         end
         sigma_N[2] = max(sigma_N[2],sigma_b); # Check if the Unloading-Reloading Stress Value Exceeds the Stress Boundary
     end
     #
     if eps_t[2] >= eps_t_max # Virgin Load Case
        H = K_s-K_s*(-2*omega/pi)^n_c;
        sigma_zero = (s*(sigma_c+sigma_t)+sqrt(s^2*(sigma_c+sigma_t)^2+4*k2*(s^2+k1*alpha*c^2)))/(2*(s^2+c^2*alpha*k1));
        tau_0 = sigma_zero*c*sqrt(alpha);
        sigma_T_b = tau_0*exp(-H/tau_0*max(eps_t[2]-tau_0/(alpha*E),0));
        sigma_T[2] = min(alpha*E*eps_t[2],sigma_T_b);
     else
         sigma_T[2] = max(sigma_T[1]+alpha*E*(eps_t[2]-eps_t[1]),0);
         H = K_s-K_s*(-2*omega/pi)^n_c;
         sigma_zero = (s*(sigma_c+sigma_t)+sqrt(s^2*(sigma_c+sigma_t)^2+4*k2*(s^2+k1*alpha*c^2)))/(2*(s^2+c^2*alpha*k1));
         tau_0 = sigma_zero*c*sqrt(alpha);
         sigma_T_b = tau_0*exp(-H/tau_0*max(eps_t[2]-tau_0/(alpha*E),0));
         sigma_T[2] = min(sigma_T[2],sigma_T_b);
     end
     sigma_eff[2] = sqrt(sigma_N[2]^2+(sigma_T[2]^2)/alpha); # Calculate The Effective Stress Based on the Compressive and Shear Stress Evaluated
end
return sigma_eff,sigma_N,sigma_T,eps_
end
