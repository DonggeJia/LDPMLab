const G_Steel = 79.3*10^3 #Mpa
function Element_Response(eps_n,eps_t,eps_,E,sigma_N,sigma_T,sigma_eff,eps_pos_max,eps_neg_max,eps_n_max,eps_t_max,K_t,K_s,eps_vol,area,len,el_disp,B)
# ELEMENT_RESPONSE Calculate element nodal forces (Calls Constitutive Law)
    eps_n[2] = B.BN2'*el_disp[7:12]-B.BN1'*el_disp[1:6]; # Normal Strain
    eps_l = B.BL2'*el_disp[7:12]-B.BL1'*el_disp[1:6]; # Strain Parallel to l-axis
    eps_m = B.BM2'*el_disp[7:12]-B.BM1'*el_disp[1:6]; # Strain Parallel to m-axis
    eps_t[2] = sqrt(eps_l^2+eps_m^2); # Tangential Strain

    ### elastic_problem
    # sigma_N[2]=eps_n[2]*E
    # sigma_T[2]=eps_t[2]*E*alpha
    # sigma_eff[2]=sqrt(sigma_N[2]^2+sigma_T[2]^2/alpha)
    # eps_[2] = sqrt(eps_n[2]^2+alpha*eps_t[2]^2);                                  # Compute Effective Strain

    
    ####nonlinear behavior, return sigma_eff, sigma_N, sigma_T, eps_ 
    sigma_eff, sigma_N, sigma_T, eps_ = Constitutive_Law(eps_n, eps_t, eps_, E, sigma_eff, sigma_N, sigma_T, eps_pos_max, eps_neg_max, eps_n_max, eps_t_max, K_t, K_s, eps_vol); # Evaluate Stresses from the Constitutive Law
    
    
    if eps_t[2] != 0 # Avoids Division by 0
        sigma_l = sigma_T[2]*eps_l/eps_t[2]; # m-component of the calculated stress (Principle of Virtual Power)
        sigma_m = sigma_T[2]*eps_m/eps_t[2]; # n-component of the calculated stress (Principle of Virtual Power)
    else
        sigma_l = 0;
        sigma_m = 0;
    end
    f1 = -area*len*(sigma_N[2]*B.BN1+sigma_l*B.BL1+sigma_m*B.BM1); # Node_I Forces in Global Coordinates
    f2 = area*len*(sigma_N[2]*B.BN2+sigma_l*B.BL2+sigma_m*B.BM2); # Node_J Forces in Global Coordinates
    qi = [f1;f2]; # Reconstruct internal forces
    return qi,eps_,eps_n,eps_t,sigma_eff,sigma_N,sigma_T
end

function Element_Response_steel(eps_n,eps_t,eps_,E,sigma_N,sigma_T,sigma_eff,eps_pos_max,eps_neg_max,eps_n_max,eps_t_max,K_t,K_s,eps_vol,area,len,el_disp,B)
    # ELEMENT_RESPONSE Calculate element nodal forces (Calls Constitutive Law)
        eps_n[2] = B.BN2'*el_disp[7:12]-B.BN1'*el_disp[1:6]; # Normal Strain
        eps_l = B.BL2'*el_disp[7:12]-B.BL1'*el_disp[1:6]; # Strain Parallel to l-axis
        eps_m = B.BM2'*el_disp[7:12]-B.BM1'*el_disp[1:6]; # Strain Parallel to m-axis
        eps_t[2] = sqrt(eps_l^2+eps_m^2); # Tangential Strain
    
        ### elastic_problem
        # sigma_N[2]=eps_n[2]*E
        # sigma_T[2]=eps_t[2]*E*alpha
        # sigma_eff[2]=sqrt(sigma_N[2]^2+sigma_T[2]^2/alpha)
        # eps_[2] = sqrt(eps_n[2]^2+alpha*eps_t[2]^2);                                  # Compute Effective Strain
    
        
        ####nonlinear behavior, return sigma_eff, sigma_N, sigma_T, eps_ 
        
        sigma_N[2] = Constitutive_Law_steel(eps_n[2])
        sigma_T[2] = eps_t[2]*2*G_Steel
        
        if eps_t[2] != 0 # Avoids Division by 0
            sigma_l = sigma_T[2]*eps_l/eps_t[2]; # m-component of the calculated stress (Principle of Virtual Power)
            sigma_m = sigma_T[2]*eps_m/eps_t[2]; # n-component of the calculated stress (Principle of Virtual Power)
        else
            sigma_l = 0;
            sigma_m = 0;
        end
        f1 = -area*len*(sigma_N[2]*B.BN1+sigma_l*B.BL1+sigma_m*B.BM1); # Node_I Forces in Global Coordinates
        f2 = area*len*(sigma_N[2]*B.BN2+sigma_l*B.BL2+sigma_m*B.BM2); # Node_J Forces in Global Coordinates
        qi = [f1;f2]; # Reconstruct internal forces
        return qi,eps_,eps_n,eps_t,sigma_eff,sigma_N,sigma_T
    end

    
    function Element_Response_bond(eps_n,eps_t,eps_,E,sigma_N,sigma_T,sigma_eff,eps_pos_max,eps_neg_max,eps_n_max,eps_t_max,K_t,K_s,eps_vol,area, len,el_disp,B)
        # ELEMENT_RESPONSE Calculate element nodal forces (Calls Constitutive Law)
            eps_n[2] = B.BN2'*el_disp[7:12]-B.BN1'*el_disp[1:6]; # Normal Strain
            eps_l = B.BL2'*el_disp[7:12]-B.BL1'*el_disp[1:6]; # Strain Parallel to l-axis
            eps_m = B.BM2'*el_disp[7:12]-B.BM1'*el_disp[1:6]; # Strain Parallel to m-axis
            eps_t[2] = sqrt(eps_l^2+eps_m^2); # Tangential Strain
        
            ### elastic_problem
            # sigma_N[2]=eps_n[2]*E
            # sigma_T[2]=eps_t[2]*E*alpha
            # sigma_eff[2]=sqrt(sigma_N[2]^2+sigma_T[2]^2/alpha)
            # eps_[2] = sqrt(eps_n[2]^2+alpha*eps_t[2]^2);                                  # Compute Effective Strain
        
            
            ####nonlinear behavior, return sigma_eff, sigma_N, sigma_T, eps_ 
            
            sigma_N[2] = eps_n[2]*2*E[1]
            sigma_T[2] = Constitutive_Law_bond(eps_t[2]*len)
            
            if eps_t[2] != 0 # Avoids Division by 0
                sigma_l = sigma_T[2]*eps_l/eps_t[2]; # m-component of the calculated stress (Principle of Virtual Power)
                sigma_m = sigma_T[2]*eps_m/eps_t[2]; # n-component of the calculated stress (Principle of Virtual Power)
            else
                sigma_l = 0;
                sigma_m = 0;
            end
            f1 = -area*len*(sigma_N[2]*B.BN1+sigma_l*B.BL1+sigma_m*B.BM1); # Node_I Forces in Global Coordinates
            f2 = area*len*(sigma_N[2]*B.BN2+sigma_l*B.BL2+sigma_m*B.BM2); # Node_J Forces in Global Coordinates
            qi = [f1;f2]; # Reconstruct internal forces
            return qi,eps_,eps_n,eps_t,sigma_eff,sigma_N,sigma_T
    end
           