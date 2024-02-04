function eps_update!(eps_,eps_n,eps_t,eps_pos_max,eps_neg_max,eps_n_max,eps_t_max)
# EPS_UPDATE Update Maximum and Minimum Epsilon
    if eps_n > eps_pos_max # Check if the Actual Value of the Normal Strain is Greater Than the Maximum Previously Attained
        eps_pos_max = eps_; # (Effective Strain Value!)
        eps_n_max = eps_n; # (Normal Strain Value!)
    elseif eps_n < eps_neg_max # Check if the Actual Value of the Normal Strain is Greater Than the Maximum Previously Attained (Normal Strain Value!)
        eps_neg_max = eps_n; # Normal Strain Update
    end
    if eps_t > eps_t_max
        eps_t_max = eps_t; # If the actual value of the shear strain exceeds the previously attained maximum value, update
    end
    return eps_pos_max, eps_n_max, eps_neg_max, eps_t_max
end
