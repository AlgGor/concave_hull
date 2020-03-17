function [slope_coeff, phi_left, phi_right] = key_phi_line(x_key_left, x_key_right, x_k_1, u_k_1, x_k, u_k)

       slope_coeff = (u_k- u_k_1)/(x_k-x_k_1); 
       phi_left = u_k - slope_coeff * (x_k - x_key_left);
       phi_right = u_k - slope_coeff * (x_k - x_key_right);

end