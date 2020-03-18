function ind_z_k = find_z_k(ind_phi_left, ind_phi_right, x_prev, u_prev, BETA)
    
    ind_phi_right = ind_phi_right - ind_phi_left + 1;
    x_1 = x_prev(ind_phi_right); x_2 = x_prev(end);
    u_1= u_prev(ind_phi_right); u_2 = u_prev(end);
    
    x_cur_vec = x_prev((ind_phi_right+1) : end);
    phi_cur_vec = u_2 - (u_2 - u_1)./(x_2 - x_1) .* (x_2 - x_cur_vec);
    ind_w_k = find(u_prev((ind_phi_right+1) : end) - phi_cur_vec, 1);
    ind_z_k = ind_phi_left + find(x_prev - x_prev(ind_w_k)/BETA, 1) - 1;
    
end