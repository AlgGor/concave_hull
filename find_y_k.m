function ind_y_k = find_y_k(ind_left_key_pt, ind_phi_left, x_prev, u_prev, BETA)
    
    ind_phi_left = ind_phi_left - ind_left_key_pt + 1;
    x_1 = x_prev(ind_phi_left); x_2 = x_prev(end);
    u_1= u_prev(ind_phi_left); u_2 = u_prev(end);
    
    x_cur_vec = x_prev((ind_left_key_pt+1) : ind_phi_left);
    phi_cur_vec = u_2 - (u_2 - u_1)./(x_2 - x_1) .* (x_2 - x_cur_vec);
    ind_w_k = find(u_prev((ind_left_key_pt+1) : ind_phi_left) - phi_cur_vec, 1, 'last');
    ind_y_k = ind_left_key_pt + find(x_prev(ind_w_k)*BETA - x_vec, 1, 'last') - 1;
    
end