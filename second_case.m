function u_vec = second_case(ind_left_key_pt, ind_phi_left, x_prev, u_prev, BETA)
    
    ind_phi_left = ind_phi_left - ind_left_key_pt + 1;
    x_2 = x_prev(end);
    u_2 = u_prev(end);
    
    x_cur_vec = x_prev((ind_phi_left+1) : end);
    u_cur_vec = interp1(x_prev, u_prev, x_cur_vec/BETA);
    
    u_vec = u_2 - (u_2- u_cur_vec)./(x_2 - x_cur_vec/BETA) .* (x_2 - x_cur_vec);
    
end