function u_vec = fourth_case(ind_left_key_pt, ind_phi_left, ind_phi_right, x_prev, u_prev, BETA)
    
    ind_phi_left = ind_phi_left - ind_left_key_pt + 1;
    ind_phi_right =  ind_phi_right - ind_left_key_pt + 1;
    
    x_cur_vec = x_prev((ind_phi_left+1) : ind_phi_right);
    u_cur_left_vec = interp1([x_prev(1)/2, x_prev], [u_prev(1), u_prev], x_cur_vec/BETA);
    u_cur_right_vec = interp1([x_prev, x_prev(end)*2], [u_prev, u_prev(end)], x_cur_vec*BETA);
    
    u_vec = u_cur_right_vec - (u_cur_right_vec- u_cur_left_vec)./(x_cur_vec*BETA - x_cur_vec/BETA) .* (x_cur_vec*BETA - x_cur_vec);
    
end