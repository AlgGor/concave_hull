function u_vec = third_case(ind_phi_left, ind_phi_right, x_prev, u_prev, BETA)
    
    ind_phi_right = ind_phi_right - ind_phi_left + 1;
    x_1 = x_prev(1);
    u_1 = u_prev(1);
    
    x_cur_vec = x_prev(2 : ind_phi_right);
    u_cur_vec = interp1([x_prev, 2*x_prev(end)], [u_prev, u_prev(end)], x_cur_vec * BETA);
    
    u_vec = u_cur_vec - (u_cur_vec - u_1)./(x_cur_vec * BETA - x_1) .* x_cur_vec * (BETA - 1);
    
end