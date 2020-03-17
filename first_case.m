function u_vec = first_case(x_prev, u_prev)
    
    x_1 = x_prev(1); x_2 = x_prev(end);
    u_1 = u_prev(1); u_2 = u_prev(end);
    
    u_vec = u_2 + (u_1 - u_2)/(x_2 - x_1) * (x_2 - x_prev(2:end));
        
end