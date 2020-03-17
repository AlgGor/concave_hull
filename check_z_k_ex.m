function  answer = check_z_k_ex(slone_coeff, x_vec, u_vec)
    
    der =  (u_vec(3) - u_vec(1)) / (x_vec(3) - x_vec(1));
    % slone_coef <0 and der<0 hence der > slone_coeff equal 
    % abs(der) < abs(slone_coef) and it is similar to the original article
    
    answer = (der < slone_coeff) && (abs(slone_coeff- der) > 2*eps);   
    
end