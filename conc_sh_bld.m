function [x_vec, u_mat, ind_spec, ind_add] = conc_sh_bld(N_STEPS, BETA) 

    N_POINTS = 100; 
    
    x_spec_vec = BETA .^ (-1:2); 

    init_grid = linspace(x_spec_vec(1), x_spec_vec(2), ceil(N_POINTS * (x_spec_vec(2) - x_spec_vec(1)))+1);
    ind_spec = 1:4;
    ind_spec(2) = numel(init_grid);
    init_grid = [init_grid(1:end-1), ...
                        linspace(x_spec_vec(2), x_spec_vec(3), ceil(N_POINTS * (x_spec_vec(3) - x_spec_vec(2)))+1)];
    ind_spec(3) = numel(init_grid);
    x_vec = [init_grid(1:end-1), ...
                    linspace(x_spec_vec(3), x_spec_vec(4), ceil(N_POINTS * (x_spec_vec(4) - x_spec_vec(3)))+1)];
    ind_spec(4) = numel(x_vec);
    clear init_grid x_spec_vec;
    u_prev_vec = double(x_vec <= 1);
    
    u_mat = zeros(N_STEPS+1, ceil(N_POINTS * (BETA^(N_STEPS+2) - x_vec(1))) + N_STEPS);
    u_mat(1, 1:numel(u_prev_vec)) = u_prev_vec;
    n_add_pts = 0;
    

    for step = 1 : N_STEPS
        u_new_vec = u_prev_vec;
        for k = 1 : (step + n_add_pts)
            n_tmp_pts = 0;
            ind_bk = k + 2;                                                          
            ind_left_key_pt = ind_spec(ind_bk - 2);                           % index of \beta^(k-2) in x_vec
            ind_phi_left = ind_spec(ind_bk - 1);                                % index of \beta^(k-1)   in x_vec     
            ind_phi_right = ind_spec(ind_bk);                                   % index of \beta^(k)   in x_vec
            ind_right_key_pt = ind_spec(ind_bk + 1);                        % index of \beta^(k+1) in x_vec

            [slone_coeff, phi_left, phi_right] = key_phi_line(  x_vec(ind_left_key_pt+1), x_vec(ind_right_key_pt), ...
                                                                                    x_vec(ind_phi_left), u_prev_vec(ind_phi_left), ...
                                                                                    x_vec(ind_phi_right), u_prev_vec(ind_phi_right));
            % is key_line less on the left side with taking into account the numeric error  
            is_phi_less_left = (phi_left < u_prev_vec(ind_left_key_pt+1)) && (abs(u_prev_vec(ind_left_key_pt+1) - phi_left)> 2*eps);     
            % is key_line less on the right side with taking into account the numeric error  
            is_phi_less_right =(phi_right < u_prev_vec(ind_right_key_pt)) && (abs(u_prev_vec(ind_right_key_pt) - phi_right)> 2*eps);              
            situation_type = sit_type_determ(is_phi_less_left, is_phi_less_right);

            switch situation_type 
                case 1
                    u_new_vec( (ind_phi_left+1) : ind_phi_right) = first_case( x_vec(ind_phi_left : ind_phi_right), ...
                                                                                                         u_prev_vec(ind_phi_left : ind_phi_right));
                    % disp(['step=',num2str(step),', k=', num2str(k), ', simple 1 type']);
                case 2
                    % y'(b^(k-1)-0) > slone_coeff ?
                    is_y_k_exists = check_y_k_ex(slone_coeff, ...
                                                            x_vec( (ind_phi_left - 2) : ind_phi_left), ...
                                                            u_prev_vec( (ind_phi_left - 2) : ind_phi_left));  

                    if ~is_y_k_exists
                        u_tmp = second_case( ind_left_key_pt, ind_phi_left,  ...
                                                        x_vec(ind_left_key_pt : ind_phi_right), ... 
                                                        u_prev_vec(ind_left_key_pt : ind_phi_right), BETA);
                        u_new_vec( (ind_phi_left+1) : ind_phi_right) = u_tmp;
                        %disp(['step=',num2str(step),', k=', num2str(k), ', simple 2 type']);
                    else        % need to find y_k
                        %disp('lol_2');
                        %ind_y_k = find_y_k();
                    end
                case 3      
                    % y'(b^(k)+0) < slone_coeff ?
                    is_z_k_exists = check_z_k_ex(slone_coeff, ...
                                                            x_vec( ind_phi_right : (ind_phi_right+2)), ...
                                                            u_prev_vec( ind_phi_right : (ind_phi_right+2)));    
                    if ~is_z_k_exists
                        u_tmp = third_case( ind_phi_left, ind_phi_right,  ...
                                                        x_vec(ind_phi_left : ind_right_key_pt), ... 
                                                        u_prev_vec(ind_phi_left : ind_right_key_pt), BETA);
                        u_new_vec( (ind_phi_left+1) : ind_phi_right) = u_tmp;
                        %disp(['step=',num2str(step),', k=', num2str(k), ', simple 3 type']);
                    else
                        %disp('lol_3');
                    end
                case 4
                    %disp('lol_4');
            end

        end
        
        tmp_n_pts =  ceil(N_POINTS * (BETA^(step+2) - x_vec(end)));
        u_mat(step+1,1:ind_phi_right) =  u_new_vec(1:ind_phi_right);
        x_vec = cat(2, x_vec(1:end-1), linspace(x_vec(end), BETA^(step+2), tmp_n_pts+1));
        u_prev_vec = cat(2, u_new_vec, 0 * (1 : tmp_n_pts));         
        
        ind_spec = cat(2, ind_spec, numel(x_vec));

        % не забудь сделать корректировку индексов массива ключевых точек с
        % учетом добавленных на каждом шаге
        n_add_pts = n_add_pts + n_tmp_pts;
    end
    
    u_mat(:, ind_phi_right + 1:end) = [];
    x_vec(ind_phi_right +1:end) = [];
    ind_spec = ind_spec(1:end-1);
    ind_add = [];
    
end






