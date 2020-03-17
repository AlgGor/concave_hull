%% Recurrence Equations diploma
%   All abreviations were used according to A.Y.Zanockin and S.N.Smirnov
%   article "Guaranteed deterministic approach to superhedging: 
%               properties of binary european option"

%% Params
clear;
BETA = 2;                       % greater than 1 
N_POINTS = 1000;          % points per 1 unit
ALPHA = 1/BETA;           % less than 1  

%% Creating initial grids
x_spec_vec = BETA .^ (-1:2); 

init_grid = linspace(x_spec_vec(1), x_spec_vec(2), ceil(N_POINTS * (x_spec_vec(2) - x_spec_vec(1)))+1);
ind_spec = 1:4;
ind_spec(2) = numel(init_grid);
init_grid = [init_grid(1:end-1), ...
                    linspace(x_spec_vec(2), x_spec_vec(3), ceil(N_POINTS * (x_spec_vec(3) - x_spec_vec(2)))+1)];
ind_spec(3) = numel(init_grid);
init_x_vec = [init_grid(1:end-1), ...
                linspace(x_spec_vec(3), x_spec_vec(4), ceil(N_POINTS * (x_spec_vec(4) - x_spec_vec(3)))+1)];
ind_spec(4) = numel(init_x_vec);
clear init_grid x_spec_vec;
init_ind_spec = ind_spec; 
init_u_vec = double(init_x_vec <= 1);


%% Concave shell construction

N_STEPS = 5;

u_prev_vec = init_u_vec;
x_vec = init_x_vec;
ind_spec = init_ind_spec;
n_add_pts = 0;
disp([newline, newline]);

close all;
graph_raw(x_vec, u_prev_vec, BETA, 0, ind_spec);
 
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
                disp(['step=',num2str(step),', k=', num2str(k), ', simple 1 type']);
                u_new_vec((ind_phi_left+1) : ind_phi_right) = first_case( x_vec(ind_phi_left : ind_phi_right), ...
                                                                                                     u_prev_vec(ind_phi_left : ind_phi_right));
            case 2
                % y'(b^(k-1)-0) > slone_coeff ?
                is_y_k_exists = check_y_k_ex(slone_coeff, ...
                                                        x_vec( (ind_phi_left - 2) : ind_phi_left), ...
                                                        u_prev_vec( (ind_phi_left - 2) : ind_phi_left));
                      
                if ~is_y_k_exists
                    disp(['step=',num2str(step),', k=', num2str(k), ', simple 2 type']);
                    u_tmp = second_case( ind_left_key_pt, ind_phi_left,  ...
                                                    x_vec(ind_left_key_pt : ind_phi_right), ... 
                                                    u_prev_vec(ind_left_key_pt : ind_phi_right), BETA);
                    u_new_vec((ind_phi_left+1) : ind_phi_right) = u_tmp;
                else        % need to find y_k
                    disp(['step=',num2str(step),', k=', num2str(k), ', SOPHISTICATED 2 type']);
                    
                    ind_y_k = find_y_k(ind_left_key_pt, ind_phi_left, ...
                                                    x_vec(ind_left_key_pt : ind_phi_right), ...
                                                    u_prev_vec(ind_left_key_pt : ind_phi_right), BETA);                
                    u_tmp = second_case(  ind_left_key_pt, ind_y_k,  ...
                                                        x_vec(ind_left_key_pt : ind_phi_right), ... 
                                                        u_prev_vec(ind_left_key_pt : ind_phi_right), BETA); 
                    u_new_vec((ind_phi_left+1) : ind_y_k ) = u_tmp;
                    u_new_vec((ind_y_k+1) : ind_phi_right) = first_case( x_vec(ind_y_k : ind_phi_right ), ...
                                                                                                    u_prev_vec(ind_y_k : ind_phi_right));
                end

            case 3      
                % y'(b^(k)+0) < slone_coeff ?
                is_z_k_exists = check_z_k_ex(slone_coeff, ...
                                                        x_vec( ind_phi_right : (ind_phi_right+2)), ...
                                                        u_prev_vec( ind_phi_right : (ind_phi_right+2)));    
                            
                if ~is_z_k_exists
                    disp(['step=',num2str(step),', k=', num2str(k), ', simple 3 type']);
                    u_tmp = third_case( ind_phi_left, ind_phi_right,  ...
                                                    x_vec(ind_phi_left : ind_right_key_pt), ... 
                                                    u_prev_vec(ind_phi_left : ind_right_key_pt), BETA);
                    u_new_vec((ind_phi_left+1) : ind_phi_right) = u_tmp;
                else
                    disp(['step=',num2str(step),', k=', num2str(k), ', SOPHISTICATED 3 type']);
                    ind_z_k = find_z_k(ind_phi_left, ind_phi_right, ...
                                                    x_vec(ind_phi_left : ind_right_key_pt), ...
                                                    u_prev_vec(ind_phi_left : ind_right_key_pt), BETA);
                     u_new_vec((ind_phi_left+1) : ind_z_k) = first_case( x_vec(ind_phi_left : ind_z_k), ...
                                                                                                        u_prev_vec(ind_phi_left : ind_z_k));
                    u_tmp = third_case(  ind_z_k, ind_phi_right,  ...
                                                        x_vec(ind_z_k : ind_right_key_pt), ... 
                                                        u_prev_vec(ind_z_k : ind_right_key_pt), BETA); 
                    u_new_vec((ind_z_k+1): ind_phi_right) = u_tmp;                                                                            
                end

            case 4
                disp(['step=',num2str(step),', k=', num2str(k), ', SOPHISTICATED 4 type']);
                is_y_k_exists = check_y_k_ex(slone_coeff, ...
                                                        x_vec( (ind_phi_left - 2) : ind_phi_left), ...
                                                        u_prev_vec( (ind_phi_left - 2) : ind_phi_left));
                is_z_k_exists = check_z_k_ex(slone_coeff, ...
                                                        x_vec( ind_phi_right : (ind_phi_right+2)), ...
                                                        u_prev_vec( ind_phi_right : (ind_phi_right+2)));   
                if is_y_k_exists
                    ind_y_k = find_y_k(ind_left_key_pt, ind_phi_left, ...
                                                    x_vec(ind_left_key_pt : ind_phi_right), ...
                                                    u_prev_vec(ind_left_key_pt : ind_phi_right), BETA);   
                end
                if is_z_k_exists
                    ind_z_k = find_z_k(ind_phi_left, ind_phi_right, ...
                                                    x_vec(ind_phi_left : ind_right_key_pt), ...
                                                    u_prev_vec(ind_phi_left : ind_right_key_pt), BETA);
                end
                if ~is_y_k_exists && ~is_z_k_exists
                    u_new_vec((ind_phi_left+1) : ind_phi_right) = fourth_case(ind_left_key_pt, ind_phi_left, ind_phi_right, ...
                                                                                                            x_vec(ind_left_key_pt : ind_right_key_pt), ...
                                                                                                            u_prev_vec(ind_left_key_pt : ind_right_key_pt));
                elseif is_y_k_exists && ~is_z_k_exists
                    
                elseif ~is_y_k_exists && is_z_k_exists
                    
                elseif is_y_k_exists && is_z_k_exists
                    ind_left = min(ind_y_k, ind_z_k); 
                    ind_right = max(ind_y_k, ind_z_k);
                    first_case_appears = ind_z_k < ind_y_k;
                    %u_new_vec((ind_phi_left+1):ind_left) = 
                    
                    
                    
                end 
        end
        
    end
    graph_raw(x_vec, u_new_vec, BETA, step, ind_spec);
    
    tmp_n_pts =  ceil(N_POINTS * (BETA^(step+2) - x_vec(end)));
    x_vec = cat(2, x_vec(1:end-1), linspace(x_vec(end), BETA^(step+2), tmp_n_pts+1));
    u_prev_vec = cat(2, u_new_vec, 0 * (1 : tmp_n_pts));         
    ind_spec = cat(2, ind_spec, numel(x_vec));
    
    % не забудь сделать корректировку индексов массива ключевых точек с
    % учетом добавленных на каждом шаге
    n_add_pts = n_add_pts + n_tmp_pts;
end

%%  Testing app version
clear;
close all;

BETA = 2;       
N_STEPS = 3;

[x_vec, u_mat, ind_spec] = conc_sh_bld(N_STEPS, BETA);

for ind_step = 1:numel(u_mat(:,1))
    graph_raw(x_vec, u_mat(ind_step,:), BETA, ind_step-1, ind_spec);
end

        




   