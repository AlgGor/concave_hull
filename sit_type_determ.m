function situation_type = sit_type_determ(is_phi_gr_left, is_phi_gr_right)

     if ~is_phi_gr_left && ~is_phi_gr_right
        situation_type = 1; 
     end  
     if is_phi_gr_left && ~is_phi_gr_right
        situation_type = 2; 
     end
    if ~is_phi_gr_left && is_phi_gr_right
        situation_type = 3; 
    end
     if is_phi_gr_left && is_phi_gr_right
        situation_type = 4; 
    end  
    
end