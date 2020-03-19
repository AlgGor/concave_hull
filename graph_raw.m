function graph_raw(x_vec, y_vec, BETA, step, ind_spec, ind_add, varargin)
    
    if nargin> 6
        app = varargin{1};
        my_ax = app.UIAxes;
    else
        figure();
        my_ax = gca;
    end 
    blackcolor = [0 0 0];
    
    leap_ind = find(y_vec, 1, 'last');
    edge_ind_spec = find(ind_spec == leap_ind, 1); 
    edge_ind_add = find(ind_add < leap_ind, 1); 
    y_vec(leap_ind) = NaN;  
    
    plot(my_ax, x_vec, y_vec, 'linewidth', 1.5, 'color', blackcolor);
    xlim(my_ax,[max(0, 2/BETA - 1) BETA^(step+1)]);
    hold(my_ax, 'on');
    plot(my_ax,  [x_vec(leap_ind)  BETA^(step+1)], [0 0], 'linewidth', 1.5, 'color', blackcolor);
    plot(my_ax, [0, BETA^(-1)], [1,1], 'linewidth', 1.5, 'color', blackcolor);
    x_leap = BETA^step;
    y_leap = (BETA/(1+BETA))^step;
    plot(my_ax, [x_leap x_leap], [0 y_leap], '--', ...
                            'color', blackcolor, 'linewidth', 0.1);
    scatter(my_ax, x_leap, y_leap, 'o', 'filled', 'MarkerFaceColor', blackcolor);
    scatter(my_ax, x_leap, 0, 'MarkerEdgeColor', blackcolor, 'MarkerFaceColor', [1 1 1]);
    scatter(my_ax, x_vec(ind_spec(1:edge_ind_spec)), y_vec(ind_spec(1:edge_ind_spec)) , 12, 'o', 'filled', 'MarkerFaceColor', blackcolor);
    if ~isempty(edge_ind_add) 
        scatter(my_ax, x_vec(ind_add(1:edge_ind_add)), y_vec(ind_add(1:edge_ind_add)) , 12, 'o', 'filled', 'MarkerFaceColor', 'r');
    end
    legend(my_ax, ['f_{',num2str(step),'}(x)'], 'location','northeast');
    hold(my_ax, 'off');
    
end