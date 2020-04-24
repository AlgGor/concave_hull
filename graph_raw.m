function graph_raw(x_vec, y_vec, BETA, step, ind_spec, ind_add, varargin)
    
    if step>0
        ZERO_PART = 0.3820;
    else
        ZERO_PART = BETA;
    end
    
    if nargin> 6
        app = varargin{1};
        my_ax = app.UIAxes;
    else
        figure();
        my_ax = gca;
    end 
    blackcolor = [0 0 0]; %rand(1,3); %[0 0 0]; 
    
    x_leap = BETA^step;
    leap_ind = find(y_vec, 1, 'last');
    edge_ind_spec = find(ind_spec == leap_ind, 1); 
    edge_ind_add = find(ind_add < leap_ind, 1); 
    y_leap = y_vec(leap_ind);
    y_vec(leap_ind) = NaN;  
    
    % main part and leap
    plot(my_ax, x_vec, y_vec, 'linewidth', 1.5, 'color', blackcolor);
    xlim(my_ax,[1/BETA BETA^step * min(1+ZERO_PART, BETA + 0.01)]); % you can make BETA + 0.05 in the end to add x axe arrow
    hold(my_ax, 'on');
    plot(my_ax,  [x_vec(leap_ind)  BETA^(step+2)], [0 0], 'linewidth', 1.5, 'color', blackcolor);
    plot(my_ax, [0, BETA^(-1)], [1,1], 'linewidth', 1.5, 'color', blackcolor);
    
    plot(my_ax, [x_leap x_leap], [0 y_leap], '--', ...
                            'color', blackcolor, 'linewidth', 0.1);
    plot(my_ax, x_vec((leap_ind-1):leap_ind), [y_vec(leap_ind-1), y_leap], 'linewidth', 1.5, 'color', blackcolor);
    
    % dots
    scatter(my_ax, x_leap, y_leap, 'o', 'filled', 'MarkerFaceColor', blackcolor);
    scatter(my_ax, x_leap, 0, 'MarkerEdgeColor', blackcolor, 'MarkerFaceColor', [1 1 1]);
    scatter(my_ax, x_vec(ind_spec(1:edge_ind_spec)), y_vec(ind_spec(1:edge_ind_spec)) , 12, 'o', 'filled', 'MarkerFaceColor', blackcolor);
    if ~isempty(edge_ind_add) 
        scatter(my_ax, x_vec(ind_add(1:edge_ind_add)), y_vec(ind_add(1:edge_ind_add)) , 12, 'o', 'filled', 'MarkerFaceColor', 'r');
    end
    legend(my_ax, ['f_{',num2str(step),'}(x)'], 'location','northeast');

    % axes meta to be pretty
    my_ax.TickLabelInterpreter = 'latex';
    my_ax.MinorGridAlpha = 0.17;
%     [left bottom width height] = my_ax.Position;
%     my_ax.FontUnits = 'normalized';
%    my_ax.FontSize = 12 - fix(step/2);
    my_ax.TickDir = 'both';
    box(my_ax,'off');
    grid(my_ax, 'on');
    grid(my_ax, 'minor');
    
    hold(my_ax, 'off');
end