function plots = plot_drone(p_x, p_theta, rotors, fig)
%PLOT_DRONE Summary of this function goes here
%   Detailed explanation goes here
    r = p_theta(1);
    p = p_theta(2);
    y = p_theta(3);

    iRb = [ cos(p)*cos(y),  -cos(r)*sin(y) + sin(r)*sin(p)*cos(y),  sin(r)*sin(y) + cos(r)*sin(p)*cos(y);
            cos(p)*sin(y),  cos(r)*cos(y) + sin(r)*sin(p)*sin(y),   -sin(r)*cos(y) + cos(r)*sin(p)*sin(y);
            -sin(p),        cos(p)*sin(r),                          cos(p)*cos(r)                           ];
    
        
    rotors_transform = zeros(12, size(rotors, 2));
    for z = 1:4 
        rotors_transform(3*z-2:3*z,:) = iRb*[ rotors([2*z-1, 2*z], :); zeros(1, size(rotors, 2)) ] + p_x(1:3);
    end
    
    plots = gobjects(4);
    
    figure(fig);
    hold on
    for z = 1:4
        path = rotors_transform(3*z-2:3*z, :);
        if z == 1
            c = 'b';
        else 
            c = 'k';
        end
        plots(z) = plot3(path(1,:), path(2,:), path(3,:), c, 'LineWidth', 2);
    end
    hold off
    

end
