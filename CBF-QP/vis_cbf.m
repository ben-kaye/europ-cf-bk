vis(x([1,2]), p_o, x(3), v_old, max_u, delta);

function vis(p, p_o, phi, v, max_turn, delta)

    vecv = [ cos(phi); sin(phi) ];

    p_dot = v * vecv;
    
    p_xo = p_o - p;
    
    sign_term = sign(p_dot' * [ 0, 1; -1, 0 ] * p_xo); %if innerprod > 0;

    Rmin = v/max_turn;
    
    % cross(p_dot, k)
    r = sign_term*[ 0, 1; -1, 0 ] * p_dot / max_turn;
    
    p_c = p + r;
    
    vel_vis = [ p, p + 2*vecv ];
    
    z = p_o - p_c;
    
    h = z'*z - (v/max_turn + delta)^2;
    
    figure(1)
    hold off
    viscircles( p_c', Rmin, 'LineStyle', '-.', 'Color', 1/255*[242, 113, 7]);


    hold on
    viscircles( p_o', delta, 'LineStyle', '-.');
    quiver( p(1),p(2),1*vecv(1),1*vecv(2), '-', 'LineWidth',1,'Color',1/255*[223, 58, 252], 'DisplayName', 'Velocity Vector');
    plot( p_c(1), p_c(2), '+', 'MarkerSize', 20, 'DisplayName', 'Manoeuvre Centre');
    plot( p_o(1), p_o(2), '+', 'MarkerSize', 20, 'DisplayName', 'Object');
    plot( p(1),p(2), '+', 'MarkerSize', 20, 'DisplayName', 'Position');
    
    title(sprintf('h = %f',h));

    axis equal
    legend()
    hold off
end
