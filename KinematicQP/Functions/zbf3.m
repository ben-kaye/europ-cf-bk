function [Abf, ubf, h1] = zbf3(x, v_last, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];

    sin2 = 2*c*s;
    cos2 = c^2 - s^2;    
    
    h1 = 2*zz - (v_dir'*z)^2 - 2*delta^2;
    
    safe_vel = 0.5;
    h2 = (v_last-safe_vel)^2;

    omeg = (z(1)^2-z(2)^2)*sin2 - 2*z(1)*z(2)*cos2;
    
    
    Lfh = [ 2*v_last*z'*v_dir; 0 ] ;
    Lgh = [ 0, omeg; 2, 0 ];

    alpha = h1;
    
    % Z(linear)BF
    ubf = gamma*[ h1; h2 ] + Lfh;
    Abf = -Lgh;
    
end