function [Abf, ubf, h] = zbf2(x, v_last, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];

    sin2 = 2*c*s;
    cos2 = c^2 - s^2;    
    
    h = 2*zz - (v_dir'*z)^2 - 2*delta^2;

    omeg = (z(1)^2-z(2)^2)*sin2 - 2*z(1)*z(2)*cos2;
    
    
    Lfh = 2*v_last*z'*v_dir;
    Lgh = [ 0, omeg ];

    alpha = h;
    
    % Z(linear)BF
    ubf = gamma*alpha + Lfh;
    Abf = -Lgh;
    
end