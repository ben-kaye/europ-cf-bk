function [Abf, ubf, h] = sf2(x, v_last, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];

    b = v_dir'*z;
%     if sign(b) > 0
%         b = max(b, delta);
%     end
  
    alpha = atan2(z(2), z(1));
    half_beta = (alpha - phi)/2;
    
    sin2 = 2*c*s;
    cos2 = c^2 - s^2;    
    
    h = zz*(2 + cos(half_beta)) - 2*delta^2;

    omeg = (z(1)^2-z(2)^2)*sin2 - 2*z(1)*z(2)*cos2;
    
    
    Lfh = 2*v_last*z'*v_dir;
    Lgh = [ 0, omeg ];
    
    % Z(linear)BF
    ubf = gamma*h + Lfh;
    Abf = -Lgh;
    
end