function [Abf, ubf, h] = augmentedBF(x,v_last, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];
    a_dir = v_last*[ -s; c ];
    
    h = 2*zz + 1/v_last^2*v_dir'*z - 2*delta^2;
    
    Lfh = 1;
    
    Lgh = [ z'*v_dir, a_dir'/zz*( 2*p/zz + z ) ];
    
    % Z(linear)BF
    ubf = gamma*h + Lfh;
    Abf = -Lgh;
    
end

