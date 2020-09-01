function [Abf, ubf, h] = augmented_zbf(x, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];
    a_dir = [ -s; c ];
    
    h = 2*zz - (v_dir'*z)^2 - 2*delta^2;
    
%     Lfh = 1;
    
%     Lgh = [ z'*v_dir, v_last*a_dir'/zz*( 2*p/zz + z ) ];
   
    
    coeff = 2*z'*v_dir;
    
    c = 2*phi;
    
    Lfh = 0;
    Lgh = coeff*[ 1, -z'*a_dir ];


    % Z(linear)BF
    ubf = gamma*h + Lfh;
    Abf = -Lgh;
    
end

