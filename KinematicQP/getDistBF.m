function [ Abf, ubf ] = getDistBF(x, p_o, delta, gamma)
    p = x([1,2]);
    phi = x(3);
       
    z = p_o - p;
    
    h = z'*z - delta^2;
    
    vdir = [ cos(phi); sin(phi) ];
    
    LgB = [ -2*z'*vdir, 0 ];
    LfB = 0;
    
    Abf = LgB;
    ubf = gamma*h - LfB;

end

