function [Abf, ubf] = getInverseB(x, v_last, max_turn, p_o, delta, gamma)
    phi = x(3);
    p = x([1,2]);

    v_dir = [ cos(phi); sin(phi) ];

    p_xo = p_o - p;
    
    sign_r = sign(v_dir' * [ p_xo(2); -p_xo(1) ]); %if innerprod > 0;

    R_min = v_last/max_turn;
    
    % cross(p_dot, k)
    r = R_min*sign_r*[ v_dir(2); -v_dir(1) ];
    z = p_xo - r;
    
    h = z'*z - (R_min + delta)^2;

    LfB = 2*z'*v_last*v_dir/h^2;
    
    LgB = 2*z'*R_min*v_dir*sign_r/h^2;
    
    Abf = [ 0, LgB ];
    ubf = gamma*h - LfB;    
end