function [x_1, r_1] = sim_xr(x, u, r, step_sz)
    phi = x(3);
    x_dot = [ cos(phi), 0; sin(phi), 0; 0, 1 ] * u;
    
    phi_r = r(3);
    v_r = r(5);
    phi_rdot = r(4);
    v_rdot = r(6);
    r_dot =  [ v_r*[ cos(phi_r); sin(phi_r) ]; phi_rdot; 0; v_rdot; 0 ];
    
    x_1 = x + step_sz*x_dot;
    r_1 = r + step_sz*r_dot;    
end

