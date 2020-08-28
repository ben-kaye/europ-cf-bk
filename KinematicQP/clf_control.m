function [v, omega] = clf_control(x, r, k1, k2)
    p = x([1,2]);
    phi = x(3);
    p_r = r([1,2]);
    phi_r = r(3);
    phi_rdot = r(4);
    v_r = r(5);
    v_rdot = r(6);
    
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p);
    
    e1 = e12(1);
    e2 = e12(2);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;

    omega = phi_rdot + 2*e2*v_r*cos(e3_aux/2-alpha) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2) + k2*sin(e3_aux/2);
    v = v_r*cos(e3) - omega*v_r*sin(e3_aux/2)/(v_r^2 + e2^2) + k1*e1;
end