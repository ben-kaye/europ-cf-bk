
[a,b]=controlFromCLF([1; 1.5], [2; 3], 3/18*pi, pi/3, 0.5, 1, 0, 3, 5)


function [ v, omega ] = controlFromCLF(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot, k1, k2)
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p_xy);

    e2 = e12(2);
    e1 = e12(1);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;

    omega = phi_rdot + 2*e2*v_r*cos(e3_aux/2-alpha) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2) + k2*sin(e3_aux/2);
    v = v_r*cos(e3) - omega*v_r*sin(e3_aux/2)/(v_r^2 + e2^2) + k1*e1;
end