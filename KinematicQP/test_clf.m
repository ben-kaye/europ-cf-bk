function [v, omega] = test_clf(x, r, k1, k4)
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
    alpha = atan(v_r*e2);
    
    e4 = -e3 - alpha;    

    omega = phi_rdot + (v_r*sin(e3) + v_rdot*e2)/(1+(v_r*e2)^2) - k4*sin(e4/2);
    v = v_r*cos(e3) + omega*v_r*sin(e4/2)/(1+(v_r*e2)^2) + k1*e1;
end