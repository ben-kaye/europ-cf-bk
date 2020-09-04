function u = cbf_qp_controller(x, r, p_o, v_last, params)
%CBF_QP_CONTROLLER for kinematic drone model
%   requires params containing k1, k2, ctrl_min(2x1), ctrl_max(2x1),
%   H(2x2), options(quadprog object), delta, gamma

    % update
    ctrls_clf = clf_control(x, r, params.k1, params.k2);
    ctrls_clf = sat_ctrls(ctrls_clf, params.ctrl_min, params.ctrl_max);
    
    f = -params.H'*ctrls_clf;   
    
    [ A, b, ~ ] = sf_constraints(x, v_last, params.ctrl_max(2), p_o, params.delta, params.gamma);
        
    % solve
    u = quadprog(H, f, A, b, [], [], [], [], [], params.options);
    u = sat_ctrls(u, params.ctrl_min, params.ctrl_max);
end

function [Abf, ubf, h] = sf_constraints(x, v_last, max_v, p_o, delta, gamma)
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];
    a_dir = [ -s; c ];

    alpha = atan2(z(2), z(1));
    c_h = cos((alpha - phi)/2);  
    
    h = zz*(2 + c_h) - 2*delta^2;
    
    a1 = z/sqrt(zz - (z'*v_dir)^2);
    
    Lfh = z'*v_dir*v_last*(2 + c_h) + a1*v_last;
    Lgh = [ 0, a1*z'*a_dir ];
    
    uv = -max_v*min(1, 1/(h + 1));
    
    ubf = [ gamma*h + Lfh; uv ];
    Abf = [ -Lgh; -1 0 ];
    
end

function ctrls_clf = clf_control(x, r, k1, k2)
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
    
    ctrls_clf = [ v; omega ];
end

function ctrls_out = sat_ctrls(ctrls_in, min_ctrl, max_ctrl)
    ctrls_out = max(ctrls_in, min_ctrl);
    ctrls_out = min(ctrls_out, max_ctrl);
end

