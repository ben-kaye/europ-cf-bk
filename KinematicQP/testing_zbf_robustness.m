%%% SIMULATION PARAMETERS %%%

sim_time = 3; % {s}
step_size = 1e-3; % {s}
Ts = 1e-1; % {s}

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

x = [ 0.1; 0.1; 40/180*pi ]; % [ p_x, p_y, phi ]
r = [ 0.1; 0.1; 45/180*pi; -0.5; 1; 0];  % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
p_o = [ 1; 0.5 ]; % {m, m}

delta = 0.1; % {m}
v_min = 0; % {ms-1}
v_max = 3; % {ms-1}
omeg_max = 1.5; % {rads-1}
q_gamma = 2.3;

k1 = 5;
k2 = 5;

%%% QP SETUP %%%

lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

[ Abf, ubf] = augmented_zbf(x, ctrls(1), p_o, delta, q_gamma);

H = diag([1, 1]);
f = -H'*ctrls;

%%% SIMULATION %%%

x_t = NaN*ones(2, Ns);
r_t = NaN*ones(2, Ns);
u_t = NaN*ones(2, Ns);
h_t = NaN*ones(1, Ns);

options =  optimset('Display','off');

errc = 0;
errs = logical(zeros(1,Ns));
for e = 1:Ns
    %%% SOLVE %%%
    
    ctrl_sol = quadprog(H, f, Abf, ubf, [], [], lx, ux, [], options);
    
    if isempty(ctrl_sol)
        %%% INFEASIBLE %%%
        fprintf('err\n');
        errs(e) = 1;
    else
        %%% SOLVED %%%
        ctrls = ctrl_sol;
    end
    
    %%% SUB-STEP FORWARD EULER INTEGRATION
    
    for s = 1:floor(N/Ns)
        [x, r] = sim_xr(x, ctrls, r, step_size, 1);
        
        dist2 = sum((x([1,2]) - r([1,2])).^2);
        tau = scale_time(step_size, dist2);
        if e < 6
            r(4) = r(4) - 1.1*step_size;
        else
            if e < 10
                r(4) = -0.9;
                r(6) = -0.05;
            else 
                r(4) = 0;
                r(6) = 0;
            end
        end
    end
    
    %%% STORE RESULTS %%%
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    %%% UPDATE %%%
    
    [ Abf, ubf, h] = augmented_zbf(x, ctrls(1), p_o, delta, q_gamma);
    
    h_t(e) = h;

    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end
errc = sum(errs);

%%% PLOT %%%

plot_ctrls
plot_res
% animate_pos

% override sim_xr
function [x_1, r_1] = sim_xr(x, u, r, step_sz, scale)
    phi = x(3);
    x_dot = [ cos(phi), 0; sin(phi), 0; 0, 1 ] * u;
    
    phi_r = r(3);
    v_r = r(5);
    phi_rdot = r(4);
    v_rdot = r(6);
    r_dot =  [ v_r*[ cos(phi_r); sin(phi_r) ]; phi_rdot; 0; v_rdot; 0 ];
    
    x_1 = x + step_sz*x_dot;
    
    if scale
        dist2 = sum((x([1,2]) - r([1,2])).^2);
        tau = scale_time(step_sz, dist2);
    end
    
    r_1 = r + tau*r_dot;    
end

function [Abf, ubf, h] = augmented_zbf(x, v_last, p_o, delta, gamma)
    
    p = x([1,2]);
    phi = x(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];
    a_dir = [ -s; c ];
    
    h = 2*zz - (v_dir'*z)^2 - 2*delta^2;
    
%     coeff = 2*z'*v_dir;
    
    b = 2*phi;
    omeg = (z(1)^2-z(2)^2)*sin(b) - 2*z(1)*z(2)*cos(b);
    
    
    Lfh = 2*v_last*z'*v_dir;
    Lgh = [ 0, omeg ];

    alpha = h + h^3;
    
    % Z(linear)BF
    ubf = gamma*alpha + Lfh;
    Abf = -Lgh;
    
end



