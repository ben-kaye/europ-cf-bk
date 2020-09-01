%%% SIMULATION PARAMETERS %%%

sim_time = 3; % {s}
step_size = 1e-3; % {s}
Ts = 1e-2; % {s}

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

x = [ 0.1; 0.1; 40/180*pi ]; % [ p_x, p_y, phi ]
r = [ 0.1; 0.1; 45/180*pi; -0.5; 1; 0];  % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
p_o = [ 1; 0.5 ]; % {m, m}

delta = 0.1; % {m}
v_min = 0; % {ms-1}
v_max = 3; % {ms-1}
% omeg_max = 1.5; % {rads-1}
omeg_max = 3;
q_gamma = 100;

k1 = 5;
k2 = 5;

%%% QP SETUP %%%

lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

[ Abf, ubf] = rbf2(x, ctrls(1), p_o, delta, q_gamma);

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
        errs(e) = 1;
    else
        %%% SOLVED %%%
        ctrls = ctrl_sol;
    end
    
    %%% SUB-STEP FORWARD EULER INTEGRATION
    
    for s = 1:floor(N/Ns)
        [x, r] = sim_xr(x, ctrls, r, step_size);
        r = ref_lookup(r,Ts*e,step_size,1);
    end
    
    %%% STORE RESULTS %%%
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    %%% UPDATE %%%
    
    [ Abf, ubf, h] = rbf2(x, ctrls(1), p_o, delta, q_gamma);
    
    h_t(e) = h;

    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end
errc = sum(errs);
fprintf('Simulation complete. %d errors.\n', errc);

%%% PLOT %%%

plot_ctrls
plot_res
% animate_pos




