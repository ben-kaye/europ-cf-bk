%%% SIMULATION PARAMETERS %%%

sim_time = 10; % {s}
step_size = 1e-3; % {s}
Ts = 1e-1; % {s}

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

x = [ 3; 2; -3*pi/18 ]; % [ p_x, p_y, phi ]
r = [ 2; 3; pi/3; 0.5; 1; 0]; % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
p_o = [ -1; 6 ]; % {m, m}

delta = 0.8; % {m}
v_min = 1e-3; % {ms-1}
v_max = 4; % {ms-1}
omeg_max = 1.5; % {rads-1}
q_gamma = 3;

k1 = 5;
k2 = 5;

%%% QP SETUP %%%

lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

[ Abf, ubf] = augmented_zbf(x, p_o, delta, q_gamma);

H = diag([1, 1]);
f = -H'*ctrls;

%%% SIMULATION %%%

x_t = NaN*ones(2, Ns);
r_t = NaN*ones(2, Ns);
u_t = NaN*ones(2, Ns);
h_t = NaN*ones(1, Ns);

options =  optimset('Display','off');

errc = 0;

for e = 1:Ns
    %%% SOLVE %%%
    
    ctrl_sol = quadprog(H, f, Abf, ubf, [], [], lx, ux, [], options);
    
    if isempty(ctrl_sol)
        %%% INFEASIBLE %%%
        fprintf('err\n');
        errc = ercc + 1;
    else
        %%% SOLVED %%%
        ctrls = ctrl_sol;
    end
    
    %%% SUB-STEP FORWARD EULER INTEGRATION
    
    for s = 1:floor(N/Ns)
        [x, r] = sim_xr(x, ctrls, r, step_size);
    end
    
    %%% STORE RESULTS %%%
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    %%% UPDATE %%%
    
    [ Abf, ubf, h] = augmented_zbf(x, p_o, delta, q_gamma);
    
    h_t(e) = h;

    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end

%%% PLOT %%%

plot_res
animate_pos

t = 0:step_size*N/Ns:sim_time;
t = t(1:end-1);
figure(3)
plot(t, h_t);
yline(0, 'LineStyle','-.','Color',[1,0,0])