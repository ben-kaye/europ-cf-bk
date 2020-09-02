function [x_t, u_t, r_t, h_t, errs] = bf_qp(BF, sim_time, step_size, Ts, x0, r0, path_id, p_o, delta, v_min, v_max, omeg_max, gamma, k1, k2)

%%% SIMULATION PARAMETERS %%%

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

x = x0; % [ p_x, p_y, phi ]
r = r0;  % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
% p_o = settings.po; % {m, m}

%%% QP SETUP %%%

lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

[ Abf, ubf] = BF(x, ctrls(1), p_o, delta, gamma);

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
        r = ref_lookup(r, Ts*e, step_size, path_id);
    end
    
    %%% STORE RESULTS %%%
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    %%% UPDATE %%%
    
    [ Abf, ubf, h] = BF(x, ctrls(1), p_o, delta, gamma);
    
    h_t(e) = h;

    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end
errc = sum(errs);
fprintf('Simulation complete. %d errors.\n', errc);

end