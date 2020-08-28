%%% using quadprog %%%

sim_time = 10;
step_size = 1e-4;
Ts = 1e-3;

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);


x = [ 1; 1; -3*pi/18 ]; % { p_x, p_y, phi }
r = [ 2; 3; pi/3; 1; 1; 0]; % { r_x, r_y, phi_r, phi_rdot, v_r, v_rdot }
p_o = [ -1; 6 ];
delta = 0.4;
v_min = 1e-3;
v_max = 4;
omeg_max = 1.5;


lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

% (x, v_last, max_turn, p_o, delta, gamma)
[ Arcbf, urcbf ] = getInverseB(x, ctrls(2), max_omeg, p_o, delta, gamma);


H = diag([1, 1]);
f = -H'*ctrls;

x_t = zeros(2, Ns);
r_t = zeros(2, Ns);
u_t = zeros(2, Ns);

for e = 1:Ns
    ctrls = quadprog(H, f, Arcbf, urcbf, [], [], lx, ux);
    
    for s = 1:floor(N/Ns)
        [x, r] = sim_xr(x, ctrls, r, step_size);
    end
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    [ Arcbf, urcbf ] = getInverseB(x, ctrls(2), max_omeg, p_o, delta, gamma);

    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end

