%%% using quadprog %%%

sim_time = 10;
step_size = 1e-2;
Ts = 1e-1;

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);


x = [ 3; 2; -3*pi/18 ]; % { p_x, p_y, phi }
r = [ 2; 3; pi/3; 0.5; 1; 0]; % { r_x, r_y, phi_r, phi_rdot, v_r, v_rdot }
p_o = [ -1; 6 ];
delta = 0.4;
v_min = 1e-3;
v_max = 4;
omeg_max = 1.5;
q_gamma = 0.5;

k1 = 5;
k2 = 5;


lx = [ v_min; -omeg_max ];
ux = [ v_max; omeg_max ];

[v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
ctrls = [v_ctrl; omeg_ctrl];
ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);

% (x, v_last, max_turn, p_o, delta, gamma)
[ Abf, ubf] = augmentedBF(x, ctrls(1), p_o, delta, q_gamma);

H = diag([1, 1]);
f = -H'*ctrls;

x_t = zeros(2, Ns);
r_t = zeros(2, Ns);
u_t = zeros(2, Ns);
h_t = zeros(1, Ns);

options =  optimset('Display','off');

for e = 1:Ns
    ctrl_sol = quadprog(H, f, Abf, ubf, [], [], lx, ux, [], options);
    if isempty(ctrl_sol)
        fprintf('err\n');
    else
        ctrls = ctrl_sol;
    end
    
    for s = 1:floor(N/Ns)
        [x, r] = sim_xr(x, ctrls, r, step_size);
    end
    
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);
    u_t(:,e) = ctrls;
    
    [ Abf, ubf, h] = augmentedBF(x, ctrls(1), p_o, delta, q_gamma);
    
    if h > 50
        h_t(e) = NaN;
    else 
        h_t(e) = h;
    end
    [v_ctrl, omeg_ctrl] = clf_control(x, r, k1, k2);
    ctrls = [v_ctrl; omeg_ctrl];
    ctrls = sat_ctrls(ctrls, [v_min; -omeg_max], [v_max; omeg_max]);
    
    f = -H'*ctrls;
end

plot_res

t = 0:step_size*N/Ns:sim_time;
t = t(1:end-1);
figure(3)
plot(t, h_t);

