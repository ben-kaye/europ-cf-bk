% function [x_t, u_t, r_t, h_t] = bf_qp(BF, sim_time, step_size, Ts, x0, r0, path_id p_o, delta, v_min, v_max, omeg_max, gamma, k1, k2)
clear 

sim_time = 3; % {s}
step_size = 1e-3; % {s}
Ts = 1e-2; % {s}

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

x0 = [ 0.1; 0.1; 40/180*pi ]; % [ p_x, p_y, phi ]
r0 = [ 0.1; 0.1; 45/180*pi; -0.5; 1; 0];  % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
p_o = [ 1; 0.5 ]; % {m, m}

path_id = 1;

delta = 0.9; % {m}
v_min = 0; % {ms-1}
v_max = 3; % {ms-1}
% omeg_max = 1.5; % {rads-1}
omeg_max = 3;
gamma = 3;

k1 = 5;
k2 = 5;

BF = @zbf2;

[x_t, u_t, r_t, h_t, errs] = bf_qp(BF, sim_time, step_size, Ts, x0, r0, path_id, p_o, delta, v_min, v_max, omeg_max, gamma, k1, k2);

plot_ctrls;
plot_res;
animate_pos;
