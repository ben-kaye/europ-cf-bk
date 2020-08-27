% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                              Uses OSQP                              *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear

%%% SIMULATION PARAMETERS %%%

sim_time = 10; % {s}
step_size = 1e-4; % {s}
Ts = 1e-3; % {s} f = 100Hz
p_o = [ -1; 6]; % {m; m} object position
p_0 = [ 3; 2 ]; % {m; m} initial drone position
phi_0 = -1; % {rad} initial drone bearing
x = [p_0; phi_0]; % {m; m; rad} state

%%% CBF-CLF-QP PARAMETERS %%%

delta = 0.4; % {m} avoidance distance
v = 0.5; % {ms-1}
max_u = 1.5; % {rads-1}
max_v = 12; % {ms-1}
min_v = 0.1; % {ms-1}
relax_weight = 10; % relax penalty
R = diag([ 1, 1 ]); % control penalty
gamma = 100;
K1 = 5;
K2 = 5;

%%% REFERENCE SIGNAL GENERATION %%%

p_r = [ 2; 3 ];
x_r = p_r;
phi_r = 1;
phi_rdot = 1/2;
v_r = 1;
v_rdot = 0;

% compute CLF input to track ref
[ ctrl_lf1, ctrl_lf2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot,K1,K2);
ctrl = [ ctrl_lf1; ctrl_lf2 ];
ctrl = constrain_u(ctrl, min_v, max_v, max_u); % {ms-1; rads-1} control vector

%%% QP SETUP %%%
solver = osqp;

[ Abf, ubf, h ] = getBFconstraint(x([1,2]),x(3),p_o,ctrl(1),max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 
lbf = -inf;

% P = blkdiag(R,relax_weight);
P = R;
% q = [ -2*R'*ctrl; 0 ];
q = -2*R'*ctrl;

Av = [ 1, 0 ];
uv = max_v;
lv = min_v; % might rather, do 0

Au = [ 0, 1 ];
uu = max_u;
lu = -max_u;

% l = [ lclf; lv; lu ];
% A = [ Aclf; Av; Au ];
% u = [ uclf; uv; uu ];

l = [ lbf; lv; lu ];
A = [ Abf; Av; Au ];
u = [ ubf; uv; uu ];

solver.setup(P,q,Abf,lbf,ubf,'warm_start',true,'verbose',false);

%%% SIMULATION %%%

N = floor(sim_time/step_size);
Ns = floor(sim_time/Ts);

if (Ns>N)
    error('Sample_Time < Simulation_Step');
end

x_t = zeros(4,N);
x_rt = zeros(2,N);
h_t = zeros(1,N);

errc=0;

switched = false;
switched2 = false;
for e = 1:Ns
    % solve
    res = solver.solve();
    if res.info.status_val ~= 1
        %%% ERROR %%%
        errc = errc+1;
        error('Primal infeasibility')
    else
        resx = res.x;
        ctrl = resx([1,2]);
    end
    
    
    %simulate substep
    for f = 1:floor(N/Ns)
        x = simulate_state(step_size, x, ctrl);
        [ x_r, phi_r, v_r ] = simulate_reference(step_size, x_r, phi_r, phi_rdot, v_r, v_rdot);
    end
    
    %store
    x_t(1:3, e) = x;
    x_t([4,5], e) = ctrl;
    x_rt(:, e) = x_r;


    %update
    % [Lf, Lg, h] = getCBF(x([1,2]), x(3), p_o, v, max_u, delta, gamma);


    v_old = ctrl(1);

    % (p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    [ ctrl_lf1, ctrl_lf2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot, K1,K2);
    ctrl = [ ctrl_lf1; ctrl_lf2 ];    
    ctrl = constrain_u(ctrl, min_v, max_v, max_u);
    
    %%% USING v_old instead of ctrl_lf1
    [ Abf, ubf, h ] = getBFconstraint(x([1,2]),x(3),p_o,v_old,max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 

    if h < 0 
        % error('Safety violated')
    end
    
    h_t(e) = h;
    if e*step_size > 10  && ~switched
        phi_rdot = -phi_rdot;
        v_rdot = -0.01;
        switched = true;
        
    end
    
    if e*step_size > 15 && ~switched2
       v_rdot = 0;
       phi_rdot = 0;
       switched2 = true;
    end
    

%     l = [ lbf; lv; lu ];
%     A = [ Abf; Av; Au ];
%     u = [ ubf; uv; uu ];



    % q = [ -2*R'*ctrl; 0 ];
    q = -2*R'*ctrl;
    
    solver.update('Ax',Abf,'l',lbf,'u',ubf,'q',q);

    
end

%%% PLOTTING %%%

t = 0:step_size:sim_time;
z = e*N/Ns;

t = t(1:e);
h_t = h_t(1:e);
x_t = x_t(:,1:e);
x_rt = x_rt(:,1:e);

u_t = x_t([4,5],:);

plot_path

%%% FUNCTIONS %%%

function u = constrain_u(u,min_v,max_v,max_omega)
    max_vals = [ max_v; max_omega ];
    u = sign(u) .* min(abs(u), max_vals);
    u(1) = max(u(1),min_v);    
end

function x_k_1 = simulate_state(time_step, x_k, u)
    phi = x_k(3);

    x_dot = [ cos(phi), 0; sin(phi), 0; 0, 1 ]*u;
    
    x_k_1 = x_k + x_dot*time_step;
    
end

function [x_rk1, phi_rk1, v_rk1] = simulate_reference(time_step, x_rk, phi_rk, phi_rdot, v_rk, v_rdot)
    z = simulate_state(time_step, [x_rk; phi_rk], [v_rk; phi_rdot]);

    x_rk1 = z([1,2]);
    phi_rk1 = z(3);

    v_rk1 = v_rk + time_step*v_rdot;
    v_rk1 = max(v_rk1, 0);
end

function [ Abf, ubf, h ] = getBFconstraint(p_xy, phi, p_o, v, max_turn, delta, gamma)
    p_dot = v * [ cos(phi); sin(phi) ];

    p_xo = p_o - p_xy;
    
    sign_term = sign(p_dot' * [ 0, 1; -1, 0 ] * p_xo); %if innerprod > 0;


    % cross(p_dot, k)
    r = sign_term*[ 0, 1; -1, 0 ] * p_dot / max_turn;

    p_c = p_xy + r;
    
    z = p_o - p_c;
    
    h = z'*z - (v/max_turn + delta)^2;
    
    t = exp(1/h^2);
    sig = (t-1)/(t+1);

    
    Abf = [ 0, -1 ];
    ubf = sig*sign_term*max_turn;

end


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