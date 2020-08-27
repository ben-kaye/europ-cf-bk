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

%%% REFERENCE SIGNAL GENERATION %%%

p_r = [ 2; 3 ];
x_r = p_r;
phi_r = 1;
phi_rdot = 1/2;
v_r = 1;
v_rdot = 0;

% compute CLF input to track ref
[ ctrl1, ctrl2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
ctrl = [ ctrl1; ctrl2 ];
ctrl = constrain_u(ctrl, min_v, max_v, max_u); % {ms-1; rads-1} control vector

%%% QP SETUP %%%
solver = osqp;

[ Acbf, ucbf, h ] = getCBFconstraints(x([1,2]),x(3),p_o,ctrl(1),max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 
[ Aclf, uclf ] = getCLFconstraints(x([1,2]),x_r,x(3),phi_r,phi_rdot,v_r,v_rdot); %(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
lclf = -inf;
lcbf = -inf;

P = blkdiag(R,relax_weight);
q = [ -2*R'*ctrl; 0 ];

Av = [ 1, 0, 0 ];
uv = max_v;
lv = min_v; % might rather, do 0

Au = [ 0, 1, 0 ];
uu = max_u;
lu = -max_u;

% l = [ lclf; lv; lu ];
% A = [ Aclf; Av; Au ];
% u = [ uclf; uv; uu ];

l = [ lcbf; lclf; lv; lu ];
A = [ Acbf; Aclf; Av; Au ];
u = [ ucbf; uclf; uv; uu ];

solver.setup(P,q,A,l,u,'warm_start',true,'verbose',false);

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
    [ ctrl_lf1, ctrl_lf2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
    ctrl = [ ctrl_lf1; ctrl_lf2 ];    
    ctrl = constrain_u(ctrl, min_v, max_v, max_u);
    
    
    %%% USING v_old instead of ctrl_lf1
    [ Acbf, ucbf, h ] = getCBFconstraints(x([1,2]),x(3),p_o,v_old,max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 
    [ Aclf, uclf ] = getCLFconstraints(x([1,2]),x_r,x(3),phi_r,phi_rdot,v_r,v_rdot); %(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    
    if h < 0 
        error('Safety violated')
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
    
%     l = [ lclf; lv; lu ];
%     A = [ Aclf; Av; Au ];
%     u = [ uclf; uv; uu ];

    l = [ lcbf; lclf; lv; lu ];
    A = [ Acbf; Aclf; Av; Au ];
    u = [ ucbf; uclf; uv; uu ];
    

    q = [ -2*R'*ctrl; 0 ];

    
    solver.update('Ax',A,'l',l,'u',u,'q',q);

    
end

if(errc > 0)
    fprintf('Sim complete. %d errors encountered\n', errc);
else
    fprintf('Sim complete.');
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

function [ Acbf, ucbf, h ] = getCBFconstraints(p_xy, phi, p_o, v, max_u, delta, gamma)
    p_dot = v * [ cos(phi); sin(phi) ];

    p_xo = p_o - p_xy;
    
    sign_term = sign(p_dot' * [ 0, 1; -1, 0 ] * p_xo); %if innerprod > 0;


    % cross(p_dot, k)
    r = sign_term*[ 0, 1; -1, 0 ] * p_dot / max_u;

    p_c = p_xy + r;
    
    z = p_o - p_c;
    
    h = z'*z - (v/max_u + delta)^2;

    LfBF = 2*z'*p_dot/h/(1+h);
    
    LgBF = sign_term*LfBF/max_u;
    
    BF = -log(h/(1+h));
    
   
    
    
    % %Trying something new here:
    % Acbf = [ LfBF/v, LgBF, 0 ];
    % ucbf = gamma/BF;

    Acbf = [ 0, LgBF, 0 ];
    ucbf = gamma/BF - LfBF; 
end

function [ Aclf, uclf ] = getCLFconstraints(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p_xy);
    e2 = e12(2);
    e1 = e12(1);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;
    
    V = 1/2*e1^2 + 1/2*e2^2 + 2*(1-cos(e3_aux/2));

    LfV = v_r*(e1*cos(e3) - e2*sin(alpha)) +  sin(e3_aux/2)*(phi_rdot + 2*e2*v_r*cos(e3_aux/2 - alpha) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2));
    LgV = [ -e1, -sin(e3_aux/2)*(1 + e1*v_r/(v_r^2+e2^2)) ];
    
    Aclf = [ LgV, -1 ];
%     uclf = -LfV;
    uclf = -V - LfV;

end

function [ u1, u2 ] = controlFromCLF(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p_xy);

    e2 = e12(2);
    e1 = e12(1);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;

    k1 = 5;
    k2 = 5;

    u2 = phi_rdot + 2*e2*v_r*cos(e3_aux/2-alpha) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2) + k2*sin(e3_aux/2);
    u1 = v_r*cos(e3) - u2*v_r*sin(e3_aux/2)/(v_r^2 + e2^2) + k1*e1;


%     alpha = atan(v_r*e2);
%     e3_aux = alpha + e3;
%     
%     u2 = phi_rdot + 2*e2*v_r*cos(e3_aux/2-alpha) + (e2*v_r*sin(e3) + e2*v_rdot)/(1 + (v_r * e2)^2) + k2*sin(e3_aux/2);
%     u1 = v_r*cos(e3) - u2*v_r*sin(e3_aux/2)/(1 + (v_r * e2)^2) + k1*e1;
end