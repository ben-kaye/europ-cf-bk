
sim_time = 20; % {s}
step_size = 1e-3; % {s}
delta = 0.2; % {m} avoidance distance
gamma = 1;
v = 0.5; % {ms-1}
max_u = 0.3;

p_o = [ 3; 3 ];

phi_0 = -.3; % {rad}
p_0 = [ 1; 2 ]; % {m; m}

x = [p_0; phi_0]; % {m; m; rad} state

p_r = [ 2; 3 ];
phi_r = 1;
phi_rdot = 1/2;
v_r = 1;
v_rdot = 0;

u = zeros(2,1); % {ms-1; rads-1} control vec
[ u1, u2 ] = controlFromCLF(x([1,2]), p_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
u = [ u1; u2 ];


N = floor(sim_time/step_size);

x_t = zeros(4,N);
x_rt = zeros(2,N);
h_t = zeros(1,N);

for e = 1:N
    % solve

    %simulate 
    x = simulate_state(step_size, x, u);
    [ x_r, phi_r, v_r ] = simulate_reference(step_size, x_r, phi_r, phi_rdot, v_r, v_rdot);

    %store
    x_t(1:3, e) = x;
    x_t([4,5], e) = u;
    x_rt(:, e) = x_r;


    %update
    % [Lf, Lg, h] = getCBF(x([1,2]), x(3), p_o, v, max_u, delta, gamma);
%     h_t(e) = h;

    % (p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    [ u1, u2 ] = controlFromCLF(x([1,2]), p_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
    u = [ u1; u2 ];
    
    % u = -Lf/Lg;
end

figure(1)
t = 0:step_size:sim_time;
t = t(1:e);
h_t = h_t(1:e);
x_t = x_t(:,1:e);
plot(t,h_t)

figure(2)
px = x_t(1,:);
py = x_t(2,:);
plot(px,py,'.');
hold on
viscircles(p_o', delta);
axis equal

function x_k_1 = simulate_state(time_step, x_k, u)
    p_xy = x_k([1, 2]);
    phi = x_k(3);

    v = u(1);
    phi_dot = u(2);

    % truncate u if necessary
    % u = sign(u) * max(abs(u), u_max);

    p_xydot = u(1) * [ cos(phi); sin(phi) ];
    x_dot = [ p_xydot; u(2) ];
    
    x_k_1 = x_k + x_dot*time_step;
    
end

function [x_rk1, phi_rk1, v_rk1] = simulate_reference(time_step, x_rk, phi_rk, phi_rdot, v_rk, v_rdot)
    z = simulate_state(time_step, [x_rk; phi_rk], [v_rk; phi_rdot]);

    x_rk1 = z([1,2]);
    phi_rk1 = z(3);

    v_rk1 = v_rk + time_step*v_rdot;
end

function [ LfBF, LgBF, h ] = getCBF(p_xy, phi, p_o, v, max_u, delta, gamma)
    p_dot = v * [ -sin(phi); cos(phi) ];

    p_xo = p_o - p_xy;
    
    innerprod = p_dot' * [ 0, 1; -1, 0 ] * p_xo; %if innerprod > 0;


    % cross(p_dot, k)
    r = sign(innerprod)*[ 0, 1; -1, 0 ] * p_dot / max_u;

    p_c = p_xy + r;
    
    z = p_o - p_c;
    
    h = z'*z - delta^2;

    LfBF = 2*z'*p_dot/h/(1+h);
    LgBF = sign(innerprod)*LfBF/max_u;
    
    BF = -log(h/(1+h));
    
    Acbf = [ LgBF, 0 ];
    ucbf = gamma/BF - LfBF;    
end

function [ u1, u2 ] = controlFromCLF(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p_xy);

    e2 = e12(2);
    e1 = e12(1);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;

    k1 = 2;
    k2 = 2;

    u2 = phi_rdot + 2*e2*v_r*cos(e3/2) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2) + k2*sin(e3_aux/2);
    u1 = v_r*cos(e3) - u2*v_r*sin(e3_aux/2)/(v_r^2 + e2^2) + k1*e1;
end