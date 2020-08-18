
sim_time = 10; % {s}
step_size = 1e-3; % {s}
delta = 0.2; % {m} avoidance distance
gamma = 1;
v = 0.5; % {ms-1}
max_u = 0.3;

p_o = [ 3; 3 ];

phi_0 = -1.2; % {rad}
p_0 = [ 1; 2 ]; % {m; m}

x = [p_0; phi_0]; % {m; m; rad} state
u = 0;


N = floor(sim_time/step_size);

x_t = zeros(4,N);
h_t = zeros(1,N);

for e = 1:N
    % solve

    %simulate 
    x = simulate_state(step_size, x, u, v);

    %store
    x_t(1:3, e) = x;
    x_t(4, e) = u;

    %update
    [Lf, Lg, h] = getCBF(x([1,2]), x(3), p_o, v, max_u, delta, gamma);
    h_t(e) = h;
    
    u = -Lf/Lg;
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

function x_k_1 = simulate_state(time_step, x_k, u, v)
    p_xy = x_k([1, 2]);
    phi = x_k(3);

    % truncate u if necessary
    % u = sign(u) * max(abs(u), u_max);

    p_xydot = v * [ -sin(phi); cos(phi) ];
    x_dot = [ p_xydot; u ];
    
    x_k_1 = x_k + x_dot*time_step;
    
end

function [ LfBF, LgBF, h ] = getCBF(p_xy, phi, p_o, v, max_u, delta, gamma)
    p_dot = v * [ -sin(phi); cos(phi) ];
    
    % cross(p_dot, k)
    r = -[ 0, 1; -1, 0 ] * p_dot / max_u;

    p_c = p_xy - r;
    
    z = p_o - p_c;
    
    h = z'*z - delta^2;

    LfBF = 2*z'*p_dot/h/(1+h);
    LgBF = -LfBF;
    
    BF = -log(h/(1+h));
    
    Acbf = [ LgBF, 0 ];
    ucbf = gamma/BF - LfBF;    
end