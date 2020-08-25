clear 


sim_time = 20; % {s}
step_size = 1e-3; % {s}
delta = 0.2; % {m} avoidance distance
gamma = 1;
v = 0.5; % {ms-1}
max_u = 1.5;
max_v = 4;
relax_weight = 100;

p_o = [ 5; 3 ];

phi_0 = -1; % {rad}
p_0 = [ 3; 2 ]; % {m; m}


x = [p_0; phi_0]; % {m; m; rad} state

p_r = [ 2; 3 ];
x_r = p_r;
phi_r = 1;
phi_rdot = 1/2;
v_r = 1;
v_rdot = 0;

% u = ones(2,1); % {ms-1; rads-1} control vec
[ u1, u2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
u = [ u1; u2 ];


%%% QP SETUP %%%
solver = osqp;

P = diag([0, 0, relax_weight]);
q = zeros(3,1);

[ Acbf, ucbf, h ] = getCBFconstraints(x([1,2]),x(3),p_o,u(1),max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 
[ Aclf, uclf ] = getCLFconstraints(x([1,2]),x_r,x(3),phi_r,phi_rdot,v_r,v_rdot); %(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
lclf = -inf;
lcbf = -inf;

Av = [ 1, 0, 0 ];
uv = max_v;
lv = 0;

Au = [ 0, 1, 0 ];
uu = max_u;
lu = -max_u;

l = [ lclf; lv; lu ];
A = [ Aclf; Av; Au ];
u = [ uclf; uv; uu ];

l = [ lcbf; lclf; lv; lu ];
A = [ Acbf; Aclf; Av; Au ];
u = [ ucbf; uclf; uv; uu ];

solver.setup(P,q,A,l,u,'warm_start',true,'verbose',false);

%%% SIMULATION %%%

N = floor(sim_time/step_size);

x_t = zeros(4,N);
x_rt = zeros(2,N);
h_t = zeros(1,N);

errc=0;

switched = false;
switched2 = false;
for e = 1:N
    % solve
    res = solver.solve();
    if res.info.status_val ~= 1
        %%% ERROR %%%
        errc = errc+1;
    else
        resx = res.x;
        u_act = resx([1,2]);
    end
    
    
    %simulate 
    x = simulate_state(step_size, x, u_act);
    [ x_r, phi_r, v_r ] = simulate_reference(step_size, x_r, phi_r, phi_rdot, v_r, v_rdot);

    %store
    x_t(1:3, e) = x;
    x_t([4,5], e) = u_act;
    x_rt(:, e) = x_r;


    %update
    % [Lf, Lg, h] = getCBF(x([1,2]), x(3), p_o, v, max_u, delta, gamma);


    % (p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)
    [ u1, u2 ] = controlFromCLF(x([1,2]), x_r, x(3), phi_r, phi_rdot, v_r, v_rdot);
    u_act = [ u1; u2 ];
    
%     truncate u
    u_max = [ max_v; max_u ];
    u_act = sign(u_act) .* min(abs(u_act), u_max);
    u_act(1) = max(u_act(1),0);
    
    % u = -Lf/Lg;
    
    
    [ Acbf, ucbf, h ] = getCBFconstraints(x([1,2]),x(3),p_o,u_act(1),max_u,delta,gamma); %(p_xy, phi, p_o, v, max_u, delta, gamma) 
    [ Aclf, uclf ] = getCLFconstraints(x([1,2]),x_r,x(3),phi_r,phi_rdot,v_r,v_rdot); %(p_xy, p_r, phi, phi_r, phi_rdot, v_r, v_rdot)

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
    
    l = [ lclf; lv; lu ];
    A = [ Aclf; Av; Au ];
    u = [ uclf; uv; uu ];

    l = [ lcbf; lclf; lv; lu ];
    A = [ Acbf; Aclf; Av; Au ];
    u = [ ucbf; uclf; uv; uu ];
    
    solver.update('Ax',A,'l',l,'u',u);

    
end



t = 0:step_size:sim_time;
t = t(1:e);
h_t = h_t(1:e);
x_t = x_t(:,1:e);
x_rt = x_rt(:,1:e);

u_t = x_t([4,5],:);

u_t(:,1) = [nan;nan];
figure(1)
plot(t,u_t,'-','LineWidth', 1.5)
hold on 
plot(t,h_t)
legend('v','\omega')

figure(2)
px = x_t(1,:);
py = x_t(2,:);
plot(px,py,'-','LineWidth',1.5, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
hold on
plot(x_rt(1,:), x_rt(2,:), '-.', 'LineWidth',1, 'Color',[1, 0, 0],'DisplayName', 'Reference');
viscircles(p_o', delta);
axis equal

function x_k_1 = simulate_state(time_step, x_k, u)
    phi = x_k(3);
    
    % truncate u if necessary
    % 


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
    p_dot = v * [ -sin(phi); cos(phi) ];

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
    
    %Trying something new here:
    
    
    %Test here
    Acbf = [ LfBF/v, LgBF, 0 ];
    ucbf = gamma/BF;

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