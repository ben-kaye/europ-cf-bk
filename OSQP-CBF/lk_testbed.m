% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                              Uses OSQP                              *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                 Implementing paper by A. Ames et al                 *
% *            https://ieeexplore.ieee.org/document/7782377             *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

%The LQR formulation is definitely not right. Works but barely.

clear

%%% CONSTANTS AND MODEL PARAMETERS %%%
%--------------------------------------------------------------------------
g = 9.81; % {ms-2}
M = 1650; % {kg}
Cf = 1.33e5; % {N rad-1}
Cr = 9.88e4; % {N rad-1}
Iz = 2315.3; % {kg m2}
y_max = 0.9; % {m}
a = 1.11; % {m}
b = 1.59; % {m}
R = 100; % {m}
a_max = 0.3 * g; % {ms-2}
v0 = 10; % {ms-1}
rd = v0/R; % {rad s-1}
c = 10;
gam = 1;
psc = 10; % // 100 recommended
alph = 3;

step_size = 5e-4; % {s}
sim_time = 25; % {s}
N = floor(sim_time/step_size);

state = zeros(4,1);
state = [ 0.1; 0.2; 0.1; 0.1 ];
u_lat = 0;

CONTROLLER_TYPE = 0; % {1: LYAPUNOV, 0: LQR}

%%% LQR SETUP %%%
%--------------------------------------------------------------------------

% QP program min xT P x + qT x st Ax <= b
% x in Re2
% x(1) = u
% x(2) = delta

Kp = 5;
Kd = 0.4;

% xdot = At x + B u + E rd;
At = [ zeros(4,1), [ [ 1, v0, 0 ];...
        [ -1/M/v0 * (Cf + Cr), 0, 1/M/v0 * (b*Cr - a*Cf) - v0 ]; ...
        [ 0, 0, 1];...
        [ 1/Iz/v0 * ( b*Cr - a*Cf ), 0, -1/Iz/v0 * (a^2*Cf + b^2*Cr)] ] ];
    
B = [ 0; Cf/M; 0; a*Cf/Iz ];

% E =[ 0; 0; -1; 0 ];
E = [ 0; 0; -1; 0 ];


C = [ 1, 0, 20, 0];

% paper expects Kd*C'At'*At*C but dims do not comply
Q = Kp*(C'*C) + Kd*(At'*(C'*C)*At);
R = 600;



K = lqr(At, B, Q, R);

%%% OSQP SETUP %%%
%--------------------------------------------------------------------------

P = diag([1, psc]);
q = [0;0];

Afc = [ 1 0 ];
[lfc, ufc] = getulfc(state(2),state(4),a,b,M,Cf,Cr,v0,rd,a_max);

[ Acbf_lk, ucbf_lk, hF ] = getCBFconstr(state(1),state(2),state(3),state(4),y_max,a_max,M,Cf,Cr,a,b,v0,rd,gam);

lcbf_lk = -inf;

Ak = [ 1, -1 ];
uk = -K*(state - [0;0;0;1]*rd);
lk = uk;

% y,v,yaw,r,rd,m,a,b,v0,C_r,C_f
[Aclf, uclf] = getCLFconstr(state(1),state(2),state(3),state(4),rd,M,a,b,v0,Cr,Cf, alph);
lclf = -inf;

if CONTROLLER_TYPE
    %%% CLF %%%
    A = [Afc; Acbf_lk; Aclf];
    u = [ufc; ucbf_lk; uclf];
    l = [lfc; lcbf_lk; uclf];
else
    %%% LQR %%%
    A = [Afc; Acbf_lk; Ak];
    u = [ufc; ucbf_lk; uk];
    l = [lfc; lcbf_lk; lk];
end


solver = osqp;
solver.setup(P,q,A,l,u, 'warm_start',true,'verbose',false)


%%% SIMULATION - FORWARD EULER %%%
%--------------------------------------------------------------------------
err_count = 0;
state_hist = zeros(9,N);
t_thresh1 = 10;
switched1 = false;
t_thresh2 = 15;
switched2 = false;

if (hF < 0)
    error('Outside of safety bound at sim start')
end

for e = 1:N
    
    
    res = solver.solve();
    if( res.info.status_val ~= 1 )
        %%% SOLVER FAILED %%%
        err_count = err_count + 1;
        break
    else
        %%% SOLVED %%%
        w = res.x;
        u_lat = w(1);
        delta = w(2);
    end
    
    %%% IMPLEMENT STATE TRANSITION VIA EULER FORWARD METHOD %%%
    state_dot = At*state + B*u_lat + E*rd;
    state = state + step_size * state_dot;
    
    state_hist(1:4,e) = state;
    state_hist(5, e) = state_dot(2);
    state_hist(6,e) = u_lat;
    state_hist(7,e) = delta;
    state_hist(8,e) = hF;
    state_hist(9,e) = rd;
    
    if e*step_size > t_thresh1 && ~switched1
        rd = -rd;
        switched1 = true;
    end
    
    if e*step_size > t_thresh2 && ~switched2
        R = 50;
        rd = +v0/R;
        switched2 = true;
    end
    
    
    
    %%% UPDATE CONSTRAINTS %%%
    % (v,r,a,b,M,Cf,Cr,v0,rd,a_max)
    [lfc, ufc] = getulfc(state(2),state(4),a,b,M,Cf,Cr,v0,rd,a_max);

    
    % (y,v,yaw,r,y_max,a_max,M,C_f,C_r,a,b,v_0,r_d,gam)
    [ Acbf_lk, ucbf_lk, hF ] = getCBFconstr(state(1),state(2),state(3),state(4),y_max,a_max,M,Cf,Cr,a,b,v0,rd,gam);

    
    uk = -K*(state - [0; 0; 0; rd]);
    lk = uk;

    [Aclf, uclf] = getCLFconstr(state(1),state(2),state(3),state(4),rd,M,a,b,v0,Cr,Cf,alph);
    lclf = -inf;

    
    
    if CONTROLLER_TYPE
        %%% CLF %%%
        A = [Afc; Acbf_lk; Aclf];
        u = [ufc; ucbf_lk; uclf];
        l = [lfc; lcbf_lk; uclf];
    else
        %%% LQR %%%
        A = [Afc; Acbf_lk; Ak];
        u = [ufc; ucbf_lk; uk];
        l = [lfc; lcbf_lk; lk];
    end


    
    solver.update('Ax',nonzeros(A),'u',u,'l',l);
end


%%% PLOTTING %%%
%--------------------------------------------------------------------------
t = 0:step_size:sim_time-step_size;
e=e-1;
t = t(1:e);

y = state_hist(1,1:e);
vel = state_hist(2,1:e);
yaw = state_hist(3,1:e);
yaw_rate = state_hist(4,1:e);
u_in = state_hist(6,1:e);
accel = state_hist(5,1:e);
dels = state_hist(7,1:e);
hFt = state_hist(8,1:e);
r_dt = state_hist(9,1:e);

figure(2);

subplot(3,1,1)

plot(t, y,'DisplayName','y', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
yline(-y_max, '-.', 'DisplayName', 'y_{min}','LineWidth',1,'Color', [1 0 0])
yline(y_max, '-.', 'DisplayName', 'y_{max}','LineWidth',1,'Color', [1 0 0])
title('Lateral Displacement (m)')

subplot(3,1,2)

plot(t,accel,'DisplayName','$\ddot{y}$', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252])
yline(-a_max, '-.', 'DisplayName', '$\ddot{y}_{min}$','LineWidth',1,'Color', [1 0 0])
yline(a_max, '-.', 'DisplayName', '$\ddot{y}_{max}$','LineWidth',1,'Color', [1 0 0])
% axis([0 t(end) -a_max-0.5 a_max+0.5]);

title('Acceleration (ms^{-1})')

L=legend();
set(L,'Interpreter','latex')

subplot(3,1,3)

plot(t, dels,'DisplayName','\delta','LineWidth', 1.5, 'Color', 1/255*[0, 157, 252])
hold on 
plot(t, hFt, 'DisplayName', 'h_F', 'LineWidth', 1.5, 'Color', [1 0 0])
hold off
legend();
title('Relaxation parameter')
%chatters if sample time is too high T>=O(1e-2)

figure(1)
subplot(3,1,1)
maxu = max(abs(u_in));
plot(t,u_in,'DisplayName','u', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252])
axis([ 0, t(e), -maxu, maxu]);
title('Wheel input');
subplot(3,1,2)
a_real = v0*(yaw_rate+r_dt);
maxa = max(abs(a_real));
plot(t, a_real, 'DisplayName','a', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252])
axis([ 0, t(e), -maxa, maxa]);
title('Centripetal Acceleration');
subplot(3,1,3)
plot(t, r_dt, 'DisplayName','r_d', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252])
title('Road Curvature Rate');


%%% FUNCTIONS %%%
%--------------------------------------------------------------------------
function [lfc, ufc] = getulfc(v,r,a,b,M,Cf,Cr,v0,rd,a_max)
    o = 1/Cf*( Cf/v0*( v + a*r ) + Cr/v0*( v - b*r ) + M*v0*rd );
    lfc = -M*a_max/Cf + o;
    ufc = M*a_max/Cf + o;
end

function [ Acbf_lk, ucbf_lk, h_F ] = getCBFconstr(y,v,yaw,r,y_max,a_max,M,C_f,C_r,a,b,v_0,r_d,gam)
%One shot implementation of the CBF constraint

    y_dot = v + v_0 * yaw;
    
    h_F = y_max - y * sign(y_dot) - y_dot^2/2/a_max;
    
    Q = 2*a_max*sign(y_dot) + C_r/M/v_0 * (b*r - v) - C_f/M/v_0 * (a*r + v) - v_0*r_d;
    
    LfBF = y_dot / 2 / a_max / M / h_F / (1+h_F) * Q;
    LgBF = y_dot * C_f / 2 / a_max / M / h_F / (1+h_F);
    
    BF = -log( h_F/(1+h_F) );
    
    Acbf_lk = [ LgBF, 0 ];
    ucbf_lk = -LfBF + gam/BF;
end

function [Aclf, uclf] = getCLFconstr(y,v,yaw,r,rd,m,a,b,v0,C_r,C_f, alpha)
    ydot = v + yaw*v0;
    V = y^2 + ydot^2;
    
    LgV = 2*ydot*C_f/m;
    LfV = 2*ydot*(y + C_r/m/v0*(b*r - v) - C_f/m/v0*(a*r + v) -rd*v0);
    
    Aclf = [LgV, -1];
    uclf = -alpha * V - LfV;
    
end