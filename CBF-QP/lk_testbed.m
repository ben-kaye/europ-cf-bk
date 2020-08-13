clear

%%% CONSTANTS AND MODEL PARAMETERS %%%
%--------------------------------------------------------------------------
g = 9.81; % {ms-2}
M = 1650; % {kg}
Cf = 1.33e5; % {N rad-1}
Cr = 9.88e4; % {N rad-1}
Iz = 2315.3; % {kg m2}
a_max = 0.3 * g; % {ms-2}
psc = 100;
v0 = 27.7; % {ms-1}
c = 10;
y_max = 0.9; % {m}
a = 1.11; % {m}
b = 1.59; % {m}
rd = 10; % {m}
gam = 1;

step_size = 1e-2; % {s}
run_time = 10; % {s}
N = floor(run_time/step_size);

state = zeros(4,1);
u_lat = 0;

%%% LQR SETUP %%%
%--------------------------------------------------------------------------

Kp = 5;
Kd = 0.4;

% xdot = At x + B u + E rd;
At = [ zeros(4,1), [ [ 1, v0, 0 ];...
        [ -1/M/v0 * (Cf + Cr), 0, 1/M/v0 * (b*Cr - a*Cf) - v0 ]; ...
        [ 0, 0, 1];...
        [ 1/Iz/v0 * ( b*Cr - a*Cf ), 0, -1/Iz/v0 * (a^2*Cf + b^2*Cr)] ] ];
B = [ 0; Cf/M; 0; a*Cf/Iz ];
E =[ 0; 0; -1; 0 ];
C = [ 1, 0, 20, 0];

Q = Kp*(C'*C) + Kd*C*(At'*At)*C';
R = 600;

K = lqr(At, B, Q, R);

%%% OSQP SETUP %%%
%--------------------------------------------------------------------------

P = diag([1, psc]);
q = [0;0];

Afc = [ 1 0 ];
[lfc, ufc] = getulfc(state(2),state(4),a,b,M,Cf,Cr,v0,rd,a_max);

Acbf_lk = getAcbf_lk(state(1),state(2),y_max,a_max,M,Cf);
ucbf_lk = getucbf_lk(state(1),state(2),state(4),y_max,a_max,M,Cf,Cr,a,b,v0,rd,gam);
lcbf_lk = -inf;

Ak = [ 1, -1 ];
%I think should be [0; 0; rd; 0];
uk = -K*(state - [0; 0; 0; rd]);
lk = uk;

A = [Afc; Acbf_lk; Ak];
u = [ufc; ucbf_lk; uk];
l = [lfc; lcbf_lk; lk];

solver = osqp;
solver.setup(P,q,A,l,u, 'warm_start',true,'verbose',false)


%%% SIMULATION - FORWARD EULER %%%
%--------------------------------------------------------------------------
err_count = 0;
for e = 1:N
    res = solver.solve();
    if( res.info.status_val ~= 1 )
        %%% SOLVER FAILED %%%
        err_count = err_count + 1;
    else
        %%% SOLVED %%%
        w = res.x;
        u_lat = w(1);
        delta = w(2);
    end
    
    %%% IMPLEMENT STATE TRANSITION VIA EULER FORWARD METHOD %%%
    state_dot = At*state + B*u_lat + E*rd;
    state = state + step_size * state_dot;
    
    %%% UPDATE CONSTRAINTS %%%
    [lfc, ufc] = getulfc(state(2),state(4),a,b,M,Cf,Cr,v0,rd,a_max);

    Acbf_lk = getAcbf_lk(state(1),state(2),y_max,a_max,M,Cf);
    ucbf_lk = getucbf_lk(state(1),state(2),state(4),y_max,a_max,M,Cf,Cr,a,b,v0,rd,gam);

    uk = -K*(state - [0; 0; 0; rd]);
    lk = uk;

    
    A = [Afc; Acbf_lk; Ak];
    u = [ufc; ucbf_lk; uk];
    l = [lfc; lcbf_lk; lk];
    
    solver.update('Ax',nonzeros(A),'u',u,'l',l);
end


%%% PLOTTING %%%
%--------------------------------------------------------------------------

%%% FUNCTIONS %%%
%--------------------------------------------------------------------------
function [lfc, ufc] = getulfc(v,r,a,b,M,Cf,Cr,v0,rd,a_max)
    o = Cf/v0 * (v + a*r) + Cr/v0 * (v -b*r) + M*v0*rd;
    lfc = -M*a_max + o;
    ufc = M*a_max + o;
end

function Acbf_lk = getAcbf_lk(y, v, y_max, a_max, M, Cf)
    ydot = v; % AGAIN CHECK THIS
    hF = y_max - y * sign(ydot) - 1/2/a_max*ydot^2;
    LgBf = Cf/( a_max*M*hF * (1+hF) );
    Acbf_lk = [ LgBf, 0 ];
end

function ucbf_lk = getucbf_lk(y,v,r,y_max,a_max,M,Cf,Cr,a,b,v0,rd,gam)
    ydot = v; % TEMPORARY MAY BE WRONG?? as ydot = v + v0*yaw BUT v is lat vel and y is lat displ
    Q = ydot * sign(ydot) - 1/M/a_max/v0 * ( Cf*(v + a*r) + Cr*(v + b*r) ) -v0*rd/a_max;
    
    hF = y_max - y * sign(ydot) - 1/2/a_max*ydot^2;
    
    Bf = -log(hF/(1+hF));
    LfBf = Q/hF/(hF+1);
    ucbf_lk = gam/Bf - LfBf;    
end