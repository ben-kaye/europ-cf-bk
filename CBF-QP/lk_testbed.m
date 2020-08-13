clear

%%% CONSTANTS AND MODEL PARAMETERS %%%
%--------------------------------------------------------------------------
g = 9.81;
M = 1650;
Cf = 1.33e5;
Cr = 9.88e4;
Iz = 2315.3;
a_max = 0.3 * g;
psc = 100;
v0 = 27.7;
c = 10;
ymax = 0.9;
a = 1.11;
b = 1.59;

%%% OSQP SETUP %%%
%--------------------------------------------------------------------------

P = diag(1, psc);
q = zeros(1,2);

Afc = [ 1 0 ];
[lfc, ufc] = getulfc();

Acbf_lk = getAcbf_lk();
ucbf_lk = getucbf_lk();
lcbf_lk = -inf;

A = [Afc; Acbf_lk];
u = [ufc; ucbf_lk];
l = [lfc; lcbf_lk];

solver = osqp;
solver.setup()


%%% SIMULATION - FORWARD EULER %%%
%--------------------------------------------------------------------------

% xdot = At x + B u + E rd;
At = [ zeros(4,1), [ [ 1, v0, 0 ];...
        [ -1/M/v0 * (Cf + Cr), 0, 1/M/v0 * (b*Cr - a*Cf) - v0 ]; ...
        [ 0, 0, 1];...
        [ 1/Iz/v0 * ( b*Cr - a*Cf ), 0, -1/Iz/v0 * (a^2*Cf + b^2*Cr)] ] ];
B = [ 0; Cf/M; 0; a*Cf/Iz ];
E =[ 0; 0; -1; 0 ];

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
    [lfc, ufc] = getulfc();

    Acbf_lk = getAcbf_lk();
    ucbf_lk = getucbf_lk();

    A = [Afc; Acbf_lk];
    u = [ufc; ucbf_lk];
    l = [lfc; lcbf_lk];
    
    solver.update('Ax',nonzeros(A),'u',u,'l',l);
end


%%% PLOTTING %%%
%--------------------------------------------------------------------------

%%% FUNCTIONS %%%
%--------------------------------------------------------------------------
function [lfc, ufc] = getulfc(v,r,a,b,M,Cf,Cr,v0,a_max)
    o = Cf/v0 * (v + a*r) + Cr/v0 * (v -b*r) + M*v0*rd;
    lfc = -M*a_max + o;
    ufc = M*a_max + o;
end

function Acbf_lk = getAcbf_lk()
    ydot = v; % AGAIN CHECK THIS
    hF = y_max - y * sign(ydot) - 1/2/a_max*ydot^2;
    LgBf = Cf/( a_max*M*hF * (1+hF) );
    Acbf_lk = [ LgBf, 0 ];
end

function ucbf_lk = getucbf_lk()
    ydot = v; % TEMPORARY MAY BE WRONG?? as ydot = v + v0*yaw BUT v is lat vel and y is lat displ
    Q = ydot * sign(ydot) - 1/M/a_max/v0 * ( Cf*(v + a*r) + Cr*(v + b*r) ) -v0*rd/a_max;
    
    hF = y_max - y * sign(ydot) - 1/2/a_max*ydot^2;
    
    Bf = -log(hF/(1+hF));
    LfBf = Q/hF/(hF+1);
    ucbf_lk = gam/Bf - LfBf;    
end