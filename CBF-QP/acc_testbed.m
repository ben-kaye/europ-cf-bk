% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                              Uses OSQP                              *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                 Implementing paper by A. Ames et al                 *
% *            https://ieeexplore.ieee.org/document/7040372             *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear

%%% CONSTANTS AND MODEL PARAMETERS %%%
%--------------------------------------------------------------------------
g = 9.81; % {ms-2}
m = 1650; % {kg}
f0 = 0.1; % {N}
f1 = 5; % {Ns/m}
f2 = 0.25; % {Ns+2/m}
vd = 22; % {ms-1} desired velocity
eps = 5; % model param // suggested 10
gam = 1; % model param // suggested 1
ca = 0.25; % accel lim l
cd = 0.25; % accel lim u

psc = exp(-5); % TO BE OVERRIDDEN
psc_no_acc_lim = 5e1; % big number needed to keep delta down
psc_acc_lim = 500;

T = 1e-3; % {s} sample period
sim_step = 1e-3; % {s} simulation step size
sim_time = 20; % {s}
N = sim_time/T; % sim length
z0 = 150; % {m} initial distance
vstart = 18; % {m} initial velocity)
% vlstart = 10; % {ms-1} initial speed
v0 = 10; % {ms-1} initial lead car velocity
w0 = [ vstart; v0;  z0 ]; % {m} {ms-1} {m}
w = w0; % state :: w(1) = v; w(2) = vl; w(3) = z (dist between lead and car)
mu = 0;
a_L = 0; % {ms-2} lead car accel


FCBF_ON = 0; %{1:ON 0:OFF} Whether to include hard constraint on accel, and appropriate CBF

if FCBF_ON
    psc = psc_acc_lim;
else 
    psc = psc_no_acc_lim;
end

%%% STATE SPACE %%%
% dx = At x + B u + D(x)

Fr = getFr(w(1), f0, f1, f2);

At = [ zeros(2,3); [-1, 1, 0] ];
B = [ 1/m; 0; 0 ];
D = [ -Fr/m; a_L; 0 ];


%%% FORMULATE INITIAL OSQP MATRICES %%%
%--------------------------------------------------------------------------


P = 2 * diag( [ 1/m^2, psc ]);
q = -2/m * [ Fr; 0 ]; % needs to be updated on each sim

lcc = -cd * m * g;
ucc = ca * m * g;
Acc = [ 1 0 ];

lclf = -inf;
[Aclf, uclf] = getCLFconstr(m, w(1), Fr, vd, eps );

lcbf = -inf;
[Acbf, ucbf] = getCBFconstr(m, w(1), w(2), w(3), Fr, gam);

lfcbf = -inf;
[Afcbf, ufcbf] = getFCBFconstr(m, w(1), w(2), w(3), Fr, cd, g);

if FCBF_ON
    l = [ lclf; lcbf; lcc; lfcbf ];
    u = [ uclf; ucbf; ucc; ufcbf ];
    A = [ Aclf; Acbf; Acc; Afcbf ];
else
    l = [ lclf; lcbf ];
    u = [ uclf; ucbf ];
    A = [ Aclf; Acbf ];
end

solver = osqp;
solver.setup(P, q, A, l, u, 'warm_start', true, 'verbose', false);

%%% SIMULATION %%%
%--------------------------------------------------------------------------

pastW = zeros(6, N*floor(T/sim_step));
p = 0;
errCount = 0;
u_in = 0;
del = 0;

a_L = 3;

for e = 1:N

%     if (w(2) > 1.1*vd)
%         fprintf('Desired velocity exceeded, terminating\n')
%         break
%     end
    
    res = solver.solve();
    if ~strcmp(res.info.status, 'solved')
%         fprintf('ERROR\n')
        errCount = errCount + 1;
%         break
    else
        x_opt = res.x;
        u_in = x_opt(1); % wheel force
        del = x_opt(2); % relaxation param
    end 
    
    %%% forwards euler integration
    
%     for k = 1:floor(T/sim_step)
    p = p + 1;


    pastW(1,p) = w(1); % v
    pastW(3,p) = w(2); % v_L
    pastW(4,p) = w(3); % z
    
    pastW(5,p) = del;

    Fr = getFr(w(1), f0, f1, f2);
    
    D = [ -Fr/m; a_L; 0 ];
    wdot =  At*w + B*u_in + D;
    w = w + wdot * sim_step;
    
    pastW(2, p) = wdot(1); %vdot
    pastW(6, p) = res.info.obj_val;
%     end
    

    %%% UPDATE CONSTRAINTS %%%
    
    Fr = getFr(w(1), f0, f1, f2);
    
    % update q
    q = -2 * Fr / m * [ 1; 0 ];
    
    [Aclf, uclf] = getCLFconstr(m, w(1), Fr, vd, eps);
    [Acbf, ucbf] = getCBFconstr(m, w(1), w(2), w(3), Fr, gam);
    [Afcbf, ufcbf] = getFCBFconstr(m, w(1), w(2), w(3), Fr, cd, g);
    
    if FCBF_ON
        l = [ lclf; lcbf; lcc; lfcbf ];
        u = [ uclf; ucbf; ucc; ufcbf ];
        A = [ Aclf; Acbf; Acc; Afcbf ];
    
    else
        l = [ lclf; lcbf ];
        u = [ uclf; ucbf ];
        A = [ Aclf; Acbf ];
    end
    
    solver.update('Ax', nonzeros(A), 'l', l, 'u', u, 'q', q);    
end

%%% PLOT RESULTS %%%
%--------------------------------------------------------------------------
t = 0:sim_step:N*T;
t = t(1:p);
v = pastW(1,1:p);
a = pastW(2,1:p);
v_L = pastW(3,1:p);
z = pastW(4,1:p);
del = pastW(5,1:p);
obj = pastW(6,1:p);
figure(1)

subplot(2,2,1)
plot(t,v, 'DisplayName','v', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
hold on
plot(t,z, 'DisplayName','z', 'LineWidth',1.5,'Color', 1/255*[252, 169, 0]);
plot(t,v_L, 'DisplayName','V_L', 'LineWidth',1.5,'Color', [ 1 0 0]);
yline(vd,'k-.','DisplayName','v_d');
% axis([0, 0.99*p*T, min(v)-1, max(z)+1])
xlabel('t (s)')
ylabel('(ms^{-1}) or (m)');
title('Velocity, Separation')
legend()
hold off

subplot(2,2,2)
plot(t,z-1.8*v, 'DisplayName','Safe Heading (Hard)', 'LineWidth',1.5,'Color', 1/255*[252, 169, 0]);
hold on
plot(t, vd-v, 'DisplayName', 'Soft Constraint', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
yline(0, '-.', 'DisplayName', 'Constraint Violated', 'LineWidth',1,'Color', [1 0 0]);
hold off
% axis([0, 0.99*p*T, min(vd-v)-1, max(z-1.8*v)+1])
title('Constraints')
xlabel('t (s)')
legend()

subplot(2,2,3)
plot(t, del, 'DisplayName','\delta','LineWidth',1.5,'Color', 1/255*[0, 157, 252])
title('Relaxation Param')
xlabel('t (s)')
legend()
mindel = min(del);
maxdel = max(del);
% axis([0, 0.99*p*T, min(del)-1, max(del)+1])
sgtitle(sprintf('Relaxation weight %f\nEpsilon %f\nGamma %f',psc,eps,gam));

subplot(2,2,4)
% plot(t, obj, 'DisplayName', 'J');
% title('Objective function')

plot(t,a, 'DisplayName','a','LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
yline(-cd*g, '-.', 'DisplayName', 'a_{min}','LineWidth',1,'Color', [1 0 0])
yline(ca*g, '-.', 'DisplayName', 'a_{max}','LineWidth',1,'Color', [1 0 0])
% axis([0, 0.99*p*T, -3.3, 3.3])
title('Acceleration')
xlabel('t (s)')
ylabel('a (ms^{-2})')
legend()


%%% FUNCTIONS %%% TAKE CARE ON SIGNS
%--------------------------------------------------------------------------
function Fr = getFr(v, f0, f1, f2)
    Fr = f0 + f1 * v + f2 * v^2;
end

function [Aclf, uclf] = getCLFconstr(m, v, Fr, vd, eps)
    y = (v - vd);

    LfV = -2*y*Fr/m;
    LgV = 2*y/m;
    
    V = y^2;
    
    Aclf = [ LgV, -1 ];
    uclf = -eps*V - LfV;
end

function [Afcbf, ufcbf] = getFCBFconstr(m, v, vL, z, Fr, cd, g)
    
    zdot = vL - v;
    
    h = z - 1.8*v;
    
    hF = h - 1/2/cd/g*zdot^2;
    
    LgBF = (zdot - 1.8*cd*g)/cd/g/m/hF^2;
    
    LfBF = -( m*cd*g*zdot + Fr*(zdot - 1.8*cd*g) )/cd/g/m/hF^2;
    
    Afcbf = [ LgBF, 0 ];
    
    ufcbf = hF - LfBF;
end


function [Acbf, ucbf] = getCBFconstr(m, v, vL, z, Fr, gam)
    h = z - 1.8*v;
    
    LgB = 1.8/m/h/(1 + h);   
    LfB = -1/m/h/(1 + h)*( 1.8*Fr + m*(vL - v) );
    
    B = -log( h/(1 + h) );
    
    ucbf = gam/B - LfB ;
    
    Acbf = [ LgB, 0 ];
end