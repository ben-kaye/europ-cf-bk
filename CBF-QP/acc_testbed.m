% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *        Implementing Automatic Cruise-Control with CBF-CLF-QP        *
% *                       using paper by A. Ames                        *
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
vd = 24; % {ms-1} ideal velocity
v0 = 13.89; % {ms-1} lead car velocity
eps = 10; % model param
gam = 1; % model param
ca = 0.3; % accel lim l
cd = 0.3; % accel lim u
psc = 1e-5; % relaxation weight

% psc = 1e-1;

T = 1e-2; % {s} sample period
sim_step = 1e-2; % {s} simulation step size
sim_time = 30; % {s}
N = sim_time/T; % sim length
z0 = 100; % {m} initial distance
xstart = 900;
vstart = 9; % {ms-1} initial speed
w0 = [ xstart; vstart;  z0 ]; % {m} {ms-1} 
w = w0; % state variable of car and lead
mu = 0;

%%% FORMULATE INITIAL OSQP MATRICES %%%
%--------------------------------------------------------------------------
Fr = getFr(w(2), f0, f1, f2);

P = 2 * diag( [ 1/m^2, psc ]);
q = -2/m * [ Fr; 0 ]; % needs to be updated on each sim
% 
lcc = -ca * m * g;
ucc = cd * m * g;
Acc = [ 1 0 ];

lclf = -inf;
uclf = getbclf(m, w(2), vd, eps, Fr);
Aclf = getAclf(m, w(2), vd);

lcbf = -inf;
ucbf = getbcbf(m, w(2), w(3), v0, Fr, gam);
Acbf = getAcbf(m, w(2), w(3));

lfcbf = -inf;
ufcbf = getbfcbf();
Afcbf = getAfcbf();

l = [ lclf; lcbf ];
u = [ uclf; ucbf ];
A = [ Aclf; Acbf ];

%%% using hard constraint or FCBF delete appropriately %%%
l = [ l; lcc ];
u = [ u; ucc ];
A = [ A; Acc ];
% 
% l = [ l; lfcbf ];
% u = [ u; ufcbf ];
% A = [ A; Afcbf ];
% %% %%%

solver = osqp;
solver.setup(P, q, A, l, u, 'warm_start', true, 'verbose', false);

%%% SIMULATION %%%
%--------------------------------------------------------------------------

pastW = zeros(5, N*floor(T/sim_step));
p = 0;
errCount = 0;
u_in = 0;

for e = 1:N
    
    if (w(2) > 1.1*vd)
        fprintf('Desired velocity exceeded, terminating\n')
        break
    end
    
    res = solver.solve();
    if ~strcmp(res.info.status, 'solved')
%         fprintf('ERROR\n')
        errCount = errCount + 1;
%         break
    else
        x = res.x;
        u_in = x(1); % control input relative
        del = x(2);
    end
    
    %%% forwards euler integration
    
    for k = 1:floor(T/sim_step)
        p = p + 1;
        pastW([1 2],p) = w(1:2);
        pastW(4,p) = w(3);
        pastW(5,p) = del;
    
        Fr = getFr(w(2), f0, f1, f2);
        
        wdot =  [ w(2); -Fr/m; v0 - w(2) ] + ...
                [ 0; 1/m; 0 ] * u_in;
            
        w = w + wdot * sim_step;
        pastW(3, p) = wdot(2);
    end
    
    Fr = getFr(w(2), f0, f1, f2);
    
    %%% UPDATING q
    q = -2 * Fr / m * [ 1; 0 ];
    
    
    lclf = -inf;
    uclf = getbclf(m, w(2), vd, eps, Fr);
    Aclf = getAclf(m, w(2), vd);
 
    lcbf = -inf;
    ucbf = getbcbf(m, w(2), w(3), v0, Fr, gam);
    Acbf = getAcbf(m, w(2), w(3));

    lfcbf = -inf;
    ufcbf = getbfcbf();
    Afcbf = getAfcbf();
    
    l = [ lclf; lcbf ];
    u = [ uclf; ucbf ];
    A = [ Aclf; Acbf ];

    %%% using hard constraint or FCBF delete appropriately %%%
    l = [ l; lcc ];
    u = [ u; ucc ];
    A = [ A; Acc ];

%     % l = [ l; lfcbf ];
%     % u = [ u; ufcbf ];
%     % A = [ A; Afcbf ];
    %%% %%%
    
    solver.update('Ax', nonzeros(A), 'l', l, 'u', u, 'q', q);    
end

%%% PLOT RESULTS %%%
%--------------------------------------------------------------------------
t = 0:sim_step:N*T;
x = [pastW(1,:), NaN];
v = [pastW(2, :), NaN];
a = [pastW(3, :), NaN];
z = [pastW(4, :), NaN];
del = [pastW(5,:), NaN];
figure(1)

subplot(2,2,1)
% plot(t,x, 'DisplayName','x');
hold on
plot(t,v, 'DisplayName','v');
plot(t,a, 'DisplayName','a');
plot(t,z, 'DisplayName','z');
yline(vd,'k-.','DisplayName','v_d');
axis([0 0.99*p*T 0 max(z)])
xlabel('t (s)')
title('Velocity, Accel, Separation')
legend()
hold off

subplot(2,2,2)
plot(t,z-1.8*v, 'DisplayName','safe heading');
hold on
plot(t, vd-v, 'DisplayName', 'Soft Constraint');
axis([0 0.99*p*T min(vd-v) max(z-1.8*v)])
title('Constraints')
xlabel('t (s)')
legend()

subplot(2,2,3)
plot(t, del, 'DisplayName','\delta')
title('Relaxation Param')
xlabel('t (s)')
legend()
axis([0 0.99*p*T min(del) max(del)])


%%% FUNCTIONS %%% TAKE CARE ON SIGNS
%--------------------------------------------------------------------------
function Fr = getFr(v, f0, f1, f2)
    Fr = f0 + f1 * v + f2 * v^2;
end

function Aclf = getAclf(m, v, vd)
    psi1 = 2/m * (v - vd);  
    
    Aclf = [ psi1, -1 ];
end

function bclf = getbclf(m, v, vd, eps, Fr)
    h = (v - vd); 
    psi0 = -2/m * h * Fr + eps * h^2;
    
    bclf = -psi0;
    
end

function Afcbf = getAfcbf() % will complete later
    Lg = 0;
    Afcbf = [ Lg, 0 ];
end

function bfcbf = getbfcbf() % will complete later
    Lf = 0;
    bfcbf = 0;
end

function bcbf = getbcbf(m, v, z, v0, Fr, gam)
    h = z - 1.8*v;
    
    B = -log( h / (1 + h) );
    
    Lf = -1/m*( 1.8*Fr + m*(v0 - v) )/( h * (1 + h) );
    
    bcbf = -Lf + gam/B;
end

function Acbf = getAcbf(m, v, z)
    h = z - 1.8*v;
    Lg = 1.8/m * 1/(h*(1 + h));
    
    Acbf = [ Lg, 0 ];
end