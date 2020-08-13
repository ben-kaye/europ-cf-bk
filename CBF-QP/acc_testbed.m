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
% psc = exp(-5); % relaxation weight

psc = 1; % TO BE OVERRIDDEN
psc_no_acc_lim = 5e1; % big number needed to keep delta down
psc_acc_lim = 200;

T = 1e-3; % {s} sample period
sim_step = 1e-3; % {s} simulation step size
sim_time = 30; % {s}
N = sim_time/T; % sim length
z0 = 100; % {m} initial distance
xstart = 900; % {m} initial x-coord a bit unused)
vstart = 9; % {ms-1} initial speed
w0 = [ xstart; vstart;  z0 ]; % {m} {ms-1} {m}
w = w0; % state :: w(1) = x; w(2) = v; w(3) = z (dist between lead and car)
mu = 0;

FCBF_ON = 1; %{1:ON 0:OFF} Whether to include hard constraint on accel, and appropriate CBF

if FCBF_ON
    psc = psc_acc_lim;
else 
    psc = psc_no_acc_lim;
end

%%% FORMULATE INITIAL OSQP MATRICES %%%
%--------------------------------------------------------------------------
Fr = getFr(w(2), f0, f1, f2);

P = 2 * diag( [ 1/m^2, psc ]);
q = -2/m * [ Fr; 0 ]; % needs to be updated on each sim

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
ufcbf = getbfcbf(m,w(2),w(3),v0,cd,g,Fr);
Afcbf = getAfcbf(m,w(2),w(3),v0,cd,g);

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
        pastW(6, p) = res.info.obj_val;
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
    ufcbf = getbfcbf(m,w(2),w(3),v0,cd,g,Fr);
    Afcbf = getAfcbf(m,w(2),w(3),v0,cd,g);

    
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
x = [pastW(1,:), NaN];
v = [pastW(2, :), NaN];
a = [pastW(3, :), NaN];
z = [pastW(4, :), NaN];
del = [pastW(5,:), NaN];
obj = [pastW(6,:), NaN];
figure(1)

subplot(2,2,1)
% plot(t,x, 'DisplayName','x');

plot(t,v, 'DisplayName','v', 'LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
hold on
plot(t,z, 'DisplayName','z', 'LineWidth',1.5,'Color', 1/255*[252, 169, 0]);
yline(vd,'k-.','DisplayName','v_d');
yline(v0, 'k-.','DisplayName','v_0');
axis([0, 0.99*p*T, min(v)-1, max(z)+1])
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
axis([0, 0.99*p*T, min(vd-v)-1, max(z-1.8*v)+1])
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
axis([0, 0.99*p*T, min(del)-1, max(del)+1])
sgtitle(sprintf('Relaxation weight %f',psc));

subplot(2,2,4)
% plot(t, obj, 'DisplayName', 'J');
% title('Objective function')

plot(t,a, 'DisplayName','a','LineWidth',1.5,'Color', 1/255*[0, 157, 252]);
yline(-cd*g, 'DisplayName', 'a_{min}')
yline(ca*g, 'DisplayName', 'a_{max}')
title('Acceleration')
xlabel('t (s)')
ylabel('a (ms^{-2})')
legend()


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

function Afcbf = getAfcbf(m,v, z, v0, cd, g) % will complete later
    hF = -1.8*v - 1/2*( v0 - v )^2 /cd/g + z;
    
    Lg = -1/m/hF^2 * (  (v0 - v)/cd/g - 1.8 );
    
    Afcbf = [ Lg, 0 ];
end

function bfcbf = getbfcbf(m, v, z, v0, cd, g, Fr) % will complete later
    hF = -1.8*v - 1/2*( v0 - v )^2 /cd/g + z;
    Lf = 1/hF^2 * ( Fr/m * ( ( v0 - v )/cd/g - 1.8 ) - (v0 -v) );
    bfcbf = hF - Lf;
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