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
T = 1e-1; % {s} sample period
sim_step = 1e-2; % {s} simulation step size
N = 100; % sim length
z0 = 10; % {m} initial distance
w0 = [ 0; 10;  z0 ]; % {m} {ms-1} 
w = w0; % state variable of car and lead
mu = 0;

%%% FORMULATE INITIAL OSQP MATRICES %%%
%--------------------------------------------------------------------------
Fr = getFr(w(2), f0, f1, f2);

P = 2 * diag( [ 1/m^2, psc ]);
q = -2 * [ Fr; 0 ]; % needs to be updated on each sim

x0 = [ 0; 0 ];

% leq = x0;
% ueq = x0;
% Aeq = eye(2);

lcc = ca * m * g;
ucc = cd * m * g;
Acc = [ 1 0 ];

lclf = -inf;
uclf = -getPsi0(m, w(2), vd, eps, Fr);
Aclf = [ getPsi1(m, w(2), vd), -1 ];

lcbf = -inf;
ucbf = getCBFfgam(m, w(2), w(3), v0, Fr, gam);
Acbf = [ getCBFg(m, w(2), w(3)), 0];

lfcbf = -inf;
ufcbf = getFCBFf();
Afcbf = [ getFCBFg() 0 ];

l = [ lclf; lcbf ];
u = [ uclf; ucbf ];
A = [ Aclf; Acbf ];

%%% using hard constraint or FCBF delete appropriately %%%
l = [ l; lcc ];
u = [ u; ucc ];
A = [ A; Acc ];

% l = [ l; lfcbf ];
% u = [ u; ufcbf ];
% A = [ A; Afcbf ];
%%% %%%

solver = osqp;
solver.setup(P, q, A, l, u, 'warm_start', true, 'verbose', false);

%%% SIMULATION %%%
%--------------------------------------------------------------------------

pastW = zeros(4, N*floor(T/sim_step));
p = 0;
errCount = 0;
u_in = 0;

for e = 1:N
    res = solver.solve();
    if ~strcmp(res.info.status, 'solved')
%         fprintf('ERROR\n')
        errCount = errCount + 1;
    else
        x = res.x;
        u_in = x(1); % control input relative
    end
    
    %%% forwards euler integration
    
    for k = 1:floor(T/sim_step)
        p = p + 1;
        pastW([1 2],p) = w(1:2);
        pastW(4,p) = w(3);
        %%% CHECK IF CAUSES WEIRDNESS
        Fr = getFr(w(2), f0, f1, f2);
%         u = Fr + m*mu;
%         u = mu;
        
        wdot =  [ w(2); -Fr; v0 - w(2) ] + ...
                [ 0; 1/m; 0 ] * u_in;
            
        w = w + wdot * sim_step;
        pastW(3, p) = wdot(2);
    end
    
%     % update inequalities %%% NOT OKAY
%     leq = x;
%     ueq = x;
    
    Fr = getFr(w(2), f0, f1, f2);

    lclf = -inf;
    uclf = -getPsi0(m, w(2), vd, eps, Fr);
    Aclf = [ getPsi1(m, w(2), vd), -1 ];

    lcbf = -inf;
    ucbf = getCBFfgam(m, w(2), w(3), v0, Fr, gam);
    Acbf = [ getCBFg(m, w(2), w(3)), 0];

    lfcbf = -inf;
    ufcbf = getFCBFf();
    Afcbf = [ getFCBFg() 0 ];
    
    l = [ lclf; lcbf ];
    u = [ uclf; ucbf ];
    A = [ Aclf; Acbf ];

    %%% using hard constraint or FCBF delete appropriately %%%
    l = [ l; lcc ];
    u = [ u; ucc ];
    A = [ A; Acc ];

    % l = [ l; lfcbf ];
    % u = [ u; ufcbf ];
    % A = [ A; Afcbf ];
    %%% %%%
    
    solver.update('Ax', nonzeros(A), 'l', l, 'u', u);    
end

%%% PLOT RESULTS %%%
%--------------------------------------------------------------------------
t = 0:sim_step:N*T;
x = [pastW(1,:), NaN];
v = [pastW(2, :), NaN];
a = [pastW(3, :), NaN];
z = [pastW(4, :), NaN];

figure(1)
hold on
plot(t,x, 'DisplayName','x');
plot(t,v, 'DisplayName','v');
plot(t,a, 'DisplayName','a');
plot(t,z, 'DisplayName','z');




%%% FUNCTIONS %%% TAKE CARE ON SIGNS
%--------------------------------------------------------------------------
function Fr = getFr(v, f0, f1, f2)
    Fr = f0 + f1 * v + f2 * v^2;
end

function psi1 = getPsi1(m, v, vd)
    psi1 = 2/m * (v - vd);    
end

function psi0 = getPsi0(m, v, vd, eps, Fr)
    h = (v - vd); 
    psi0 = -2/m * h * Fr + eps * h^2;
end

function fcbfG = getFCBFg()
    fcbfG = 0;
end

function fcbfF = getFCBFf()
    fcbfF = 0;
end

function cbfFgam = getCBFfgam(m, v, z, v0, Fr, gam)
    h = z - 1.8*v;
    
    B = - log( h / (1 + h) );
    
    cbfFgam =   + gam / B ...
                - 1/m*( 1.8*Fr + m*(v0 - v) )/( h * (1 - 1.8*v*z) );
end

function cbfG = getCBFg(m, v, z)
    h = z - 1.8*v;
    cbfG = 1.8/m / (1 + h) / h;
end