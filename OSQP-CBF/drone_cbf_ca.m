% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                      Ben Kaye, (c) 2020                             *
% *           Sequential QP with RCBF constraint for CA                 *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                       Using OSQP solver                             *
% *                 github.com/oxfordcontrol/osqp                       *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear

%%% PARAMETERS %%%


T_controller = 1e-3;
Ad = [1,0,0.000999995725421902,0,0,4.89519300035470e-06;0,1,0,0.000999995725421901,-4.89519300035470e-06,0;0,0,0.999987182106712,0,0,0.00978058201074243;0,0,0,0.999987182106712,-0.00978058201074243,0;0,0,0,0.00261084855543965,0.994003677759774,0;0,0,-0.00261084855543965,0,0,0.994003677759774];
% d_system.B([1:2, 4:5, 7:8], [4, 5])
Bd = [4.27457809838634e-09,0;0,4.27457809838634e-09;1.28178932883944e-05,0;0,1.28178932883944e-05;0,-0.00261084855543965;0.00261084855543965,0];


Au = Ad([3,4], :);
Bu = Bd([3,4], :);

p_o = [ 1.2; 1 ];

gamma = 1; % {}
g = 9.81; % {ms-2}
a_max = 1 * g; % {ms-2} Maximum decceleration allowed (design parameter)
delta = 0.3; % {m} Distance constraint
% T = 1e-3; % {s} Simulation step size
T = T_controller; % MUST BE 0.1s as dynamic system discretised to
sim_time = 10; % {s}


N = floor(sim_time/T);

p0 = [ 0; 0 ];
x_r = [ 2; 2; 0; 0; 0; 0 ];

Q = diag([ 1e10; 1e10; 0; 0; 0; 0 ]);
R = diag([ 1; 1 ]);

P = blkdiag(Q,R);

q = [ -2*Q'*x_r; 0; 0 ];

x = [ p0; 1e-5; 1e-5; 1e-5; 1e-5 ];

Aeq = [ eye(6), -Bd ];
leq = Ad*x;
ueq = leq;

[ Acbf, ucbf, ~ ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma);
lcbf = -inf;

Acc = [ zeros(2,6), eye(2) ];
lcc = -1 *  ones(2,1);
ucc = ones(2,1);


%%% SOLVER SETUP %%%

A = [ Aeq; Acbf; Acc ];
u = [ ueq; ucbf; ucc ];
l = [ leq; lcbf; lcc ];

solver = osqp;

solver.setup(P, q, A, l, u, 'warm_start',true,'verbose',false);

%%% SIMULATION %%%

errors = 0;
u_act = [ 0; 0 ];
p_t = zeros(2,N);
h_t = zeros(1,N);

for e = 1:N

    res = solver.solve();
    if res.info.status_val ~=1
        %%% INFEASIBILITY %%%
        errors = errors + 1;
        break
    else
        %%% SOLVED %%%
        x_star = res.x;
        u_act = x_star([7, 8]);
    end

    %%% SIMULATE TIME STEP %%%
    x = Ad*x + Bd * u_act;
    p_t(:,e) = x([1,2]);

    %%% UPDATE SOLVER %%%
    leq = Ad*x;
    ueq = leq;
    [ Acbf, ucbf, h ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma);
    h_t(e) = h;

    A = [ Aeq; Acbf; Acc ];
    u = [ ueq; ucbf; ucc ];
    l = [ leq; lcbf; lcc ];


    solver.update('Ax', nonzeros(A), 'u', u, 'l', l)
end

%%% PLOTTING %%%
t = 1:T:N;
t = t(1:e);
p_t = p_t(:,1:e);
h_t = h_t(1:e);

plot(p_t(1,:), p_t(2,:), '.');
viscircles(p_o', delta);



%%% FUNCTIONS %%%

function [ Acbf, ucbf, h ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma)
%Reciprocal Control Barrier Function constraint
    p = x([1, 2]);
    pdot = x([3, 4]);
    p_L = p_o - p;

    h = (p_L'*p_L) - ( delta + (p_L'*pdot)^2/2/a_max )^2;

    Bf = -log(h/(1+h));

    coeff = 2*(p_L'*pdot)/a_max*( delta + (p_L'*pdot)^2/2/a_max );
    LfBF = 1/h/(1+h)*( 2*(p_L'*pdot) + coeff*( p_L'*Au*x - (pdot'*pdot) ) );
    LgBF = 1/h/(1+h)*coeff*( p_L'*Bu );

    Acbf = [ zeros(1,6), LgBF ];
    ucbf = gamma/Bf - LfBF;        
end