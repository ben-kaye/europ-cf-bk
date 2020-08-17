%$sign
%Ben Kaye, (c) 2020
%Sequential QP with RCBF constraint for CA
%$br
%Using OSQP solver
%github.com/oxfordcontrol/osqp
%$endsign

%%% PARAMETERS %%%

Ad = [      1,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
            0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
            0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
            0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
            0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
            0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   
];

% d_system.B([1:2, 4:5, 7:8], [4, 5])
Bd =    [   0.00370875152761323,    0;
            0,                      0.00370875152761323;
            0.105706677250003,      0;
            0,                      0.105706677250003;
            0,                      -0.193245482770890;
            0.193245482770890,      0                       
];


Au = Ad([1,2], :);
Bu = Bd([1,2], :);

p_o = [ 1; 1 ];

gamma = 1; % {}
g = 9.81; % {ms-2}
a_max = 0.2 * g; % {ms-2} Maximum decceleration allowed (design parameter)
delta = 0.3; % {m} Distance constraint
T = 1e-3; % {s} Simulation step size
sim_time = 10; % {s}


N = floor(sim_time/T);

p0 = [ 0; 0 ];
x_r = [ 2; 2; 0; 0; 0; 0 ];

Q = [ 1; 1; 0; 0; 0; 0 ];
R = [ 1; 1 ];

P = blkdiag(Q,R);

q = [ -2*Q'*x_r; 0; 0 ];

x = [ p0; 1e-5; 1e-5; 1e-5; 1e-5 ];

Aeq = [ eye(6), Bd ];
leq = Ad*x;
ueq = leq;

[ Acbf, ucbf ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma);
lcbf = -inf;

%%% SOLVER SETUP %%%

A = [ Aeq; Acbf ];
u = [ ueq; ucbf ];
l = [ leq; lcbf ];

solver = osqp;

solver.setup(P, q, A, l, u);

%%% SIMULATION %%%

errors = 0;
u_act = [ 0; 0 ];

for e = 1:N

    res = solver.solve();
    if res.info.status ~= -1
        %%% INFEASIBILITY %%%
        errors = errors + 1;
    else
        %%% SOLVED %%%
        x_star = res.x;
        u_act = x_star([7, 8]);
    end

    %%% SIMULATE TIME STEP %%%
    x = Ad*x + Bd * u_act;

    %%% UPDATE SOLVER %%%
    leq = Ad*x;
    [ Acbf, ucbf ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma);

    A = [ Aeq; Acbf ];
    u = [ ueq; ucbf ];
    l = [ leq; lcbf ];

    solver.update('Ax', nonzeros(A), 'u', u, 'l', l)
end

%%% FUNCTIONS %%%

function [ Acbf, ucbf ] = getCBF(Au, Bu, x, p_o, a_max, delta, gamma)
%Reciprocal Control Barrier Function constraint
    p = x([1, 2]);
    pdot = x([3, 4]);
    p_L = p_o - p;

    h = (p_L'*p_L) - ( delta + (p_L'*pdot)^2/2/a_max )^2;

    Bf = -log(h/(1+h));

    coeff = 2*(p_L'*pdot)/a_max*( delta + (p_L'*pdot)^2/2/a_max ) 
    LfBF = 1/h/(1+h)*( 2*(p_L'*pdot) + coeff*( p_L'*Au*x - (pdot'*pdot) ) );
    LgBF = 1/h/(1+h)*coeff*( p_L'*Bu );

    Acbf = [ zeros(1,6), LgBF ];
    ucbf = gamma/Bf - LfBF;        
end