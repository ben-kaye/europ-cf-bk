% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *                 Program by Ben Kaye (c) 2020                  *
% *            Using model provided by Aren Karapetyan            *
% *                         EUROP Project                         *
% *                                                               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *                    Implements OSQP solver                     *
% *                                                               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear,clc

% System State-Space

% d_system.A([1:2, 4:5, 7:8], [1:2, 4:5, 7:8])
Ad = [   1,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
        0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
        0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
        0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
        0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
        0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   ];

% d_system.B([1:2, 4:5, 7:8], [4, 5])
Bd = [   0.00370875152761323,    0;
        0,                      0.00370875152761323;
        0.105706677250003,      0;
        0,                      0.105706677250003;
        0,                      -0.193245482770890;
        0.193245482770890,      0                       ];

[nx, nu] = size(Bd);    
    
% initial conditions
ref = [ 3; 0; 0; 0; 0; 0 ];
x0 = [ -4; 3; 0; 0; 0; 0 ];
obj_xy = [ 0; 1 ];

% predictive horizon
N = 10;

% constraints
umin = [ -1; -1 ];
umax = [ 1; 1 ];
xmin = [ -inf; -inf; -2; -2; -inf; -inf ];
xmax = [ inf; inf; +2; +2; inf; inf ];


% LQ objective
Q = diag([1 1 0.1 0.1 0 0]);
QN = Q;
R = 0.01 * eye(2);

% Augmented Quadratic Objective
P = blkdiag( kron(speye(N), Q ), QN, kron(speye(N), R) );
% Augmented Linear Objective
q = [ repmat( -Q*ref, N, 1 ); -QN*ref; zeros(N*nu, 1) ];
% q = zeros( nx*(N+1) + N*nu, 1 );


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );
Aeq = [ Ax, Bu ];
leq = [ -x0; zeros(N*nx, 1) ];
ueq = leq;

% input and state constraints
Aineq = speye( (N+1)*nx + N*nu );
lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1) ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1) ];

% osqp constraints
A = [ Aeq; Aineq ];
l = [ leq; lineq ];
u = [ ueq; uineq ];

% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate
figure(1)
plot(x0(1), x0(2), 'b+')
hold on
plot(ref(1), ref(2), 'b+')
simtime = 100;
for i = 1 : simtime
    
    res = prob.solve();
    if ~strcmp(res.info.status, 'solved')
        error('OSQP could not solve')
    end
    
    ctrl = res.x( (N+1)*nx+1 : (N+1)*nx + nu );
    x0 = Ad * x0 + Bd * ctrl;
    
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l',l, 'u',u);
    
    hold on
    plot(x0(1),x0(2), 'ro')
    
end










