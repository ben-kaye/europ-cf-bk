% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *                 Program by Ben Kaye (c) 2020                  *
% *                         EUROP Project                         *
% *                                                               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *        Using CrazyFlie model provided by Aren Karapet         *
% *                        and OSQP Solver                        *
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
ref = [ 0.6; 0.2; 0; 0; 0; 0 ];
x0 = [ -0.2; -0.5; 0; 0; 0; 0 ];
obj_xy = [ 0; 1 ];

% predictive horizon
N = 10;
static_states = 2;

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

% append dynamics of object = stationary
% Aeq = blkdiag( Aeq, eye(2) );

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
simtime = 10;
for i = 1 : simtime
    
    res = prob.solve();
    if ~strcmp(res.info.status, 'solved')
        error('OSQP could not solve')
    end
    
    xN = res.x;
    ctrls = xN( (N+1)*nx + 1 : end );
    ctrls = reshape(ctrls, [nu, N]);
    
    x0 = Ad * x0 + Bd * ctrls(:,1);
    hold on
    plot(x0(1),x0(2), 'r.', 'MarkerSize', 20)
%     axis([-1 1 -1 1])
    axis equal
    
    x1 = x0;
    for j = 2:4
        x1 = Ad * x1 + Bd * ctrls(:, j);
        hold on
        zscat = scatter(x1(1),x1(2), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
        zscat.MarkerFaceAlpha = .2;
        zscat.MarkerEdgeAlpha = .2;
        axis equal
    end
    
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l',l, 'u',u);
    
end






