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
ref = [ 0.6; 0.2; 0; 0; 0; 0 ];
x0 = [ -0.2; -0.5; 0; 0; 0; 0 ];
obj_xy = [ -0.9; -0.9 ];
min_dist = 0.1;

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
P = blkdiag( kron(speye(N), Q ), QN, kron(speye(N), R), zeros(2));
% Augmented Linear Objective
q = [ repmat( -Q*ref, N, 1 ); -QN*ref; zeros(N*nu, 1); zeros(2,1) ];
% q = zeros( nx*(N+1) + N*nu, 1 );


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);

Bu = kron( [ sparse(1, N); speye(N) ], Bd );
Aeq = [ Ax, Bu ];

% append dynamics of object = stationary
Aeq = blkdiag( Aeq, speye(2) );

leq = [ -x0; zeros(N*nx, 1); -obj_xy ];
ueq = leq;

% input and state constraints

% augmented state no.
Nstates = (N+1)*nx + N*nu + 2;

Aineq = speye( (N+1)*nx + N*nu + 2 );
Aineq = [ Aineq; zeros( N+1, Nstates )];

lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1); -inf*ones((N+1)+2, 1) ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1); inf*ones((N+1)+2, 1) ];

% [ Aineq, lineq ] = update_dist_constraint(Aineq, lineq, N, nx, nu, min_dist, [x0; zeros(Nstates - nx -nu -2, 1)]);

% osqp constraints
A = [ Aeq; Aineq ];
l = [ leq; lineq ];
u = [ ueq; uineq ];

% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate
figure(1)
plot(x0(1), x0(2), 'bx', 'MarkerSize', 30)
hold on
plot(ref(1), ref(2), 'bx', 'MarkerSize', 30)
plot(obj_xy(1), obj_xy(2), 'k+')
viscircles(obj_xy', min_dist)
simtime = 100;
for i = 1 : simtime
    
    res = prob.solve();
    if ~strcmp(res.info.status, 'solved')
        error('OSQP could not solve')
    end
    
    xN = res.x;
    ctrl = xN( (N+1)*nx+1 : (N+1)*nx + nu );
    x0 = Ad * x0 + Bd * ctrl;
    
    
%     [ Aineq, lineq ] = update_dist_constraint(Aineq, lineq, N, nx, nu, min_dist, xN);
%     A = [ Aeq; Aineq];
%     l = [ leq; lineq];
    
    
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    prob.update('l',l, 'u',u, 'Ax',A);
    
    hold on
    plot(x0(1),x0(2), 'r.', 'MarkerSize', 20)
%     axis([-1 1 -1 1])
    axis equal
    
end

function [ Aineq, lineq ] = update_dist_constraint(Aineq, lineq, N, nx, nu, min_dist, x)
   
%     min_dist = 0.3;
    
    start_idx = (N+1)*nx + N*nu + 2; % i = start_idx + k
    obj_xy = x(end-1:end);
    
%     Aineq = blkdiag(speye( (N+1)*nx + N*nu ), zeros(2));
    % add in N+1 object dist constraints
    Aineq(start_idx+1:end,:) = zeros(N+1, start_idx);
    
%     lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1); zeros(2,1) ];
    lineq(start_idx+1:end) = min_dist * ones( N+1, 1 );
   
%     uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1); zeros(2,1) ];
%     uineq = [ unieq; inf * ones(N+1, 1) ];
    
    
    for k = 1:N+1
        
        state_k = x( (k-1)*nx + 1 : k*nx );
        x_ij = state_k(1:2) - obj_xy;
        x_norm = norm(x_ij,2);
        eta = x_ij / x_norm;
        
        
        T = eta' * [ eye(2), -eye(2) ];
        Aineq( start_idx + k, [ (k-1)*nx + 1 : (k-1)*nx + 2, end-1 : end ] ) = T;
%         lineq( start_idx + k ) = lineq( start_idx + k ) - x_norm; % + 0
%         lineq( start_idx + k ) = lineq( start_idx + k );
    end
end










