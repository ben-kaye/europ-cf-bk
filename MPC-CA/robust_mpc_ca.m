% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                           EUROP Project:                            *
% *                               MPC-CA                                *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *             CrazyFlie controller model by Aren Karapet              *
% *                          Using OSQP Solver                          *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear,clc

% PARAMETERS
VIS_ON = 1;
N = 15;
N_OBJ_STATES = 2;
MIN_DIST = 0.05;
IDEAL_DIST = 0.7;
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 20;

% initial conditions
ref = [ 2; 0.5; 0; 0; 0; 0 ];
x0 = [ 0.1; 0.1; 0; 0; 0; 0 ];
obj_xy = [ 1; 0.5 ];

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
% augmented state no.
Nstates = (N+1)*nx + N*nu + N_OBJ_STATES;    
    
   
    


% predictive horizon


% constraints
umin = [ -1; -1 ];
umax = [ 1; 1 ];
xmin = [ -inf; -inf; -2; -2; -inf; -inf ];
xmax = [ inf; inf; +2; +2; inf; inf ];


% LQ objective
Q = diag([1 1 0.1 0.1 0 0]);
QN = Q;
R = 0.01 * eye(nu);

% Augmented Quadratic Objective
P = blkdiag( kron(speye(N), Q ), QN, zeros(N_OBJ_STATES), kron(speye(N), R));
% Augmented Linear Objective
q = [ repmat( -Q*ref, N, 1 ); -QN*ref; zeros(N_OBJ_STATES,1); zeros(N*nu, 1)];
% q = zeros( nx*(N+1) + N*nu, 1 );


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

% append dynamics of object = stationary
Ax = blkdiag(Ax, -speye(N_OBJ_STATES));
Bu = [Bu; zeros(N_OBJ_STATES, size(Bu, 2))];

Aeq = [ Ax, Bu ];

leq = [ -x0; zeros(N*nx, 1); -obj_xy ];
ueq = leq;
xN = leq;


% input and state constraints

Aineq = speye( Nstates );

lineq = [ repmat( xmin, N+1, 1); -inf*ones(N_OBJ_STATES, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); +inf*ones(N_OBJ_STATES, 1); repmat( umax, N, 1) ];

Aconstr = [ kron(speye(N+1), [ 1 1 0 0 0 0]), -ones(N+1, N_OBJ_STATES), zeros(N+1, N*nu) ];

% osqp constraints
A = [ Aeq; Aineq; Aconstr ];
l = [ leq; lineq; -inf*ones((N+1), 1) ];
u = [ ueq; uineq; inf*ones((N+1), 1)];


% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate

if (VIS_ON)
    figure(1)
    plot(x0(1), x0(2), 'bx', 'MarkerSize', 30, 'DisplayName', 'Start')
    hold on
    plot(ref(1), ref(2), 'rx', 'MarkerSize', 30, 'DisplayName', 'End')
    plot(obj_xy(1), obj_xy(2), 'r+', 'DisplayName', 'Radius of Avoidance')
    viscircles(obj_xy', MIN_DIST)

    xlabel('x (m)')
    ylabel('y (m)')
    title('Drone at time kT, T=0.1s')
end

error_count = 0;

x_hist = zeros((N+1)*nx, SIM_TIME);

% solve unconstrained problem
res = prob.solve();
xN = res.x;

vis_predicted_path
for e = EPSILON
    
    relaxation_params = e * ones(N+1, 1);
    
    dists = MIN_DIST * ( ones(N+1,1) - relaxation_params ) + IDEAL_DIST * relaxation_params;
    
    [ A, l ] = compute_linear_constrs(A, l, N, nx, nu, N_OBJ_STATES, dists, xN);
    
    prob.update('l',l, 'Ax', nonzeros(A));
    
    res = prob.solve();
    xN = res.x;
    
%     viscircles(obj_xy', dists(1));
    vis_predicted_path    
end


for i = 1 : SIM_TIME
    
    res = prob.solve();
    % record result if solved
    if ~strcmp(res.info.status, 'solved')
            error_count = error_count + 1;                
    else 
        error_count = 0;
        xN = res.x;
    end 
   
    % get control actions
    ctrls = xN( (N+1)*nx + N_OBJ_STATES + 1 : end );
    ctrls = reshape(ctrls, [nu, N]);
    
    % in event of unsolvable path, carry on with previously calculated
    % optimal path
    ctrl = ctrls(:, 1 + error_count);
    
    % activate control
    x0 = Ad * x0 + Bd * ctrl;
    
    % store state
    x_hist(1:nx,i) = x0;
    
    % compute linear constraints
    [ A, l ] = compute_linear_constrs(A, l, N, nx, nu, N_OBJ_STATES, dists, xN);
    
    % update initial conds and constraints
    l(1:nx) = -x0;
    u(1:nx) = -x0;
    
    prob.update('l',l, 'u',u);
    prob.update('Ax', nonzeros(A));
  
  
end

% plot results
plot(x_hist(1,:), x_hist(2,:), '.', 'MarkerSize', 20, 'Color', 1/255*[73, 146, 214], 'DisplayName', 'Path Taken')
axis equal

function [ A, l ] = compute_linear_constrs(A, l, N, nx, nu, no, min_dists, x)
    [ m, n ] = size(A);
    start_idx = m - N - 1; % i = start_idx + k

    % get obj coords
    obj_xy = x((N+1)*nx + 1:(N+1)*nx + no);

    % assign dists to <= cstr
    l(end-N:end) = min_dists;

    % for each time step
    for k = 1:N+1
        % get state at k
        state_k = x( (k-1)*nx + 1 : k*nx );

        % get vec to obj
        x_ij = state_k(1:2) - obj_xy;
        x_norm = norm(x_ij,2);

        eta = x_ij / x_norm;    

        T = eta' * [ eye(no), -eye(no) ];

        % assign to A where Ax>=l
        A( start_idx + k,...
            [ (k-1)*nx + 1 : (k-1)*nx + 2, (N+1)*nx + 1 : (N+1)*nx + no ] )...
            = T;
    end
end









