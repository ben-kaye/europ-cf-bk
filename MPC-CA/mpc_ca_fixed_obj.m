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

opts =  optimset('Display','off');

% PARAMETERS
VIS_ON = 1;
N = 5;
% N_OBJ_STATES = 2;
MIN_DIST = 0.05;
IDEAL_DIST = 0.9;
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 30;

% initial conditions
ref = [ 2; 0.5; 0; 0; 0; 0 ];
x0 = [ 0.1; 0.1; 0; 0; 0; 0 ];
p_o = [ 1; 0.5 ];

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
% Nstates = (N+1)*nx + N*nu + N_OBJ_STATES;    
    
   
    


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
P = blkdiag( kron(speye(N), Q ), QN, kron(speye(N), R));
% Augmented Linear Objective
q = [ repmat( -Q*ref, N, 1 ); -QN*ref; zeros(N*nu, 1)];
% q = zeros( nx*(N+1) + N*nu, 1 );


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

Aeq = [ Ax, Bu ];

leq = [ -x0; zeros(N*nx, 1) ];
ueq = leq;
xN = leq;


% input and state constraints

Aineq = speye((N+1)*nx + N*nu);

lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1) ];

[ A_dist, ~ ] = dist_constraints(N, nx, nu, p_o, zeros(N, 1), -leq);
l_dist = -inf*ones(N, 1);
u_dist = inf*ones(N,1);

% osqp constraints
% A = [ Aeq; Aineq; A_dist ];
% l = [ leq; lineq; l_dist ];
% u = [ ueq; uineq; u_dist ];


% osqp init and setup
% prob = osqp;
% prob.setup( P, q, A, l, u, 'warm_start', false, 'verbose', false );

% simulate

if (VIS_ON)
    figure(1)
    plot(x0(1), x0(2), 'bx', 'MarkerSize', 30, 'DisplayName', 'Start')
    hold on
    plot(ref(1), ref(2), 'rx', 'MarkerSize', 30, 'DisplayName', 'End')
    plot(p_o(1), p_o(2), 'r+', 'DisplayName', 'Radius of Avoidance')
    viscircles(p_o', IDEAL_DIST, 'LineStyle', '-.')

    xlabel('x (m)')
    ylabel('y (m)')
    title('Drone at time kT, T=0.1s')
end

error_count = 0;

x_hist = zeros((N+1)*nx, SIM_TIME);

% solve unconstrained problem
% res = prob.solve();
% xN = res.x;

xN = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);

% vis_predicted_path
for e = EPSILON
    
    relaxation_params = e * ones(N, 1);
    
    dists = MIN_DIST * ( ones(N,1) - relaxation_params ) + IDEAL_DIST * relaxation_params;
    
    [ A_dist, u_dist ] = dist_constraints(N, nx, nu, p_o, dists, xN);
    
%     A(end-N+1:end) = A_dist;
%     u(end-N+1:end) = u_dist;
    
%     prob.update('u', u, 'Ax', nonzeros(A));
    
%     res = prob.solve();
%     xN = res.x;

    xN = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);
    
%     viscircles(obj_xy', dists(1));
%     vis_predicted_path    
end


for i = 1 : SIM_TIME
    
%     res = prob.solve();
    xN = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);
    % record result if solved
%     if ~strcmp(res.info.status, 'solved')
%             error_count = error_count + 1;                
%     else 
%         error_count = 0;
%         xN = res.x;
%     end 
   
    % get control actions
    ctrls = xN( (N+1)*nx + 1 : end );
    ctrls = reshape(ctrls, nu, []);
    
    % in event of unsolvable path, carry on with previously calculated
    % optimal path
    ctrl = ctrls(:, 1 + error_count);
    
    % activate control
    x0 = Ad * x0 + Bd * ctrl;
    
    % store state
    x_hist(1:nx,i) = x0;
    
    % compute linear constraints
    [ A_dist, u_dist ] = dist_constraints(N, nx, nu, p_o, dists, xN);
    
    % update initial conds and constraints
%     l(1:nx) = -x0;
%     u(1:nx) = -x0;
    
    leq(1:nx) = -x0;
    ueq(1:nx) = -x0;

%     A(end-N+1:end) = A_dist;
%     u(end-N+1:end) = u_dist;
%     
%     prob.update('l',l, 'u',u);
%     prob.update('Ax', nonzeros(A));
  
  
end

% plot results
if (VIS_ON)
    plot(x_hist(1,:), x_hist(2,:), '-^', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', 1/255*[73, 146, 214], 'DisplayName', 'Path Taken')
    axis equal
end

function [ A_dist, u_dist ] = dist_constraints(Np, nx, nu, p_o, min_dists, x_mpc)
    x_mpc = x_mpc(1:(Np+1)*nx);
    x_mpc = reshape(x_mpc, nx, []);

    % past pos at k
    p_k = x_mpc([1,2], :);
    
    p_ox = p_k(:, 2:end) - p_o;
    eta = p_ox./vecnorm(p_ox);
    
    u_dist = -min_dists - eta'*p_o;
    
    A_dist = [ zeros(Np, nx), kron(speye(Np), [ 1, 1, zeros(1, nx - 2)]), zeros(Np, nu*Np) ];
    A_dist(A_dist ~= 0) = -eta;    
end









