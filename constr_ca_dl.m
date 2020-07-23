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
Ad = [  0,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
        0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
        0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
        0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
        0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
        0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   ];

% d_system.B([1:2, 4:5, 7:8], [4, 5])
Bd = [  0.00370875152761323,    0;
        0,                      0.00370875152761323;
        0.105706677250003,      0;
        0,                      0.105706677250003;
        0,                      -0.193245482770890;
        0.193245482770890,      0                       ];

[nx, nu] = size(Bd);    
    
% initial conditions
ref = [ 1.9; 0.4; 0; 0; 0; 0 ];
x0 = [ 0; 0; 0; 0; 0; 0 ];

obj_xy = [0.9; 0.4];
min_dist = 5e-2; % Hard limit
maxmin_dist = 0.5;

% predictive horizon
N = 200;
o_states = 2;

%lin interp between min and ideal dist
relaxation_params = zeros(N+1, 1);

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
P = blkdiag( kron(speye(N), Q ), QN, zeros(o_states), kron(speye(N), R));
% Augmented Linear Objective
q = [ repmat( -Q*ref, N, 1 ); -QN*ref; zeros(o_states,1); zeros(N*nu, 1)];
% q = zeros( nx*(N+1) + N*nu, 1 );


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

% append dynamics of object = stationary
Ax = blkdiag(Ax, -speye(o_states));
Bu = [Bu; zeros(o_states,size(Bu, 2))];

Aeq = [ Ax, Bu ];

% should it be -ve or +ve
leq = [ -x0; zeros(N*nx, 1); -obj_xy ];
ueq = leq;
Nx = leq;

% input and state constraints

% augmented state no.
Nstates = (N+1)*nx + N*nu + o_states;

Aineq = speye( Nstates );

lineq = [ repmat( xmin, N+1, 1); -inf*ones(o_states, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); +inf*ones(o_states, 1); repmat( umax, N, 1) ];

% [ Aineq, lineq ] = update_dist_constraint(Aineq, lineq, N, nx, nu, min_dist, [x0; zeros(Nstates - nx -nu -2, 1)]);


Aplus = [ kron(speye(N+1), [ 1 1 0 0 0 0]), -ones(N+1, o_states), zeros(N+1, N*nu) ];

% osqp constraints
A = [ Aeq; Aineq; Aplus ];
l = [ leq; lineq; -inf*ones((N+1), 1) ];
u = [ ueq; uineq; inf*ones((N+1), 1)];

% [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, min_dist, Nx, 1);

% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate
figure(1)
plot(x0(1), x0(2), 'bx', 'MarkerSize', 30)
hold on
plot(ref(1), ref(2), 'bx', 'MarkerSize', 30)
plot(obj_xy(1), obj_xy(2), 'k+')
% viscircles(obj_xy', min_dist);
viscircles(obj_xy', maxmin_dist);

xlabel('x (m)')
ylabel('y (m)')
title('Drone at time kT, T=0.1s')

error_count = 0;

xPast = x0;
h_bar = A(end-N:end, :);

simtime = 20;
nIt = 2;
init_its = 3;

% solve unconstrained prob
res = prob.solve();
xN = res.x;

param = [ 0.1 0.5 1 ];
vis_predicted_path
for e = param
    
    relaxation_params = e * ones(N+1, 1);
    
    dist = min_dist * ( ones(N+1,1) - relaxation_params ) + maxmin_dist * relaxation_params;
    
    [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, dist, xN, 0);
    
    prob.update('l',l, 'Ax', nonzeros(A));
    
    res = prob.solve();
    xN = res.x;
    
    viscircles(obj_xy', dist(1));
    vis_predicted_path
   
    
    
    
end

for i = 1 : simtime
%     for j = 1:nIt
        res = prob.solve();
        if ~strcmp(res.info.status, 'solved')
%             error('OSQP could not solve')
            % Try rotating linear distance constraint, failing that, make
            % first move
%             if error_count == 0                
%                 [ A, l ] = update_dist_constraint(A, l, N, nx, nu, o_states, min_dist, xN, 1);
%                 prob.update('l',l, 'u',u);
%                 prob.update('Ax', nonzeros(A));
%                 res = prob.solve();
%             end
            
            if ~strcmp(res.info.status, 'solved')
                error_count = error_count + 1;
            else
                error_count = 0;
                xN = res.x;
            end
            
            
            
        else 
            error_count = 0;
            xN = res.x;
        end 
        
        
        
%         [ A, l ] = update_dist_constraint(A, l, N, nx, nu, o_states, min_dist, xN, 0);
    
%         l(1:nx) = -x0;
%         u(1:nx) = -x0;
    %     prob.update('l',l, 'u',u, 'Ax',A); updating A throws an error
%         prob.update('l',l, 'u',u);
%         prob.update('Ax', nonzeros(A));
        
        
        
%     end
    
    
    
    ctrls = xN( (N+1)*nx + o_states + 1 : end );
    ctrls = reshape(ctrls, [nu, N]);
    
    %in event of unsolvable path, carry on with previously calculated
    %optimal path
    ctrl = ctrls(:, 1 + error_count);
    
    x0 = Ad * x0 + Bd * ctrl;
    
    
    xPast = [xPast, x0];
    
    x1 = x0;
    
    
    h_bar = A(end-N:end, :);
    [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, dist, xN, 0);
    h_bark = A(end-N:end, :);
    
    
    l(1:nx) = -x0;
    u(1:nx) = -x0;
%     prob.update('l',l, 'u',u, 'Ax',A); updating A throws an error
    prob.update('l',l, 'u',u);
%     
    
    
    %if h_bar is unchanged, and eta and delta x are parallel for all k
%     rotate constrainst by 30 deg
    if norm(full(h_bar - h_bark), 2) <= 1e-1
        temp = reshape(xN(1:(N+1)*nx), nx, []);
        x_ak = temp(1:2, :);
        diff = [-eye(N), zeros(N,1)]+ [zeros(N,1), eye(N)];
        xdif = x_ak * diff';
        cp = norm(cross([xdif; zeros(1,N)], [etas(:,2:end); zeros(1,N)]),2);
        if cp <= 2e-1
           % rotate
            [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, dist, xN, 1);
        end
    end
    
    prob.update('Ax', nonzeros(A));
    
        
        
    hold on
    plot(x0(1),x0(2), 'r.', 'MarkerSize', 20)
    
%     for j = 2:3
%         x1 = Ad * x1 + Bd * ctrls(:, j);
%         zscat = scatter(x1(1),x1(2), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%         zscat.MarkerFaceAlpha = .2;
%         zscat.MarkerEdgeAlpha = .2;
%     end
    axis equal
    
end

function [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, no, min_dist, x, rotate)
   
%     min_dist = 0.3;
    [ m, n ] = size(A);
    start_idx = m - N - 1; % i = start_idx + k
    
    obj_xy = x((N+1)*nx + 1:(N+1)*nx + no);
    
%     Aineq = blkdiag(speye( (N+1)*nx + N*nu ), zeros(2));
    % reset last N+1 constraints
    A(end-N:end,:) = zeros(N+1, n);
    
%     lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1); zeros(2,1) ];
    l(end-N:end) = min_dist;
   
%     uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1); zeros(2,1) ];
%     uineq = [ unieq; inf * ones(N+1, 1) ];
    
    etas = zeros(2,N+1);
    for k = 1:N+1
        
        state_k = x( (k-1)*nx + 1 : k*nx );
        x_ij = state_k(1:2) - obj_xy;
        x_norm = norm(x_ij,2);
        
        eta = x_ij / x_norm;
        if rotate ~= 0
            R = [   cosd(60),   -sind(60);
                    sind(60),   cosd(60)    ];
            eta = R * eta;
        end
        etas(:, k) = eta;
        
        
        
        
        T = eta' * [ eye(no), -eye(no) ];
        A( start_idx + k, [ (k-1)*nx + 1 : (k-1)*nx + 2, (N+1)*nx + 1 : (N+1)*nx + no ] ) = T;
    end
end









