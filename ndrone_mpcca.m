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
ref0 = [ 0.1; 1; 0; 0; 0; 0 ];
ref1 = [ -0.1; 0; 0; 0; 0; 0 ];
x0 = [ 0; 0; 0; 0; 0; 0 ];
x1 = [ 0; 1; 0; 0; 0; 0 ];
init_positions = {x0, x1};

x = [x0, x1];

min_dist = 5e-2; % Hard limit
maxmin_dist = 0.1;

% predictive horizon
M = 2;
N = 4;

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
P = blkdiag( kron(speye(N), Q ), QN, kron(speye(N), R));

% Augmented Linear Objective
q = [ repmat( -Q*ref0, N, 1 ); -QN*ref0; zeros(N*nu, 1)];
q = [ q; repmat( -Q*ref1, N, 1 ); -QN*ref1; zeros(N*nu, 1)];

% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

Aeq = [ Ax, Bu ];

% initialisation
leq = [];
for k = 1:M
    xt = init_positions{k}(:);
    leq = [ leq; -xt; zeros(N*nx, 1) ];
end

ueq = leq;
% Nx = leq;

% input and state constraints

% augmented state no.
Nstates = (N+1)*nx + N*nu;

Aineq = speye( M*Nstates );

lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1) ];

% Repeat matrices M times
Aeq = kron(speye(M), Aeq);
lineq = repmat(lineq, M, 1);
uineq = repmat(uineq, M, 1);
P = kron( speye(M), P) ;
% q = repmat(q, M, 1);


% Place holder, for CA constraints
Aplus = [ kron(speye(N+1), [ 1 1 0 0 0 0]), zeros(N+1, N*nu), -kron(speye(N+1), [ 1 1 0 0 0 0]), zeros(N+1, N*nu)];

% osqp constraints
A = [ Aeq; Aineq; Aplus];
l = [ leq; lineq; -inf*ones( N+1, 1) ];
u = [ ueq; uineq; +inf*ones( N+1, 1)];

% [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, min_dist, Nx, 1);




% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate
figure(1)
plot(x0(1), x0(2), 'bx', 'MarkerSize', 30)
hold on
plot(ref0(1), ref0(2), 'bx', 'MarkerSize', 30)
% plot(obj_xy(1), obj_xy(2), 'k+')
% viscircles(obj_xy', min_dist);
% viscircles(obj_xy', maxmin_dist);

xlabel('x (m)')
ylabel('y (m)')
title('Drone at time kT, T=0.1s')

error_count = 0;

xPast = x0;
h_bar = A(end-N:end, :);

simtime = 15;
nIt = 2;
init_its = 3;

% solve unconstrained prob
res = prob.solve();
xN = res.x;

% param = [ 0.1 0.5 1 ];
% vis_predicted_path
% for e = param
%     
%     relaxation_params = e * ones(N+1, 1);
%     
%     dist = min_dist * ( ones(N+1,1) - relaxation_params ) + maxmin_dist * relaxation_params;
%     
%     [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, dist, xN, 0);
%     
%     prob.update('l',l, 'Ax', nonzeros(A));
%     
%     res = prob.solve();
%     xN = res.x;
%     
% %     viscircles(obj_xy', dist(1));
% %     vis_predicted_path 
% end

for k = 1 : simtime
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
    
    
    
    ctrls0 = xN((N+1)*nx+1 : (N+1)*nx+N*nu );
    ctrls1 = xN(2*(N+1)*nx + N*nu + 1: 2*((N+1)*nx + N*nu));
    
    ctrls0 = reshape(ctrls0, [nu, N]);
    ctrls1 = reshape(ctrls1, [nu, N]);
    
    %in event of unsolvable path, carry on with previously calculated
    %optimal path
    ctrl0 = ctrls0(:, 1 + error_count);
    ctrl1 = ctrls1(:, 1 + error_count);
    
%     x0 = Ad * x0 + Bd * ctrl0;
%     x1 = Ad * x1 + Bd * ctrl1;

    ctrl_act = [ ctrl0, ctrl1];
    x = Ad * x + Bd * ctrl_act;
    
%     xPast = [xPast, x0];
    
%     x1 = x0;
    
    
%     h_bar = A(end-N:end, :);
    [ A, l, ~ ] = update_dist_constraint(A, l, M, N, nx, nu, maxmin_dist*ones((N+1),1), xN, 0);
%     h_bark = A(end-N:end, :);
    
    
    l(1:nx) = -x(:,1);
    u(1:nx) = -x(:,1);
    
    l((N+1)*nx + N*nu + 1:(N+1)*nx + N*nu + nx) = -x(:,2);
    u((N+1)*nx + N*nu + 1:(N+1)*nx + N*nu + nx) = -x(:,2);
    
    prob.update('l',l, 'u',u);
%     
    
    
    % if h_bar is unchanged, and eta and delta x are parallel for all k
    % rotate constrainst by 30 deg
%     if norm(full(h_bar - h_bark), 2) <= 1e-1
%         temp = reshape(xN(1:(N+1)*nx), nx, []);
%         x_ak = temp(1:2, :);
%         diff = [-eye(N), zeros(N,1)]+ [zeros(N,1), eye(N)];
%         xdif = x_ak * diff';
%         cp = norm(cross([xdif; zeros(1,N)], [etas(:,2:end); zeros(1,N)]),2);
%         if cp <= 2e-1
%            % rotate
%             [ A, l, etas ] = update_dist_constraint(A, l, N, nx, nu, o_states, dist, xN, 1);
%         end
%     end
    
    prob.update('Ax', nonzeros(A));
    
        
        
    hold on
    plot(x(1,1),x(2,1), '.', 'MarkerSize', 20, 'Color', [0 0 0.7])
    plot(x(1,2),x(2,2), '.', 'MarkerSize', 20, 'Color', [ 0 0.7 0.1])
%     for j = 2:3
%         x1 = Ad * x1 + Bd * ctrls(:, j);
%         zscat = scatter(x1(1),x1(2), 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%         zscat.MarkerFaceAlpha = .2;
%         zscat.MarkerEdgeAlpha = .2;
%     end
    axis equal
    
end

function [ A, l, etas ] = update_dist_constraint(A, l, M, N, nx, nu, min_dist, x, rotate)
   
%     min_dist = 0.3;
%     [ m, n ] = size(A);
    start_idx = 136;
    % WORK OUT FORMULA FOR
    
    
    n_states = (N+1)*nx + N*nu;
    
%     Aineq = blkdiag(speye( (N+1)*nx + N*nu ), zeros(2));
    % reset last N+1 constraints
%     A(end-N:end,:) = zeros(N+1, n);

    %NOT TRUE should be M Comb 2 
    
    
    l(end-(N+1)*nchoosek(M,2)+1:end) = min_dist;
    
    etas = zeros(2,N+1); 
    
    %%Figure out correct Size
    
    
    for idx_i = 1:M
        x_i = x( n_states*(idx_i-1)+1 : n_states*(idx_i-1) + (N+1)*nx);
        for idx_j = 1:M
            if( idx_j > idx_i)
                x_j = x( n_states*(idx_j-1)+1 : n_states*(idx_j-1) + (N+1)*nx);
                            
                for k = 1:N+1
                    x_ik = x_i( (k-1)*nx + 1 : k*nx );
                    x_jk = x_j( (k-1)*nx + 1 : k*nx );

                    x_ji = x_ik(1:2) - x_jk(1:2);
                    x_norm = norm(x_ji,2);
                    eta = x_ji / x_norm;
                    
                    % rotate 60 deg if commanded
                    if rotate ~= 0
                        R = [   cosd(60),   -sind(60);
                                sind(60),   cosd(60)    ];
                        eta = R * eta;
                    end
                    
                    % track etas for deadlock protection
%                     etas(:, k) = eta;

                    T = eta' * [ eye(2), -eye(2) ];
                    
                    A(start_idx+k+(idx_i+idx_j-3)*(N+1),...
                        [ ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 2,...
                        ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 2]) = T;
                    
                end
            end
        end
    end
end









