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

% PARAMETERS
VIS_ON = 1; % display graphs
N = 15; % Horizon
M = 2; % Number of drones
MIN_DIST = 0.05; % Absolute min
IDEAL_DIST = 0.3; % Ideal range
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 20; % Number of iterations in simulation (linear sim)

% initial conditions
ref = [ 1,  1;
        -1, 1;
        0,  0;
        0,  0;
        0,  0;
        0,  0 ];
x0 = [  -1, -1;
        1,  -1;
        0,  0;
        0,  0;
        0,  0;
        0,  0   ];

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
 Nstates = (N+1)*nx + N*nu; 

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
q = [   repmat( -Q*ref(:,1), N, 1 ); -QN*ref(:,1); zeros(N*nu, 1); ...
        repmat( -Q*ref(:,2), N, 1 ); -QN*ref(:,2); zeros(N*nu, 1)   ];


% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

Aeq = [ Ax, Bu ];

leq = [ -x0(:,1); zeros(N*nx, 1); -x0(:,2); zeros(N*nx,1) ];
ueq = leq;
xN = leq;


lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1) ];


Aeq = kron(speye(M), Aeq);
lineq = repmat(lineq, M, 1);
uineq = repmat(uineq, M, 1);
P = kron( speye(M), P) ;


% input and state constraints

Aineq = speye( M * Nstates );

% Aconstr = [ kron(speye(N+1), [ 1 1 0 0 0 0]), -ones(N+1, N_OBJ_STATES), zeros(N+1, N*nu) ];

% osqp constraints
A = [ Aeq; Aineq ];
l = [ leq; lineq ];
u = [ ueq; uineq ];


% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate

if (VIS_ON)
    figure(1)
    plot(x0(1,1), x0(2,1), '+', 'MarkerSize', 30, 'DisplayName', 'A1 Start', 'Color', 1/255*[73, 146, 214])
    hold on
    plot(x0(1,2), x0(2,2), '+', 'MarkerSize', 30, 'DisplayName', 'A2 Start', 'Color', 1/255*[214, 88, 149])
    plot(ref(1,1), ref(2,1), 'x', 'MarkerSize', 30, 'DisplayName', 'A1 End', 'Color', 1/255*[73, 146, 214])
    plot(ref(1,2), ref(2,2), 'x', 'MarkerSize', 30, 'DisplayName', 'A2 End', 'Color', 1/255*[214, 88, 149])


    xlabel('x (m)')
    ylabel('y (m)')
    title('Drone at time kT, T=0.1s')
end

error_count = 0;

x_hist = zeros((N+1)*nx, SIM_TIME);

% solve unconstrained problem
res = prob.solve();
xN = res.x;

% for e = EPSILON
%     
%     relaxation_params = e * ones(N+1, 1);
%     
%     dists = MIN_DIST * ( ones(N+1,1) - relaxation_params ) + IDEAL_DIST * relaxation_params;
%     
%     [ A, l ] = compute_linear_constrs(A, l, N, nx, nu, N_OBJ_STATES, dists, xN);
%     
%     prob.update('l',l, 'Ax', nonzeros(A));
%     
%     res = prob.solve();
%     xN = res.x;
%     
% %     viscircles(obj_xy', dists(1));
% %     vis_predicted_path    
% end


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
    ctrls0 = xN((N+1)*nx+1 : (N+1)*nx+N*nu );
    ctrls1 = xN(2*(N+1)*nx + N*nu + 1: 2*((N+1)*nx + N*nu));
    
    ctrls0 = reshape(ctrls0, [nu, N]);
    ctrls1 = reshape(ctrls1, [nu, N]);
    
    %in event of unsolvable path, carry on with previously calculated
    %optimal path
    ctrl0 = ctrls0(:, 1 + error_count);
    ctrl1 = ctrls1(:, 1 + error_count);
    
    ctrl_act = [ ctrl0, ctrl1 ];
    
    % activate control
    x0 = Ad * x0 + Bd * ctrl_act;
    
    % store state
    x_hist(1:2*nx,i) = reshape(x0, 2*nx, []);
    
    % compute linear constraints
%     [ A, l ] = compute_linear_constrs(A, l, N, nx, nu, N_OBJ_STATES, dists, xN);
    
    % update initial conds and constraints
    
    xtemp = reshape(x0, 2*nx, []);
    l([1:nx, Nstates+1:Nstates+nx]) = -xtemp;
    u([1:nx, Nstates+1:Nstates+nx]) = -xtemp;
    
    prob.update('l',l, 'u',u);
%     prob.update('Ax', nonzeros(A));
%   
  
end

% plot results
plot(x_hist(1,:), x_hist(2,:), '.', 'MarkerSize', 20, 'Color', 1/255*[73, 146, 214], 'DisplayName', 'A1 Path Taken')
plot(x_hist(7,:), x_hist(8,:), '.', 'MarkerSize', 20, 'Color', 1/255*[214, 88, 149], 'DisplayName', 'A2 Path Taken')
axis equal

function [ A, l ] = update_lin_constrs(A, l, M, N, nx, nu, min_dists, x)
   
%     min_dist = 0.3;
%     [ m, n ] = size(A);
    start_idx = 136;
    % WORK OUT FORMULA FOR
    
    
    n_states = (N+1)*nx + N*nu;
    
%     Aineq = blkdiag(speye( (N+1)*nx + N*nu ), zeros(2));
    % reset last N+1 constraints
%     A(end-N:end,:) = zeros(N+1, n);

    %NOT TRUE should be M Comb 2 
    
    
    l(end-(N+1)*nchoosek(M,2)+1:end) = min_dists;
    
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

                    T = eta' * [ eye(2), -eye(2) ];
                    
                    A(start_idx+k+(idx_i+idx_j-3)*(N+1),...
                        [ ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 2,...
                        ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 2]) = T;
                    
                end
            end
        end
    end
end







