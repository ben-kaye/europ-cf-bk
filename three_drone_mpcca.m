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
N = 5; % Horizon
M = 3; % Number of drones
MIN_DIST = 0.05; % Absolute min
IDEAL_DIST = 0.3; % Ideal range
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 50; % Number of iterations in simulation (linear sim)

% initial conditions
ref = [ 1,  1,  -1;
        -1, 1,  0;
        0,  0,  0;
        0,  0,  0;
        0,  0,  0;
        0,  0,  0 ];
xinit = [  -1, -1, 1;
        1,  -1, 0;
        0,  0,  0;
        0,  0,  0;
        0,  0,  0;
        0,  0,  0   ];
    
if ( M == 2 )
    ref = ref(:,1:2);
    xinit = xinit(:,1:2);
end
    
x0 = xinit;

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

q = [];
for h = 1:M
q = [   q;
        repmat( -Q*ref(:,h), N, 1 ); -QN*ref(:,h); zeros(N*nu, 1) ];
end

% Augmented Dynamics
Ax = kron( speye(N+1), -speye(nx) ) +...
    kron( sparse(diag(ones(N,1), -1)), Ad);
Bu = kron( [ sparse(1, N); speye(N) ], Bd );

Aeq = [ Ax, Bu ];

leq = [];
for h = 1:M
    leq = [ leq;
            -x0(:,h); zeros(N*nx, 1) ];
end
ueq = leq;


lineq = [ repmat( xmin, N+1, 1); repmat( umin, N, 1)  ];
uineq = [ repmat( xmax, N+1, 1); repmat( umax, N, 1) ];


Aeq = kron(speye(M), Aeq);
lineq = repmat(lineq, M, 1);
uineq = repmat(uineq, M, 1);
P = kron( speye(M), P) ;


% input and state constraints

Aineq = speye( M * Nstates );

% Aconstr = [ kron(speye(N+1), [ 1 1 0 0 0 0]), -ones(N+1, N_OBJ_STATES), zeros(N+1, N*nu) ];

N_ca_constrs = (N+1)*nchoosek(M,2);

% osqp constraints

l = [ leq; lineq; -inf*ones(N_ca_constrs, 1)];
u = [ ueq; uineq; inf*ones(N_ca_constrs, 1)];

A = [ Aeq; Aineq; zeros(N_ca_constrs, size(Aeq,2)) ];
A = make_placeholder_A(A,M,N,nx,nu);

% osqp init and setup
prob = osqp;
prob.setup( P, q, A, l, u, 'warm_start', true, 'verbose', false );

% simulate



dists = IDEAL_DIST*ones((N+1)*nchoosek(M,2), 1);

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
    ctrls = zeros(N*nu, M);
    ctrl_act = zeros(nu, M);
    for g = 1:M
        ctrls(:,g) = xN( g*((N+1)*nx)+ (g-1)* N*nu+ 1: g*((N+1)*nx)+ g*N*nu);
        % if solve error, continue with previous
        ctrl_act(:, g) = ctrls( 1 + error_count*nu : (1+error_count)*nu, g);
    end
    
    % activate control
    x0 = Ad * x0 + Bd * ctrl_act;
    
    % store state
    x_hist(1:M*nx,i) = reshape(x0, M*nx, []);
    
    % compute linear constraints
    [ A, l ] = update_lin_constrs(A, l, M, N, nx, nu, dists, xN);
    
    % update initial conds and constraints
    
    xtemp = reshape(x0, M*nx, []);
    
    leq = [];
    for h = 1:M
        leq = [ leq;
                -x0(:,h); zeros(N*nx, 1) ];
    end
    
    l(1:M*(N+1)*nx) = leq;
    u(1:M*(N+1)*nx) = leq;
    
    prob.update('l',l, 'u',u);
    prob.update('Ax', nonzeros(A));
%   
  
end

% plot results

color = zeros(3, M);
color(:, 1) = 1/255*[73, 146, 214];
color(:, 2) = 1/255*[214, 88, 149];
color(:, 3) = 1/255*[108, 158, 70];

if (VIS_ON)
    figure(1)
    hold on
    for z = 1 : M        
        plot(xinit(1,z), xinit(2,z), '+', 'MarkerSize', 30, 'DisplayName', sprintf('A%d Start',z), 'Color', color(:,z));
        plot(ref(1,z), ref(2,z), 'x', 'MarkerSize', 30, 'DisplayName', sprintf('A%d End',z), 'Color', color(:,z));
        plot(x_hist(nx*(z-1)+1,:), x_hist(nx*(z-1)+2,:), '.', 'MarkerSize', 20, 'Color', color(:,z), 'DisplayName', sprintf('A%d Path Taken',z))
    end
end

axis equal
xlabel('x (m)')
ylabel('y (m)')
title('Drone at time kT, T=0.1s')
legend()


function [ A, l ] = update_lin_constrs(A, l, M, N, nx, nu, min_dists, x)
    n_states = (N+1)*nx + N*nu;
    start_idx = 2*M*(N+1)*nx + M*N*nu;
        
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

function [ A ] = make_placeholder_A(A, M, N, nx, nu) 
    start_idx = 2*M*(N+1)*nx + M*N*nu;
    for idx_i = 1:M
        for idx_j = 1:M
            if( idx_j > idx_i)                
                for k = 1:N+1
                    T = [ ones(1,2), -ones(1,2) ];
                    
                    A( start_idx+k+(idx_i+idx_j-3)*(N+1),...
                        [ ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_i-1) + (k-1)*nx + 2,...
                        ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 1 : ((N+1)*nx + N*nu)*(idx_j-1) + (k-1)*nx + 2]) = T;
                    
                end
            end
        end
    end
end





