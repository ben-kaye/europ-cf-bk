% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *                         Ben Kaye, (c) 2020                          *
% *                           EUROP Project:                            *
% *                               MPC-CA                                *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                                     *
% *             CrazyFlie controller model by Aren Karapet              *
% *                                                                     *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear,clc

opts =  optimset('Display','off');

%%% PARAMETERS %%%

VIS_ON = 1;
Np = 5; % predictive horizon
MIN_DIST = 0.05; % {m} 
IDEAL_DIST = 0.7; % {m} radius of avoidance
EPSILON = [ 0.1 0.5 1 ]; % init sequence ( dist = epsilon * delta )
SIM_TIME = 30; % {s}

%%% INITIAL CONDITIONS %%% 

ref = [ 2; 0.5; 0; 0; 0; 0 ];
x0 = [ 0.1; 0.1; 0; 0; 0; 0 ];
p_o = [ 1; 0.5 ];

% Hard coded state transition. Obtained in formulate_state_space.m
Ad = [   1,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
        0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
        0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
        0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
        0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
        0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   ];

Bd = [   0.00370875152761323,    0;
    0,                      0.00370875152761323;
    0.105706677250003,      0;
    0,                      0.105706677250003;
    0,                      -0.193245482770890;
    0.193245482770890,      0                       ];
    
[nx, nu] = size(Bd);     

%%% CONSTRAINTS %%%

umin = [ -1; -1 ];
umax = [ 1; 1 ];
xmin = [ -inf; -inf; -2; -2; -inf; -inf ];
xmax = [ inf; inf; +2; +2; inf; inf ];

%%% OBJECTIVES %%%
Q = diag([1 1 0.1 0.1 0 0]);
QN = Q;
R = 0.01 * eye(nu);

% quadratic objective
P = blkdiag( kron(speye(Np), Q ), QN, kron(speye(Np), R));
% linear objective
q = [ repmat( -Q*ref, Np, 1 ); -QN*ref; zeros(Np*nu, 1)];


% constrain predicted system to dynamic model
Ax = kron( speye(Np+1), -speye(nx) ) +...
    kron( sparse(diag(ones(Np,1), -1)), Ad);
Bu = kron( [ sparse(1, Np); speye(Np) ], Bd );

Aeq = [ Ax, Bu ];

leq = [ -x0; zeros(Np*nx, 1) ];
ueq = leq;

% input and state constraints

Aineq = speye((Np+1)*nx + Np*nu);

lineq = [ repmat( xmin, Np+1, 1); repmat( umin, Np, 1)  ];
uineq = [ repmat( xmax, Np+1, 1); repmat( umax, Np, 1) ];

[ A_dist, ~ ] = dist_constraints(Np, nx, nu, p_o, zeros(Np, 1), -leq);
l_dist = -inf*ones(Np, 1);
u_dist = inf*ones(Np,1);

% simulate

if (VIS_ON)
    figure(1)
    plot(x0(1), x0(2), 'b+', 'MarkerSize', 30, 'DisplayName', 'Start')
    hold on
    plot(ref(1), ref(2), 'rx', 'MarkerSize', 30, 'DisplayName', 'End')
    plot(p_o(1), p_o(2), 'r+', 'DisplayName', 'Radius of Avoidance')
    viscircles(p_o', IDEAL_DIST, 'LineStyle', '-.')

    xlabel('x (m)')
    ylabel('y (m)')
    title('Drone at time kT, T=0.1s')
end

error_count = 0;

x_hist = NaN*ones((Np+1)*nx, SIM_TIME);


% initialisation sequence for dist constraints
xN = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);
for e = EPSILON
    
    relaxation_params = e * ones(Np, 1);
    
    dists = MIN_DIST * ( ones(Np,1) - relaxation_params ) + IDEAL_DIST * relaxation_params;
    
    [ A_dist, u_dist ] = dist_constraints(Np, nx, nu, p_o, dists, xN);

    xN = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);  
end

%%% SIMULATE %%%

for i = 1 : SIM_TIME
    
    x_result = quadprog(P, q, A_dist, u_dist, Aeq, ueq, lineq, uineq, [], opts);
    
    if sum(isnan(x_result)) > 0
        % err
        error_count = error_count + 1;
    else
        % solved
        xN = x_result;
    end
    
    % get control actions
    ctrls = xN( (Np+1)*nx + 1 : end );
    ctrls = reshape(ctrls, nu, []);
    
    % in event of unsolvable path, carry on with previously calculated
    % optimal path
    ctrl = ctrls(:, 1 + error_count);
    
    % activate control
    x0 = Ad * x0 + Bd * ctrl;
    
    % store state
    x_hist(1:nx,i) = x0;
    
    % compute linear constraints
    [ A_dist, u_dist ] = dist_constraints(Np, nx, nu, p_o, dists, xN);
        
    leq(1:nx) = -x0;
    ueq(1:nx) = -x0; 
end

% plot results
if (VIS_ON)
    plot(x_hist(1,:), x_hist(2,:), '-^', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', 1/255*[73, 146, 214], 'DisplayName', 'Path Taken')
    axis equal
    legend();
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