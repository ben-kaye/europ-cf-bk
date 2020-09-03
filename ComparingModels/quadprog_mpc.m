clear,clc

% PARAMETERS
VIS_ON = 1;
Np = 10;
MIN_DIST = 0.05;
IDEAL_DIST = 0.9;
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 30;
STEP_SIZE = 1e-3;
Ts = 5e-2;

delta = IDEAL_DIST;

% initial conditions
ref = [ 2; 0.5; 0; 0; 0; 0 ];
x0 = [ 0.1; 0.1; 0; 0; 0; 0 ];
p_o = [ 1; 0.5 ];

xmpc = x0;
xbf = x0;

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
    
Q = diag([ 180, 180, 28, 15, 15, 7, 60, 60, 10, 5, 5, 3 ]);
R = diag([ 0.10, 0.10, 0.10, 0.10 ]);

g = 9.81;
m = 27e-3;
Jinv = diag([ 1/16, 1/16, 1/30 ]);


A = [   zeros(3),   eye(3),     zeros(3),                       zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0],   zeros(3);
        zeros(3),   zeros(3),   zeros(3),                       eye(3);
        zeros(3),   zeros(3),   zeros(3),                       zeros(3)    ];
    
B = [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        zeros(3,4);
        [ zeros(3,1), Jinv ] ];
    

[ KLQRd, ~, ~ ] = lqrd(A, B, Q, R, Ts);

[nx, nu] = size(Bd);

Ax = kron(eye(Np + 1), -eye(nx)) + kron(diag(ones(Np, 1), -1), Ad);
Bx = kron([ zeros(1, Np); eye(Np) ], Bd);
Aeq_mpc = [ Ax, Bx ];
beq_mpc = [ -x0; zeros(Np*(nx + nu), 1) ];

Qmpc = diag([1 1 0.1 0.1 0 0]);
Rmpc = 1e-1 * eye(nu);
Hmpc = blkdiag( kron(eye(Np + 1), Qmpc), kron(eye(Np), Rmpc) );
fmpc = [ repmat( -Qmpc*ref, Np + 1, 1 ); zeros(Np*nu, 1) ];

dotp_ref = zeros(2,1);


[ Ampc, bmpc ] = get_MPC_constr(beq_mpc, Np, nx, nu, p_o, delta);
UB_mpc = [ inf*ones(nx*(Np + 1), 1); max_v*ones(Np*nu, 1) ];
LB_mpc = -UB_mpc;

% (H,f,A,b,Aeq,beq,LB,UB)
mpc_u = quadprog(Hmpc, fmpc, Ampc, bmpc, Aeq_mpc, beq_mpc, LB_mpc, UB_mpc);


clf = clf_control(xmpc, ref, k1, k2);
clf = sat_ctrls(clf, ctrl_min, ctrl_max);
f = -Hmpc*clf;
Hbf = eye(2);




ctrl = get_u(xmpc, dotp_ref, m, g, KLQRd);

N = floor(SIM_TIME/Ts);
for s = 1:N
    % solve
    mpc_u = quadprog(Hmpc, fmpc, Ampc, bmpc, Aeq_mpc, beq_mpc, LB_mpc, UB_mpc);
    bf_u = quadprog(Hbf, fbf, Abf, bbf, [], [], [], []);
    
    %truncate bf_u
    %truncate mpc_u
    
    ctrl_mpc = mpc_u((Np + 1)*nx + 1:(Np + 1)*nx + nu);    
    bf_u = sat_ctrls(bf_u, ctrl_min, ctrl_max);
    
    % need to track angle some how
    phi
    
    %phi = phi_last + omega*time_step;
    
    ctrl_bf = bf_u * [ cos(phi); sin(phi) ];
    
    
    % simulate
    ctrl_a = get_u(xmpc, ctrl_mpc, m, g, KLQRd);
    ctrl_b = get_u(xbf, ctrl_bf, m, g, KLQRd);
    
    xmpc = simulate(xmpc, ctrl_a, time_step);
    xbf = simulate(xbf, ctrl_b, time_step);
    
    
    % update mpc
    beq_mpc(1:nx) = xmpc;
    [ Ampc, bmpc ] = get_MPC_constr(beq_mpc, Np, nx, nu, p_o, delta);  
    
    % update bf
    
    clf = clf_control(xmpc, ref, k1, k2);
    clf = sat_ctrls(clf, ctrl_min, ctrl_max);
    f = -Hmpc*clf;   
end

function u = get_u(x, dotp, m, g, KLQRd)
    u = [ m*g; 0; 0; 0 ] - KLQRd*(x - [ zeros(3,1); dotp; 0; zeros(6,1) ]);
end

function x = simulate(x, u, time_step)
    u
end

function [ Aub, bub ] = get_MPC_constr(x_mpc, Np, nx, nu, p_o, delta)
    x_mpc = x_mpc(1:(Np+1)*nx);
    x_mpc = reshape(x_mpc, nx, []);
    x_k = x_mpc([1,2], :);
    
    p_ox = x_k(:,2:end) - p_o;
    eta = p_ox./vecnorm(p_ox);
    
    bub = -delta*ones(Np, 1) - eta'*p_o;
    
    Aub = [ zeros(Np, nx), kron(eye(Np),[1, 1, zeros(1, nx - 2)]), zeros(Np, nu*Np) ];
    Aub(Aub~=0) = -eta;    
end
