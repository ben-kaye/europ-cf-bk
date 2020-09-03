clear, clc
% Goal here is to test the Barrier Function QP vs the MPC QP in a unified model, which neither is based on. 
% We will use the full drone dynamics with a intermediate LQR controller.


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

init_pos = [ -2; 1 ];
x0mpc = [ init_pos; 0; 0; 0; 0 ];
x0bf = [ init_pos; -30*pi/180 ];
p_o = [ 1; 0.5 ];

xmpc = x0mpc;
xbf = x0bf;

% System State-Space

% d_system.A([1:2, 4:5, 7:8], [1:2, 4:5, 7:8])

% { x, y, xdot, ydot, roll, pitch }
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

% initiliase MPC constraints
for z = EPSILON
    mpc_u = quadprog(Hmpc, fmpc, Ampc, bmpc, Aeq_mpc, beq_mpc, LB_mpc, UB_mpc);
    [ Ampc, bmpc ] = get_MPC_constr(mpc_u, Np, nx, nu, p_o, z*delta);
end

clf = clf_control(xmpc, ref, k1, k2);
clf = sat_ctrls(clf, ctrl_min, ctrl_max);
fbf = -Hbf*clf;
Hbf = eye(2);



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
    phi = xbf(3) + bf_u(2)*Ts; 

    %%% NOTE USING Ts NOT STEP_SIZE    
    %phi = phi_last + omega*time_step;
    
    ctrl_bf = bf_u * [ cos(phi); sin(phi) ];
    
    
    % simulate
    ctrl_a = get_u(x_a, ctrl_mpc, m, g, KLQRd);
    ctrl_b = get_u(x_b, ctrl_bf, m, g, KLQRd);
    
    x_a = simulate_drone(x_a, ctrl_a, STEP_SIZE);
    x_b = simulate_drone(x_b, ctrl_b, STEP_SIZE);
    ref = simulate_ref(ref, STEP_SIZE);
    
    xmpc = x_a([1,2,4,5,7,8]);    
    xbf = [ x_b([1,2]); atan2(x_b(5),x_b(4) ];

    % update mpc
    beq_mpc(1:nx) = xmpc;
    [ Ampc, bmpc ] = get_MPC_constr(beq_mpc, Np, nx, nu, p_o, delta);  
    
    % update bf
    
    clf = clf_control(xbf, ref, k1, k2);
    clf = sat_ctrls(clf, ctrl_min, ctrl_max);
    fbf = -Hbf*clf;   
end

function u = get_u(x, dotp, m, g, KLQRd)
    u = [ m*g; 0; 0; 0 ] - KLQRd*(x - [ zeros(3,1); dotp; 0; zeros(6,1) ]);
end

function x = simulate_drone(x, thrust, step_sz, J, Jinv, m, g)
    ftotal = thrust(1);


    pos_dot = x(4:6);
    euler = x(7:9);
    euler_dot = x(10:12);

    roll = euler(1);
    pitch = euler(2);
    yaw = euler(3);
    roll_dot = euler_dot(1);
    pitch_dot = euler_dot(2);
    % yaw_dot = euler_dot(3);

    % roll, pitch, yaw

    iRb = [ cos(pitch)*cos(yaw),   -cos(roll)*sin(yaw) + sin(roll)*sin(pitch)*cos(yaw),    sin(roll)*sin(yaw) + cos(roll)*sin(pitch)*cos(yaw);
            cos(pitch)*sin(yaw),   cos(roll)*cos(yaw) + sin(roll)*sin(pitch)*sin(yaw),     -sin(roll)*cos(yaw) + cos(roll)*sin(pitch)*sin(yaw);
            -sin(pitch),           cos(pitch)*sin(roll),                                   cos(pitch)*cos(roll)];

    T = [ 1,    0,          -sin(pitch);
        0,    cos(roll),   sin(roll)*cos(pitch);
        0,    -sin(roll),  cos(roll)*cos(pitch)];

    Tinv = [    1,     sin(roll)*tan(pitch),    cos(roll)*tan(pitch);
                0,     cos(roll),               -sin(roll);
                0,     sin(roll)*sec(pitch),    cos(roll)*sec(pitch)];

    Tdot = [    0,  0,                      -cos(pitch) * pitch_dot;
                0,  -sin(roll) * roll_dot,    -sin(roll)*sin(pitch)*pitch_dot + cos(roll)*cos(pitch)*roll_dot;
                0   -cos(roll) * roll_dot,    -cos(roll)*sin(pitch)*pitch_dot - sin(roll)*cos(pitch) * roll_dot];

    omega = T * euler_dot;
    omega_dot = Jinv * ( thrust(2:4) - cross(omega, J*omega) );
    euler_dotdot = Tinv * ( omega_dot - Tdot * euler_dot );

    pos_dotdot = 1/m * ( iRb * [0; 0; ftotal] - [0; 0; m*g]);

    state_dot = [   pos_dot;
                    pos_dotdot;
                    euler_dot; 
                    euler_dotdot ];
                
    x = x + state_dot*step_sz;
end

function r = simulate_ref(r, step_sz)
    phi_r = r(3);
    v_r = r(5);
    phi_rdot = r(4);
    v_rdot = r(6);
    r_dot =  [ v_r*[ cos(phi_r); sin(phi_r) ]; phi_rdot; 0; v_rdot; 0 ];
    
    r = r + step_sz*r_dot;
end

function [ Abf, bbf ] = get_BF_constr(x_bf, v_last, p_o, delta, gamma)
    Abf = 0;
    bbf = 0;
    
    p = x_bf([1,2]);
    phi = x_bf(3);
    
    z = p - p_o;
    
    zz = z'*z;
    
    c = cos(phi);
    s = sin(phi);
    
    v_dir = [ c; s ];
    
    b = v_dir'*z;    
    
    sin2 = 2*c*s;
    cos2 = c^2 - s^2;   
    
    h = 2*zz - b^2 - 2*delta^2;
    
    omeg = (z(1)^2-z(2)^2)*sin2 - 2*z(1)*z(2)*cos2;
    
    Lfh = 2*v_last*z'*v_dir;
    Lgh = [ 0, omeg ];

    bbf = gamma*h + Lfh;
    Abf = -Lgh;    
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
