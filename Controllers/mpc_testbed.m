clear,clc

opts =  optimset('Display','off');

%%% PARAMETERS %%%

VIS_ON = 1;

MIN_DIST = 0.05; % {m} 
IDEAL_DIST = 0.7; % {m} radius of avoidance
EPSILON = [ 0.1 0.5 1 ]; % init sequence ( dist = epsilon * delta )
SIM_TIME = 13; % {s}

STEP_SIZE = 1e-3;
MAX_V = 2;
DELTA = IDEAL_DIST;
delta = DELTA;

max_thrust = 1;

Ts = 1e-2;
Np = 30; % predictive horizon

% CONSTANTS

g = 9.81;
m = 0.1;
Jinv = diag([ 1/20, 1/20, 1/35 ]);
J = diag([20, 20, 35]);
layout = [ 1,       1,      1,      1;
           -0.06,  -0.06, 0.06,  0.06;
           -0.06,  0.06,  0.06,  -0.06;
           -0.006,  0.006,  -0.006, 0.006   ];
inv_layout = [  0.2500,   -4.1667,   -4.1667,  -41.6667;
                0.2500,   -4.1667,    4.1667,   41.6667;
                0.2500,    4.1667,    4.1667,  -41.6667;
                0.2500,    4.1667,   -4.1667,   41.6667 ];

%%% INITIAL CONDITIONS %%% 

r_mpc = [ 3; 1.5; 0; 0; 0; 0 ];
x_mpc = [ 0.1; 0.1; 0; 0; 0; 0 ];
p_o = [ 1, 2; 0.5, 0.7 ];
% p_o = [ 1; 0.5 ];

%%% SETUP LOW-LEVEL SIMULATION CONTROLLER %%%

A_lqr = [   zeros(3),   eye(3),     zeros(3),                       zeros(3);
            zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0],   zeros(3);
            zeros(3),   zeros(3),   zeros(3),                       eye(3);
            zeros(3),   zeros(3),   zeros(3),                       zeros(3)    ];
    
B_lqr = [   zeros(3,4);
        [   [ 0; 0; 1/m ],  zeros(3) ];
            zeros(3,4);
        [   zeros(3,1),     Jinv ] ];
    
        
Q_lqr = diag([ 150, 150, 28, 15, 15, 7, 60, 60, 10, 5, 5, 3 ]);
R_lqr = diag([ 0.10, 0.10, 0.10, 0.10 ]);

[ KLQRd, ~, ~ ] = lqrd(A_lqr, B_lqr, Q_lqr, R_lqr, Ts);
KLQRd(abs(KLQRd) < 1e-9) = 0;
Ku = KLQRd(:, 3:end);

% Hard coded state transition. Obtained in formulate_state_space.m
% discretised system of PID controllers + plant at Ts = 0.1

A_mpc = [   1,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
        0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
        0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
        0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
        0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
        0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   ];

B_mpc = [   0.00370875152761323,    0;
        0,                      0.00370875152761323;
        0.105706677250003,      0;
        0,                      0.105706677250003;
        0,                      -0.193245482770890;
        0.193245482770890,      0                       ];    
    
[nx, nu] = size(B_mpc);

Ax = kron(eye(Np + 1), -eye(nx)) + kron(diag(ones(Np, 1), -1), A_mpc);
Bx = kron([ zeros(1, Np); eye(Np) ], B_mpc);

Q_mpc = diag([1 1 0.1 0.1 0 0]);
R_mpc = 1e-1 * eye(nu);

params.H = blkdiag( kron(eye(Np + 1), Q_mpc), kron(eye(Np), R_mpc) );
params.f = [ repmat( -Q_mpc*r_mpc, Np + 1, 1 ); zeros(Np*nu, 1) ];

params.Aeq = [ Ax, Bx ];
params.beq = [ -x_mpc; zeros(Np*nx, 1) ];

params.UB = [ inf*ones(nx*(Np + 1), 1); MAX_V*ones(Np*nu, 1) ];
params.LB = -params.UB;

params.nu = nu;
params.nx = nx;
params.Np = Np;
params.options = opts;

x_hist = NaN*ones((Np+1)*nx, SIM_TIME);


% initialisation sequence for dist constraints
errc = 0;

[ ~, qp_res, errc ] = mpc_qp_controller_2(x_mpc, zeros((1+Np)*nx, 1), p_o, -inf, errc, params);
for dist_coeff = EPSILON
    [ ~, qp_res, errc ] = mpc_qp_controller_2(x_mpc, qp_res, p_o, dist_coeff*DELTA, errc, params);
end

%%% SIMULATE %%%


N_SAMPLE = floor(SIM_TIME/Ts);
N_SUBS = floor(Ts/STEP_SIZE);

state = [ x_mpc([1,2]); zeros(10, 1) ];
dotp = zeros(2, 1);

x_t = zeros(3, N_SAMPLE);
phi_t = zeros(3, N_SAMPLE);

errs = logical(zeros(1,N_SAMPLE));

for sp = 1:N_SAMPLE
    [ dotp, qp_res, errc ] = mpc_qp_controller_2(state([1,2,4,5,7,8]), qp_res, p_o, delta, errc, params);
    
    if errc > 0 
        errs(z) = 1;
    end
    
    ctrl = get_u(state, dotp, m, g, Ku);
    
    for ss = 1:N_SUBS
        state = simulate_drone(state, ctrl, STEP_SIZE, J, Jinv, layout, inv_layout, m, g, max_thrust);
    end
    
    x_t(:, sp) = state(1:3);
    phi_t(:, sp) = state(7:9);
end


figure(1)
viscircles(p_o', delta*ones(size(p_o, 2), 1), 'LineStyle', '-.');
hold on
plot(r_mpc(1), r_mpc(2), 'o');
plot(x_t(1,:), x_t(2,:))
axis equal
%%% FUNCTIONS %%%

function u = get_u(x, dotp, m, g, Ku)
    u = [ m*g; zeros(3, 1) ] - Ku*(x(3:end) - [ 0; dotp; zeros(7, 1) ]);
end

function x = simulate_drone(x, thrusts, step_sz, J, Jinv, layout, inv_layout, m, g, max_thrust)
    ftotal = thrusts(1);

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
            -sin(pitch),           cos(pitch)*sin(roll),                                   cos(pitch)*cos(roll) ];

    T = [ 1,    0,          -sin(pitch);
        0,    cos(roll),   sin(roll)*cos(pitch);
        0,    -sin(roll),  cos(roll)*cos(pitch)];

    Tinv = [    1,     sin(roll)*tan(pitch),    cos(roll)*tan(pitch);
                0,     cos(roll),               -sin(roll);
                0,     sin(roll)*sec(pitch),    cos(roll)*sec(pitch) ];

    Tdot = [    0,  0,                      -cos(pitch)*pitch_dot;
                0,  -sin(roll)*roll_dot,    -sin(roll)*sin(pitch)*pitch_dot + cos(roll)*cos(pitch)*roll_dot;
                0   -cos(roll)*roll_dot,    -cos(roll)*sin(pitch)*pitch_dot - sin(roll)*cos(pitch)*roll_dot ];

    omega = T * euler_dot;
    omega_dot = Jinv * ( thrusts(2:4) - cross(omega, J*omega) );
    euler_dotdot = Tinv * ( omega_dot - Tdot * euler_dot );

    pos_dotdot = 1/m * ( iRb * [0; 0; ftotal] - [0; 0; m*g]);

    state_dot = [   pos_dot;
                    pos_dotdot;
                    euler_dot; 
                    euler_dotdot ];
                
    x = x + state_dot*step_sz;
end