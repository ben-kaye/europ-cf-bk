
opts =  optimset('Display','off');

% PARAMETERS
VIS_ON = 1;
Np = 10;
MIN_DIST = 0.05;
IDEAL_DIST = 1.5;
EPSILON = [ 0.1 0.5 1 ];
SIM_TIME = 10;
STEP_SIZE = 1e-3;
Ts = 5e-2;

delta = IDEAL_DIST;
max_v = 3;
min_v = 0;
max_omeg = 1.5;
gamma = 0.3;
k1 = 1;
k2 = 1;


init_pos = [ -2; 1 ];
end_pos = [ 15; 9 ];
p_o = [ 2; 3 ];

dpath = end_pos - init_pos;
r_phi = atan2(dpath(2), dpath(1));

r_cbf = [ init_pos; r_phi; 0; 2/3*max_v; 0];
r_mpc = [ end_pos; zeros(4, 1) ];

x_a = [ init_pos; zeros(10, 1) ];
x_b = [ init_pos; zeros(10, 1) ];

x_mpc = [ init_pos; zeros(4, 1) ];
x_cbf = [ init_pos; r_phi ];

Np = 10; % predictive horizon;

g = 9.81;
m = 0.1;
Jinv = diag([ 1/20, 1/20, 1/35 ]);
J = diag([20, 20, 35]);


A_lqr = [   zeros(3),   eye(3),     zeros(3),                       zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0],   zeros(3);
        zeros(3),   zeros(3),   zeros(3),                       eye(3);
        zeros(3),   zeros(3),   zeros(3),                       zeros(3)    ];
    
B_lqr = [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        zeros(3,4);
        [ zeros(3,1), Jinv ] ];
    
        
Q_lqr = diag([ 180, 180, 28, 15, 15, 7, 60, 60, 10, 5, 5, 3 ]);
R_lqr = diag([ 0.10, 0.10, 0.10, 0.10 ]);

[ KLQRd, ~, ~ ] = lqrd(A_lqr, B_lqr, Q_lqr, R_lqr, Ts);
KLQRd(abs(KLQRd) < 1e-9) = 0;
Ku = KLQRd(:, 3:end);


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

mpc_params.H = blkdiag( kron(eye(Np + 1), Q_mpc), kron(eye(Np), R_mpc) );
mpc_params.f = [ repmat( -Q_mpc*r_mpc, Np + 1, 1 ); zeros(Np*nu, 1) ];

mpc_params.Aeq = [ Ax, Bx ];
mpc_params.beq = [ -x_mpc; zeros(Np*nx, 1) ];

mpc_params.UB = [ inf*ones(nx*(Np + 1), 1); max_v*ones(Np*nu, 1) ];
mpc_params.LB = -mpc_params.UB;

mpc_params.nu = nu;
mpc_params.nx = nx;
mpc_params.Np = Np;
mpc_params.options = opts;

cbf_params.H = eye(2);
cbf_params.ctrl_min = [ min_v; -max_omeg ];
cbf_params.ctrl_max = [ max_v; max_omeg ];
cbf_params.delta = delta;
cbf_params.gamma = gamma;
cbf_params.k1 = k1;
cbf_params.k2 = k2;
cbf_params.options = opts;

N_SAMPLE = floor(SIM_TIME/Ts);
N_SUBS = floor(Ts/STEP_SIZE);


mpc_err = 0;
% init mpc constraints
[ ctrl_mpc, u_last, mpc_err ] = mpc_qp_controller(x_mpc, zeros((1+Np)*nx, 1), p_o, -inf, mpc_err, mpc_params);
for dist_coeff = EPSILON
    [ ctrl_mpc, u_last, mpc_err ] = mpc_qp_controller(x_mpc, u_last, p_o, dist_coeff*delta, mpc_err, mpc_params);
end

for sp = 1:N_SAMPLE
    x_cbf = [ x_a([1,2]); atan2(x_a(5), x_a(4)) ];
    v_last = vecnorm(x_a([4,5]));
    ctrl_cbf = cbf_qp_controller(x_cbf, r_cbf, p_o, v_last, cbf_params);
    phi = x_cbf(3) + ctrl_cbf(2)*Ts;
    ctrl_a = get_u(x_a, ctrl_cbf(1)*[ cos(phi); sin(phi) ], m, g, Ku);
    
    x_mpc = x_b([1,2,4,5,7,8]);
    [ ctrl_mpc, u_last, mpc_err ] = mpc_qp_controller(x_mpc, u_last, p_o, delta, mpc_err, mpc_params);
    ctrl_b = get_u(x_b, ctrl_mpc, m, g, Ku);
    
    for ss = 1:N_SUBS
        % simulate
        x_a = simulate(x_a, ctrl_a, STEP_SIZE, J, Jinv, m, g);
        x_b = simulate(x_b, ctrl_a, STEP_SIZE, J, Jinv, m, g);
        r_cbf = simulate_ref(r, STEP_SIZE);
    end
end


function u = get_u(x, dotp, m, g, Ku)
    u = [ m*g; 0; 0; 0 ] - Ku*(x(3:end) - [ 0; dotp; zeros(7, 1) ]);
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