clear, clc
% Goal here is to test the Barrier Function QP vs the MPC QP in a unified model, which neither is based on. 
% We will use the full drone dynamics with a intermediate LQR controller.

opt =  optimset('Display','off');

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
min_v = 0.5;
max_omeg = 1.5;
gamma = 0.3;
k1 = 1;
k2 = 1;

ctrl_min = [ min_v; -max_omeg ];
ctrl_max = [ max_v; max_omeg ];

% initial conditions
end_pos = [ 15; 9 ];
init_pos = [ -2; 1 ];

r_mpc = [ end_pos; zeros(4, 1) ];

dpath = end_pos - init_pos;
r_phi = atan2(dpath(2), dpath(1));

r_bf = [ init_pos; r_phi; 0; 2/3*max_v; 0];

x0_mpc = [ init_pos; 0; 0; 0; 0 ];
x0_bf = [ init_pos; r_phi ];

x_a = [ init_pos; zeros(10,1) ];
x_b = x_a;

p_o = [ 2; 3 ];

x_mpc = x0_mpc;
x_bf = x0_bf;

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
J = diag([16, 16, 30]);


A = [   zeros(3),   eye(3),     zeros(3),                       zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0],   zeros(3);
        zeros(3),   zeros(3),   zeros(3),                       eye(3);
        zeros(3),   zeros(3),   zeros(3),                       zeros(3)    ];
    
B = [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        zeros(3,4);
        [ zeros(3,1), Jinv ] ];
    

[ KLQRd, ~, ~ ] = lqrd(A, B, Q, R, Ts);
KLQRd(abs(KLQRd) < 1e-9) = 0;
Ku = KLQRd(:, 3:end);
% remove the tracking of x,y coords


[nx, nu] = size(Bd);

Ax = kron(eye(Np + 1), -eye(nx)) + kron(diag(ones(Np, 1), -1), Ad);
Bx = kron([ zeros(1, Np); eye(Np) ], Bd);
Aeq_mpc = [ Ax, Bx ];
beq_mpc = [ -x0_mpc; zeros(Np*nx, 1) ];

Q_mpc = diag([1 1 0.1 0.1 0 0]);
R_mpc = 1e-1 * eye(nu);
H_mpc = blkdiag( kron(eye(Np + 1), Q_mpc), kron(eye(Np), R_mpc) );
f_mpc = [ repmat( -Q_mpc*r_mpc, Np + 1, 1 ); zeros(Np*nu, 1) ];

[ A_mpc, b_mpc ] = get_MPC_constr(-beq_mpc, Np, nx, nu, p_o, 0);
UB_mpc = [ inf*ones(nx*(Np + 1), 1); max_v*ones(Np*nu, 1) ];
LB_mpc = -UB_mpc;

% initiliase MPC constraints
mpc_u = quadprog(H_mpc, f_mpc, A_mpc, inf*ones(Np, 1), Aeq_mpc, beq_mpc, LB_mpc, UB_mpc, [], opt);
for z = EPSILON
    mpc_u = quadprog(H_mpc, f_mpc, A_mpc, b_mpc, Aeq_mpc, beq_mpc, LB_mpc, UB_mpc, [], opt);
    [ A_mpc, b_mpc ] = get_MPC_constr(mpc_u, Np, nx, nu, p_o, z*delta);
end

clf = clf_control(x_bf, r_bf, k1, k2);
clf = sat_ctrls(clf, ctrl_min, ctrl_max);
H_bf = eye(2);
f_bf = -H_bf*clf;
[ A_bf, b_bf ] = get_BF_constr(x_bf, clf(1), p_o, delta, gamma);



N_SAMPLE = floor(SIM_TIME/Ts);
N_SIM = floor(SIM_TIME/STEP_SIZE);

x_t = NaN*ones(4, N_SAMPLE);
r_t = NaN*ones(2, N_SAMPLE);
ctrls_t = NaN*ones(4, N_SAMPLE);
errc_mpc = 0;
errc_bf = 0;

switched = false;
for s = 1:N_SAMPLE
    % solve
    
    mpc_u_temp = quadprog(H_mpc, f_mpc, A_mpc, b_mpc, Aeq_mpc, beq_mpc, LB_mpc, UB_mpc, [], opt);
    bf_u_temp = quadprog(H_bf, f_bf, A_bf, b_bf, [], [], [], [], [],  opt);
    
    if isempty(mpc_u_temp)
        errc_mpc = errc_mpc + 1;
    else
        mpc_u = mpc_u_temp;
    end
    if isempty(bf_u_temp)
        errc_bf = errc_bf + 1;
    else
        bf_u = bf_u_temp;
    end
    
    %truncate bf_u
    %truncate mpc_u
    
    ctrl_mpc = mpc_u((Np + 1)*nx + 1:(Np + 1)*nx + nu);    
    bf_u = sat_ctrls(bf_u, ctrl_min, ctrl_max);
    
    % need to track angle some how
    phi = x_bf(3) + bf_u(2)*Ts; 

    %%% NOTE USING Ts NOT STEP_SIZE    
    %phi = phi_last + omega*time_step;
    
    ctrl_bf = bf_u(1) * [ cos(phi); sin(phi) ];
    
    
    % simulate
    ctrl_a = get_u(x_a, ctrl_mpc, m, g, Ku);
    ctrl_b = get_u(x_b, ctrl_bf, m, g, Ku);
    
    for sub_s = 1:floor(N_SIM/N_SAMPLE)
        x_a = simulate_drone(x_a, ctrl_a, STEP_SIZE, J, Jinv, m, g);
        x_b = simulate_drone(x_b, ctrl_b, STEP_SIZE, J, Jinv, m, g);
        r_bf = simulate_ref(r_bf, STEP_SIZE);
        
        dist_r = r_bf([1,2]) - end_pos;
        if dist_r'*dist_r < 5e-1 && ~switched
            ctrl_min = [ 0; -max_omeg ];
            r_bf(6) = -0.1;
            switched = true;
        end 
        if dist_r'*dist_r < 1e-2
            r_bf = [ end_pos; zeros(4,1) ];            
        end
    end
    
    if min(x_a(3), x_b(3)) < -5e-1
        error('Crash')
    end
    
    x_t(:, s) = [ x_a([1,2]); x_b([1,2]) ];
    r_t(:, s) = r_bf([1,2]);
    ctrls_t(:,s) = [ ctrl_mpc; ctrl_bf ];
    
    x_mpc = x_a([1,2,4,5,7,8]);    
    x_bf = [ x_b([1,2]); atan2(x_b(5), x_b(4)) ];

    % update mpc
    beq_mpc(1:nx) = -x_mpc;
    [ A_mpc, b_mpc ] = get_MPC_constr(mpc_u(1:(Np + 1)*nx), Np, nx, nu, p_o, delta);  
    
    % update bf
    
    clf = clf_control(x_bf, r_bf, k1, k2);
    clf = sat_ctrls(clf, ctrl_min, ctrl_max);
    f_bf = -H_bf*clf;
    [ A_bf, b_bf ] = get_BF_constr(x_bf, bf_u(1), p_o, delta, gamma);
end

fprintf('Simulation complete.\n');

figure(1);

idxs = 1:ceil(s/30):s;
viscircles(p_o', delta, 'LineStyle', '-.');
hold on
plot(p_o(1), p_o(2), 'r*', 'LineWidth', 1.5,'MarkerSize', 15, 'DisplayName', 'Object');
plot(init_pos(1), init_pos(2),'+', 'MarkerSize', 30, 'DisplayName', 'Start');
plot(end_pos(1), end_pos(2), 'x', 'MarkerSize', 30, 'DisplayName', 'Ref End');
plot(x_t(1,:), x_t(2,:), '-^', 'MarkerIndices', idxs,  'LineWidth', 1.5, 'DisplayName', 'MPC')
plot(x_t(3,:), x_t(4,:), '-^', 'LineWidth', 1.5, 'MarkerIndices', idxs, 'DisplayName', 'CBF-QP')
plot(r_t(1,:), r_t(2,:), 'k-.^', 'LineWidth', 1, 'MarkerIndices', idxs, 'DisplayName', 'CBF-ref')
axis equal
hold off
title('MPC and CBF-QP controller applied to unseen LQR Drone Model');
xlabel('x (m)')
ylabel('y (m)')
legend()

figure(2)
t = 0:Ts:SIM_TIME-Ts;
subplot(1,2,1)
plot(t, ctrls_t(1,:), 'LineWidth', 1.5, 'DisplayName', 'v_x')
hold on
plot(t, ctrls_t(2,:), 'LineWidth', 1.5,'DisplayName', 'v_y')
hold off
title('MPC Controls')
xlabel('t (s)')
legend();

subplot(1,2,2)
plot(t, ctrls_t(3,:), 'LineWidth', 1.5, 'DisplayName', 'v_x')
hold on
plot(t, ctrls_t(4,:), 'LineWidth', 1.5, 'DisplayName', 'v_y')
hold off
xlabel('t (s)')
title('CBF-QP Controls')
legend();
%%% FUNCTIONS %%%

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

function [ A_bf, b_bf ] = get_BF_constr(x_bf, v_last, p_o, delta, gamma)
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

    b_bf = gamma*h + Lfh;
    A_bf = -Lgh;    
end

function ctrls_bf = clf_control(x_bf, r, k1, k2)
    p = x_bf([1,2]);
    phi = x_bf(3);
    p_r = r([1,2]);
    phi_r = r(3);
    phi_rdot = r(4);
    v_r = r(5);
    v_rdot = r(6);
    
    e3 = phi_r - phi;
    e12 = [ cos(phi) sin(phi); -sin(phi) cos(phi) ] * (p_r - p);
    
    e1 = e12(1);
    e2 = e12(2);

    alpha = atan(e2/v_r);
    e3_aux = alpha + e3;

    omega = phi_rdot + 2*e2*v_r*cos(e3_aux/2-alpha) + (v_r^2*sin(e3) - e2*v_rdot)/(v_r^2 + e2^2) + k2*sin(e3_aux/2);
    v = v_r*cos(e3) - omega*v_r*sin(e3_aux/2)/(v_r^2 + e2^2) + k1*e1;
    
    ctrls_bf = [ v; omega ];
end

function [ A_mpc, b_mpc ] = get_MPC_constr(x_mpc, Np, nx, nu, p_o, delta)
    x_mpc = x_mpc(1:(Np+1)*nx);
    x_mpc = reshape(x_mpc, nx, []);
    x_k = x_mpc([1,2], :);
    
    p_ox = x_k(:,2:end) - p_o;
    eta = p_ox./vecnorm(p_ox);
    
    b_mpc = -delta*ones(Np, 1) - eta'*p_o;
    
    A_mpc = [ zeros(Np, nx), kron(eye(Np),[1, 1, zeros(1, nx - 2)]), zeros(Np, nu*Np) ];
    A_mpc(A_mpc~=0) = -eta;    
end

function ctrls_sat = sat_ctrls(ctrls, min_ctrl, max_ctrl)
    ctrls_sat = max(ctrls, min_ctrl);
    ctrls_sat = min(ctrls_sat, max_ctrl);
end

