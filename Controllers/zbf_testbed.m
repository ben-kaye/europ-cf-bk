opts =  optimset('Display','off');

SIM_TIME = 10; % {s}
STEP_SIZE = 1e-3; % {s}
Ts = 1e-1; % {s}



x = [ 3; 2; -3*pi/18 ]; % [ p_x, p_y, phi ]
r = [ 2; 3; pi/3; 0.5; 1; 0]; % [ r_x, r_y, phi_r, phi_rdot, v_r, v_rdot ]
p_o = [ -1; 6 ]; % {m, m}

DELTA = 1.3; % {m} %%% BREAKS AT 1.3
MIN_V = 0; % {ms-1}
MAX_V = 3; % {ms-1}
MAX_TURN = 1.5; % {rads-1}


gamma = 3;
k1 = 5;
k2 = 5;


cbf_params.H = eye(2);
cbf_params.ctrl_min = [ MIN_V; -MAX_TURN ];
cbf_params.ctrl_max = [ MAX_V; MAX_TURN ];
cbf_params.delta = DELTA;
cbf_params.gamma = gamma;
cbf_params.k1 = k1;
cbf_params.k2 = k2;
cbf_params.options = opts;

N_SAMPLE = floor(SIM_TIME/Ts);
N_SUBS = floor(Ts/STEP_SIZE);

Ns = N_SAMPLE;
delta = DELTA;
x_t = NaN*ones(2, N_SAMPLE);
r_t = NaN*ones(2, N_SAMPLE);

u = [ MAX_V; MAX_TURN ];
errs = logical(zeros(1, N_SAMPLE));
for e = 1:N_SAMPLE
    u = cbf_qp_controller(x, r, p_o, u(1), cbf_params);
    for s = 1:N_SUBS
        [x, r] = sim_xr(x, u, r, STEP_SIZE);
    end
    x_t(:,e) = x([1,2]);
    r_t(:,e) = r([1,2]);    
end

plot_res;


function [x, r] = sim_xr(x, u, r, step_sz)
    phi = x(3);
    x_dot = [ cos(phi), 0; sin(phi), 0; 0, 1 ] * u;
    
    phi_r = r(3);
    v_r = r(5);
    phi_rdot = r(4);
    v_rdot = r(6);
    r_dot =  [ v_r*[ cos(phi_r); sin(phi_r) ]; phi_rdot; 0; v_rdot; 0 ];
    
    x = x + step_sz*x_dot;
    
   
    r = r + step_sz*r_dot;    
end

