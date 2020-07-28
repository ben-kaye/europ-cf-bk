function cost = get_fitness(Q_weights, R_weights)
    step_info = get_step_info(Q_weights, R_weights);
    
    %implement formula
    cost = 0; %linear comb of step_info

end

function step_info = get_step_info(Q_weights,R_weights)

LQR_K = compute_lqr(Q_weights, R_weights);

mdl = 'quadcop_simulator';
mdlwksp = get_param(mdl,'ModelWorkspace');

Ktemp = getVariable(mdlwksp,'LQR_K');
Ktemp = LQR_K;
assignin(mdlwksp,'LQR_K',Ktemp)

sim_out = sim(mdl, 'SimulationMode','normal','AbsTol','1e-5',...
                     'SaveOutput','on');
                           
y = sim_out.yout{1}.Values.Data;
t = sim_out.tout;

start_times = [5 10 0 15];

step_info = {[],[],[],[]};


for k = 1:4    
    idx = find(t>=start_times(k),1);
    t_new = t(idx:end) - t(idx);
    y_trim = y(idx:end, k);
    step_info{k} = stepinfo(y_trim,t_new);    
end

end

function Kd = compute_lqr(Qweights, Rweights)

Q = diag(Qweights);
R = diag(Rweights);

g = 9.81;
m = 27e-3; % should use model params.. oh well
T = 1/200;

A = [   zeros(3),   eye(3),     zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0];
        zeros(3),   zeros(3),   zeros(3) ];
    
    
B =     [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        [ zeros(3,1), eye(3)] ];    
    
% State variables:
% x, y, z, xdot, ydot, zdot, gamma, beta, alpha,
% { position }, { velocity }, { euler orientation }
    
[ Kd, ~, ~] = lqrd(A, B, Q, R, 0, T);

Kd(abs(Kd)<1e-10) = 0;

end
