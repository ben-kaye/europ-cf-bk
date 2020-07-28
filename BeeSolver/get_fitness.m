function cost = get_fitness(Q_weights, R_weights)
    step_info = get_step_info(Q_weights, R_weights);
    
    Z = zeros(5,4);
    
    for h = 1:4
        step = step_info{h};
        Z(1, h) = step.RiseTime;
        Z(2, h) = step.SettlingTime;
        Z(3, h) = step.PeakTime;
        Z(4, h) = step.Overshoot;
        Z(5, h) = step.SSError;
    end
    
    weights = [ 5 5 10 1 ];
    
    M = [ 30 5 5 1000 100 ];
    
    
    %implement formula
    cost = weights*(M*Z)';
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
    step_info{k}.SSError = abs(1 - y_trim(end));
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
