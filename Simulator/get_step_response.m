clear, clc, clf

Q_weights = [100 100 400 1 1 0.111 3000 3000 50];
R_weights = [14.2 0.25 0.25 0.111];


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



[m, n] = size(y);

y_ss = y(end,:);

start_times = [5 10 0 15];
labels = {'x','y','z','yaw'};
colour = {1/255*[168, 50, 50], 1/255*[50, 168, 82], 1/255*[50, 103, 168], 1/255*[222, 177, 29]}; 
step_info = {[],[],[],[]};

step_response = {[],[],[],[]};

for k = 1:4
    
    idx = find(t>=start_times(k),1);
    t_new = t(idx:end) - t(idx);
    y_trim = y(idx:end, k);
    stopdex = find(t_new>=10, 1);
    
    step_response{k} = [ t_new(1:stopdex) y_trim(1:stopdex) ];
    
    hold on
    stairs(step_response{k}(:,1),step_response{k}(:,2), 'DisplayName',labels{k},'LineWidth',2, 'Color', colour{k});
    step_info{k} = stepinfo(y_trim,t_new);

end

qstr = [ 'Q=' sprintf('%.2f,', Q_weights) ];
rstr = [ 'R=' sprintf('%.2f,', R_weights) ];

title(sprintf('Step response for %s %s',qstr, rstr))
xlabel('t(s)')
% axis( [0 5 -0.1 1.2]);
legend();

% get velocity
% gradd = gradient(stepresponse{3}(:,2),stepresponse{3}(:,1));
% plot(stepresponse{3}(:,1),gradd);

