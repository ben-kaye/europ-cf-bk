t = 0:step_size:sim_time;
z = e*N/Ns;

t = t(1:e);
h_t = h_t(1:e);
x_t = x_t(:,1:e);
x_rt = x_rt(:,1:e);

u_t = x_t([4,5],:);

u_t(:,1) = [nan;nan];

figure(2)
plot(t,u_t,'-','LineWidth', 1.5)
hold on 
plot(t,h_t)
legend('v','\omega')

figure(1)
px = x_t(1,:);
py = x_t(2,:);
plot(px,py,'-','LineWidth',1.5, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
hold on
plot(x_rt(1,:), x_rt(2,:), '-.', 'LineWidth',1, 'Color',[0.3, 0.3, 0.3],'DisplayName', 'Reference');
viscircles(p_o', delta);
axis equal