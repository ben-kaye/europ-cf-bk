% clear
% load('result.mat')

figure(1)
Nmarkers = 20;
path = plot(x_t(1,:), x_t(2,:), '-', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:ceil(Ns/Nmarkers):Ns, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
hold on
plot(x_t(1,errs), x_t(2,errs), 'o', 'MarkerSize', 10, 'Color', [1, 0, 0], 'DisplayName', 'Infeasiblity')
plot(r_t(1,:), r_t(2,:), '-.', 'LineWidth',1, 'Marker','^', 'MarkerIndices',1:ceil(Ns/Nmarkers):Ns, 'Color',[0.3, 0.3, 0.3],'DisplayName', 'Reference');
viscircles(p_o', delta, 'LineStyle', '-.');
legend()
plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1,0,0],'DisplayName','Object')
hold off
xlabel('x (m)')
ylabel('y (m)')
axis equal
drawnow;
