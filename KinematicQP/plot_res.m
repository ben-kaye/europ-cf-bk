% clear
% load('result.mat')

figure(1)

plot(x_t(1,:), x_t(2,:), '-', 'LineWidth',1.5, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
hold on
plot(r_t(1,:), r_t(2,:), '-.', 'LineWidth',1, 'Color',[0.3, 0.3, 0.3],'DisplayName', 'Reference');
viscircles(p_o', delta, 'LineStyle', '-.');
axis equal