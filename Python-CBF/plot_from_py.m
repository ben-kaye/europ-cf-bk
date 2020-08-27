clear
load('result.mat')

figure(1)
px = x_t(:,1);
py = x_t(:,2);
plot(px,py,'-','LineWidth',1.5, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
hold on
plot(r_t(:,1), r_t(:,2), '-.', 'LineWidth',1, 'Color',[0.3, 0.3, 0.3],'DisplayName', 'Reference');
viscircles(p_o, delta, 'LineStyle', '-.');
axis equal