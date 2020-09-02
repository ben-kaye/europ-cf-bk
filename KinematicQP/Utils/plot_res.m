% clear
% load('result.mat')
indices = 1:Ns;
err_indices = indices(errs);
start_err = min(err_indices);
end_err = max(err_indices);

safe_p = x_t;
safe_p(:,start_err:end_err) = NaN;
err_p = x_t(:,max(1,(start_err-1)):(end_err+1));


figure(1)
Nmarkers = 20;
viscircles(p_o', delta, 'LineStyle', '-.', 'Color', [1, 0.6, 0]);

hold on
path = plot(safe_p(1,:), safe_p(2,:), '-', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:ceil(Ns/Nmarkers):Ns, 'Color',1/255*[64, 201, 255], 'DisplayName', 'Position');
plot(err_p(1,:), err_p(2,:), '-', 'LineWidth', 1.5, 'Color', [1, 0, 0], 'DisplayName', 'Infeasiblity')
plot(r_t(1,:), r_t(2,:), '-.', 'LineWidth',1, 'Marker','^', 'MarkerIndices',1:ceil(Ns/Nmarkers):Ns, 'Color',[0.3, 0.3, 0.3],'DisplayName', 'Reference');
legend()
plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1, 0.6, 0],'DisplayName','Object')
hold off
xlabel('x (m)')
ylabel('y (m)')
axis equal

if exist('plot_title') == 1
    title(plot_title);
end
drawnow;
