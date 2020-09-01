

fig = figure('units','normalized','outerposition',[0 0 1 1]);
axis_lims = [ -2, 4, 2, 7 ];

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])
for it = 1:Ns
    p = x_t(:,it);

    plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1,0,0],'DisplayName','Object')
    hold on
    plot(r_t(1,1:it),r_t(2,1:it), 'k-.', 'LineWidth',1)
    viscircles(p_o', delta, 'LineStyle', '-.');
    plot(p(1), p(2),'.','MarkerSize',20, 'Color',1/255*[60, 71, 163])
    
    hold off
    
    axis(axis_lims)
    axis equal
    drawnow
    
    path = sprintf('frames\\zbf_frame_%d',it);
    print('-r100',path,'-dpng')
%     pause(1e-2);
       
end
close(fig)