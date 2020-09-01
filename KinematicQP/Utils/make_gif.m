fig = figure('units','normalized','outerposition',[0 0 1 1]);

maxy = max( max(x_t(2,:)), p_o(2) + delta );
axis_lims = [ min(x_t(1,:)), max(x_t(1,:)),min(x_t(2,:)), maxy ];

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])


obj_path = plot(p(1,1), p(2,1),'.','MarkerSize',20, 'Color',1/255*[60, 71, 163]);
hold on
ref_path = plot(r_t(1,1),r_t(2,1), 'k-.', 'LineWidth',1);
plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1,0,0],'DisplayName','Object')
viscircles(p_o', delta, 'LineStyle', '-.');
hold off


for it = 2:Ns
    p = x_t(:,it);
    
    delete(obj_path);
    delete(ref_path);

    hold on
    obj_path = plot(p(1), p(2),'.','MarkerSize',20, 'Color',1/255*[60, 71, 163]);
    ref_path = plot(r_t(1,1:it),r_t(2,1:it), 'k-.', 'LineWidth',1);
    hold off 
    
    axis(axis_lims)
    axis equal
       
    drawnow
    
    file_path = sprintf('_frames\\zbf_frame_%d',it);
    print('-r100',file_path,'-dpng')

end

close(fig)