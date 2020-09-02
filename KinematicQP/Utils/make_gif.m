FILE_NAME = 'H:\\Files\\EUROP-MATLAB\\KinematicQP\\_frames\\zbf_frame_%d';

fig = figure('units','normalized','outerposition',[0 0 1 1]);

maxy = max( max(x_t(2,:)), p_o(2) + delta );
miny = min( min(x_t(2,:)), p_o(2) - delta );
maxx = max( max(x_t(1,:)), p_o(1) + delta );
minx = min( min(x_t(1,:)), p_o(1) - delta );
axis_lims = [ minx, maxx, miny, maxy ];

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])

viscircles(p_o', delta, 'LineStyle', '-.','Color', [1, 0.6, 0]);
hold on
ref_path = plot(r_t(1,1),r_t(2,1), 'k-.', 'LineWidth',1,'DisplayName','Reference');
plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1, 0.6, 0],'DisplayName','Object')
obj_path = plot(p(1,1), p(2,1),'.','MarkerSize',20, 'Color',1/255*[60, 71, 163],'DisplayName','Drone');
hold off
legend()
axis(axis_lims)
axis equal

file_path = sprintf(FILE_NAME,1);
print('-r100',file_path,'-dpng')

for it = 2:Ns
    p = x_t(:,it);
    
    delete(obj_path);
    delete(ref_path);

    hold on
    obj_path = plot(p(1), p(2), '.', 'MarkerSize', 20, 'Color',1/255*[60, 71, 163], 'DisplayName', 'Drone');
    ref_path = plot(r_t(1,1:it), r_t(2,1:it), 'k-.', 'LineWidth', 1, 'DisplayName', 'Reference');
    hold off 
    legend()
    axis(axis_lims)
    axis equal
       
    drawnow
    
    file_path = sprintf(FILE_NAME,it);
    print('-r100',file_path,'-dpng')

end

close(fig)