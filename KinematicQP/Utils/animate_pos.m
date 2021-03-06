fig = figure(3);

maxy = max( max(x_t(2,:)), p_o(2) + delta );
miny = min( min(x_t(2,:)), p_o(2) - delta );
maxx = max( max(x_t(1,:)), p_o(1) + delta );
minx = min( min(x_t(1,:)), p_o(1) - delta );
axis_lims = [ minx, maxx, miny, maxy ];

fps = 30;
Nf = ceil(Ns/sim_time/fps);

viscircles(p_o', delta, 'LineStyle', '-.', 'Color', [1, 0.6, 0]);
hold on
ref_path = plot(r_t(1,1), r_t(2,1), 'k-.', 'LineWidth', 1);
plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1, 0.6, 0],'DisplayName', 'Object')

path = plot(x_t(1,1), x_t(2,1), '.', 'MarkerSize', 20, 'Color', 1/255*[60, 71, 163]);
hold off

tic;

drawnow

for it = (1+Nf):Nf:Ns
    p = x_t(:,it);
    delete(path);
    delete(ref_path);

    hold on
    path = plot(p(1), p(2), '.', 'MarkerSize', 20, 'Color', 1/255*[60, 71, 163]);
    ref_path = plot(r_t(1,1:it), r_t(2,1:it), 'k-.', 'LineWidth',1);
    hold off 

    axis(axis_lims)
    axis equal
    drawnow
    
    frame_time = toc;
    
    if frame_time < 1/fps
        pause(1/fps - frame_time);
    end
    
    tic; 
end

draw_time = toc;
close(fig)