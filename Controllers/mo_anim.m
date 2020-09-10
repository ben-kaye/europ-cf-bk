
anim_path = "C:\\Users\\Ben\\Documents\\EUROP\\Controllers\\_fr\\fr_%d";

GIF = 1;

drone = init_graphic(0.06, 0.03, 13);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
FPS = 30;

N_draw = ceil(N_SAMPLE/FPS/SIM_TIME);

prev_plot_a = gobjects(4);

plot_a = prev_plot_a;

tic;

N_obj = size(p_o, 2);

% axis_l = [ -3, 8, -0.5, 0, 8, 0.5];
hold on
for nc = 1:N_obj
    arg = 0:2*pi/100:2*pi;
    circ = DELTA*[ cos(arg); sin(arg) ] + p_o(:, nc);
    circ = [ circ; zeros(1, size(circ, 2)) ];
    plot3(circ(1,:), circ(2,:), circ(3,:), 'r-.', 'LineWidth', 1.5);
end
hold off
view(45,25)

axis3 = [ -1, 4, -1, 3, -1, 1 ];
fr = 1;
for Y = 1:N_SAMPLE
    if mod(Y, N_draw) == 0        
        rot_sc = 1;    
                
%         for z = 1:4
%             delete(prev_plot_a(z));
%             delete(prev_plot_b(z));
%         end
        
%         plot_a = plot_drone( [ x_t([1,2],Y); 0 ], rot_sc*theta_t(1:3, Y), drone, fig);
%         for z = 1:4
%             delete(prev_plot_a(z));
%         end
        
        plot_a = plot_drone( [ x_t(:,Y); 0 ], rot_sc*phi_t(:, Y), drone, fig);
        for z = 1:4
            delete(prev_plot_a(z));
        end
        
        prev_plot_a = plot_a;
%         prev_plot_b = plot_b;
        
        axis equal
        axis(axis3);
        
%         drawnow limitrate
        
        f_time = toc;
        if GIF 
            save_path = sprintf(anim_path, fr);
            print('-r100',save_path,'-dpng');
        else            
            if f_time < 1/FPS
                pause(1/FPS - f_time);
            end
        end
        tic;
        fr = fr + 1;
    end    
end