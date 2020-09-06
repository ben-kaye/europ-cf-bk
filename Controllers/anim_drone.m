
drone = init_graphic();
fig = figure();

FPS = 40;

N_draw = ceil(N_SAMPLE/FPS/SIM_TIME);
prev_plot_a = gobjects(4);
prev_plot_b = gobjects(4);

tic;

% axis_l = [ -3, 8, -0.5, 0, 8, 0.5];
arg = 0:2*pi/100:2*pi;
circ = DELTA*[ cos(arg); sin(arg) ] + p_o;
circ = [ circ; zeros(1, size(circ, 2)) ];
plot3(circ(1,:), circ(2,:), circ(3,:), 'r-.', 'LineWidth', 1.5);



axis3 = [ -3, 8, 0, 10, -0.5, 0.5 ];

for Y = 1:N_SAMPLE
    if mod(Y, N_draw) == 0
        
        
        rot_sc = 1;
        
        
        
        plot_a = plot_drone( [ x_t([1,2],Y); 0 ], rot_sc*theta_t(1:3, Y), drone, fig);
%         plot_b = plot_drone( [ x_t([3,4],Y); 0 ], rot_sc*theta_t(4:6, Y), drone, fig);
        
        for z = 1:4
            delete(prev_plot_a(z));
%             delete(prev_plot_b(z));
        end
        
        prev_plot_a = plot_a;
%         prev_plot_b = plot_b;
        
        axis equal
        axis(axis3);
        
%         drawnow;
        
        f_time = toc;
        
        if f_time < 1/FPS
            pause(1/FPS - f_time);
        end
        
        tic;
    end    
end