

axis_lims = [ -2, 4, 2, 7 ];

for it = 1:Ns
    
    p = x_t(:,it);
    
    figure(1)
    plot(p(1), p(2),'.','MarkerSize',20)
    axis(axis_lims)
    drawnow
    pause(1e-2);
    
    
end