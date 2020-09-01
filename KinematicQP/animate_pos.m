figure(2)
axis_lims = [ -2, 4, 2, 7 ];
for it = 1:Ns
    p = x_t(:,it);

    plot(p(1), p(2),'.','MarkerSize',20)
    hold on
    plot(r_t(1,1:it),r_t(2,1:it), '-.', 'LineWidth',1)
    plot(p_o(1), p_o(2), '+', 'MarkerSize', 30, 'Color', [1,0,0],'DisplayName','Object')
    viscircles(p_o', delta, 'LineStyle', '-.');
    hold off
    
    axis(axis_lims)
    drawnow
    pause(1e-2);
       
end
close(2)