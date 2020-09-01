t = 0:step_size*N/Ns:sim_time;
t = t(1:end-1);
figure(3)
subplot(1,2,1);
plot(t, h_t,'DisplayName','h','LineWidth',1.5, 'Color', 1/255*[71, 191, 189]);
yline(0, 'LineStyle','-.','Color',[1,0,0])
xlabel('t (s)')
ylim([-1, max(h_t)+1])
title('Safety Function')
subplot(1,2,2);

plot(t, u_t(1,:),'DisplayName','v', 'LineWidth',1.5, 'Color',1/255*[102, 145, 58]);
hold on
plot(t, u_t(2,:),'DisplayName','\omega','LineWidth',1.5, 'Color',1/255*[142, 17, 184])
legend()
xlabel('t (s)')
ylabel('ms^{-1} / rads^{-1}')
title('Controls')
hold off