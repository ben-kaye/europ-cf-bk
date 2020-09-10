x_pred = qp_res(1:(Np+1)*nx);
x_pred = reshape(x_pred, nx, []);


viscircles(p_o', delta*ones(size(p_o, 2), 1), 'LineStyle', '-.');
hold on
plot(r_mpc(1), r_mpc(2), 'o');
plot(x_t(1,:), x_t(2,:))
axis equal
plot(x_pred(1, 2:end), x_pred(2, 2:end));