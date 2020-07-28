x_N = reshape(xN(1:(N+1)*nx), nx, []);

hold on
plot(x_N(1,:), x_N(2,:), '.', 'MarkerSize', 20, 'Color', [1, 0.6, 0])
axis equal