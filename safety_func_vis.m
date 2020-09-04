step = 5e-2;
eps_range = 0:step:3;
beta_range = 0:step:pi;

[ beta, epsilon ] = meshgrid(beta_range, eps_range);

h = epsilon.^2 .* (1 + sin(beta).^2) - 2;

surf(beta_range, eps_range, h)
xlabel('\beta');
ylabel('$\epsilon$', 'interpreter', 'latex');
title('h(x) visualised')