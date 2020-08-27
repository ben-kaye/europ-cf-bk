v = -1:0.01:13;  % plotting range from -5 to 5
[v, omega] = meshgrid(v);  % get 2-D mesh for x and y
conditions = (A(1,1)*v + A(1,2)*omega < u(1)) & (A(2,1)*v + A(2,2)*omega < u(2)) ...
    & (A(2,1)*v + A(2,2)*omega > l(2)) & (A(3,1)*v + A(3,2)*omega < u(3)) ...
    & (A(3,1)*v + A(3,2)*omega > l(3));
cond = zeros(length(v)); % Initialize
cond(conditions) = NaN;
surf(v, omega, cond)
xlabel('v')
ylabel('\omega')
view(0,90)