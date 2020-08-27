range = -2.5:0.01:13;  % plotting range
[v_ineq, omega_ineq] = meshgrid(range);  % get 2-D mesh for v and omega

cond = ones(length(range)); % Initialize
for it = 1:size(A,1)
    cond_i = (A(it,1)*v_ineq + A(it,2)*omega_ineq < u(it)) & (A(it,1)*v_ineq + A(it,2)*omega_ineq > l(it));
    cond = cond & cond_i;
end

polytope = zeros(length(range)); % Initialize
polytope(cond) = NaN;
surf(v_ineq, omega_ineq, polytope)
xlabel('v')
ylabel('\omega')
view(0,90)