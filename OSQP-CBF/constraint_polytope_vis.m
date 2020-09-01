minr = -2.5;
maxr = 13;

range = minr:0.01:maxr;  % plotting range
[v_ineq, omega_ineq] = meshgrid(range);  % get 2-D mesh for v and omega

cond = ones(length(range)); % Initialize
for it = 1:size(A,1)
    
    if (size(A,2) > 2) % delta exists
        if A(it,3) ~= 0
            delta_lin = [range; u(it)-A(it,1)/A(it,2)*range]';
            delta_lin(delta_lin(:,2)>maxr | delta_lin(:,2)<minr) = NaN;
        else
            cond_i = (A(it,1)*v_ineq + A(it,2)*omega_ineq < u(it)) & (A(it,1)*v_ineq + A(it,2)*omega_ineq > l(it));
        end
    else
        cond_i = (A(it,1)*v_ineq + A(it,2)*omega_ineq < u(it)) & (A(it,1)*v_ineq + A(it,2)*omega_ineq > l(it));
    end

    cond = cond & cond_i;
end

polytope = zeros(length(range)); % Initialize
polytope(cond) = NaN;
surf(v_ineq, omega_ineq, polytope, 'DisplayName', 'Feasible Region')
xlabel('v')
ylabel('\omega')
title('Feasible region')
if (size(A,2) > 2)
    hold on
    plot(delta_lin(:,1),delta_lin(:,2), 'DisplayName', '\leq \delta')
    hold off
end
view(0,90)
legend()
