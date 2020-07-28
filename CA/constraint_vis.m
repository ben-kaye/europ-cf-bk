function constraint_vis(etas, obj_xy, min_dist, k)

etak = etas(:,k);

grad = cross([0; 0; 1], [etak; 0]);
grad = grad(1:2);
m = grad(2)/grad(1);

t = obj_xy + min_dist * etak;

x_p = linspace(obj_xy(1)-2*min_dist, obj_xy(1)+2*min_dist, 8);

y_p = t(2) + m * (x_p - t(1));

hold on
plot(x_p, y_p);

end