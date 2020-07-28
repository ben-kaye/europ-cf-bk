[ A, l, etas ] = update_constrs(A, l, N, nx, nu, o_states, dist, xN, 0);
prob.update('Ax', nonzeros(A));
res = prob.solve();
xN = res.x;
vis_predicted_path;
