function [ ctrl, u, errc ] = mpc_qp_controller(x, u_last, p_o, dist, errc, params)
%MPC_QP_CONTROLLER implement sequential mpc qp controller
%   requires params with Np, nx, nu, dists(Np x 1), Aeq((Np + 1)*nx x (Np + 1)*nx + Np*nu
%   beq (Np+1 x 1), H square, f column, LB, UB
    % update
    
    [ A, b ] = dist_constraints(params.Np, params.nx, params.nu, p_o, dist, u_last);
    beq = [ -x; zeros(params.Np*params.nx, 1) ]; 
    
    % solve
    
    result = quadprog(params.H, params.f, A, b, params.Aeq, beq, params.LB, params.UB, [], params.options);
    if ~isempty(result)
        u = result;
        errc = 0;
    else
        errc = errc + 1;
        u = u_last;
    end
    
    ctrl = u((params.Np + 1)*params.nx + 1:end);
    ctrl = reshape(ctrl, params.nu, []);
    ctrl = ctrl(:, 1 + errc);    
end

function [ A_dist, b_dist ] = dist_constraints(Np, nx, nu, p_o, min_dists, x_mpc)
    x_mpc = x_mpc(1:(Np+1)*nx);
    x_mpc = reshape(x_mpc, nx, []);

    % past pos at k
    p_k = x_mpc([1,2], :);
    
    p_ox = p_k(:, 2:end) - p_o;
    eta = p_ox./vecnorm(p_ox);
    
    b_dist = -min_dists - eta'*p_o;
    
    A_dist = [ zeros(Np, nx), kron(speye(Np), [ 1, 1, zeros(1, nx - 2)]), zeros(Np, nu*Np) ];
    A_dist(A_dist ~= 0) = -eta;    
end