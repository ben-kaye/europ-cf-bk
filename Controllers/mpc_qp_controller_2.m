function [ ctrl, u, errc ] = mpc_qp_controller_2(x, u_last, p_o, dist, errc, params)
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
    if errc < params.Np
        ctrl = ctrl(:, 1 + errc); 
    else
        ctrl = zeros(2, 1);
    end
end

function [ A_dist, b_dist ] = dist_constraints(Np, nx, nu, p_o, min_dists, x_mpc)
    x_mpc = x_mpc(1:(Np+1)*nx);
    x_mpc = reshape(x_mpc, nx, []);
    
    Nobj = size(p_o, 2);
    
    % past pos at k
    p_k = x_mpc([1,2], :);

    mask_A = [ zeros(Np, nx), kron(speye(Np), [ 1, 1, zeros(1, nx - 2)]), zeros(Np, nu*Np) ];
    
    if Nobj > 0
        
        A_dist = zeros(Np*Nobj, (Np + 1)*nx + Np*nu);
        b_dist = zeros(Np*Nobj, 1);
        
        A_temp = zeros(size(mask_A));
        
        for z = 1:Nobj

            p_z = p_o(:, z);

            p_ox = p_k(:, 2:end) - p_z;
            eta = p_ox./vecnorm(p_ox);

            b_temp = -min_dists - eta'*p_z;
            A_temp(mask_A ~= 0) = -eta;    
           
            A_dist((z-1)*Np+1:z*Np, :) = A_temp;
            b_dist((z-1)*Np+1:z*Np) = b_temp;
            
        end
    end
end