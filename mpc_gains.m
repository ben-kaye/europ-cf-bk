function [ Ae, Be, Ce, PhiT_Phi, PhiT_F, PhiT_R ] = mpc_gains( A, B, C, Nc, Np)
    
    [m1, n1] = size(C);
    [~, n_in] = size(B);
    
    Ae = eye(n1+m1);
    Ae(1:n1, 1:n1) = A;
    Ae(n1+1:n1+m1, 1:n1) = C*A;
    
    Be = zeros(n1+m1, n_in);
    Be(1:n1, :) = B;
    Be(n1+1:n1+m1, :) = C*B;
    
    Ce = zeros(m1,n1+m1);
    Ce(:, n1+1:n1+m1) = eye(m1,m1);
    
%     n = n1 + m1;
    
    h(1,:) = Ce;
    F(1,:) = Ce*Ae;
    
    for p = 2:Np 
        h(p, :) = h(p-1, :) * Ae;
        F(p, :) = F(p-1, :) * Ae;
    end
    
    
    %Construct Hessian
    v = h*Be;
    Phi = zeros(Np, Nc);
    Phi(:,1) = v;
    for p = 2:Nc
        Phi(p:end, p) = v(1:Np-p+1);       
    end
    
    RsBar = ones(Np,1);
    PhiT_Phi = Phi' * Phi;
    PhiT_F = Phi' * F;
    PhiT_R = Phi' * RsBar;    
end
