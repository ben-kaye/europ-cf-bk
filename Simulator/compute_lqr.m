function Kd = compute_lqr(Qweights, Rweights)

Q = diag(Qweights);
R = diag(Rweights);

g = 9.81;
m = 27e-3; % should use model params.. oh well
T = 1/200;

A = [   zeros(3),   eye(3),     zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0];
        zeros(3),   zeros(3),   zeros(3) ];
    
    
B =     [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        [ zeros(3,1), eye(3)] ];    
    
% State variables:
% x, y, z, xdot, ydot, zdot, gamma, beta, alpha,
% { position }, { velocity }, { euler orientation }
    
[ K, ~, ~ ] = lqr(A, B, Q, R);
[ Kd, ~, ~] = lqrd(A, B, Q, R, 0, T);

K(abs(K)<1e-10) = 0;
Kd(abs(Kd)<1e-10) = 0;

end