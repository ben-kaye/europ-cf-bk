% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *                 Program by Ben Kaye (c) 2020                  *
% *            Using model provided by Aren Karapetyan            *
% *                         EUROP Project                         *
% *                                                               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

clear,clc

%System State-Space

%d_system.A([1:2, 4:5, 7:8], [1:2, 4:5, 7:8])
A = [   1,  0,  0.0962912484723868,     0,                      0,                      0.0396218015830554;
        0,  1,  0,                      0.0962912484723868,     -0.0396218015830554,    0;
        0,  0,  0.894293322749997,      0,                      0,                      0.702694931894877;
        0,  0,  0,                      0.894293322749997,      -0.702694931894877,     0;
        0,  0,  0,                      0.193245482770890,      0.452393730892869,      0;
        0,  0,  -0.193245482770890,     0,                      0,                      0.452393730892869   ];

%d_system.B([1:2, 4:5, 7:8], [4, 5])
B = [   0.00370875152761323,    0;
        0,                      0.00370875152761323;
        0.105706677250003,      0;
        0,                      0.105706677250003;
        0,                      -0.193245482770890;
        0.193245482770890,      0                       ];

%Calculate MPC algo

%initial conditions
ref = [ 3; 0; 0; 0; 0; 0 ];
x_0 = [ -4; 3; 0; 0; 0; 0 ];

figure(1);
plot(ref(1), ref(2), 'b+')
axis([-5 5 -5 5])
% axis equal


%prediction horizon
N = 10;

%sim time
T = 100;

Q = eye(6);
R = 1e-2 * eye(2);

x = x_0;




for i = 1:T
    %returns Algebraic-Ricatti iterator for k + 1
    Pplus = RicattiIter(A,B,Q,R,Q,N);
    %compute gain matrix from previous iterator
    K = -(B'*Pplus*B + R)\(B'*Pplus*A);
    u = K*(x-ref);
    
    x = A*x + B*u;
    
    hold on    
    plot(x(1), x(2), 'bo')
    axis([-5 5 -5 5])
%     axis equal
    
end

%check this is in fact getting k-1 looks to be right
function Pnmin = RicattiIter(A,B,Q,R,P,n)
    n = n - 1;
    Pnmin = P;
    if n > 0 % stop at N-1 iterations
        Pnmin = Q + A'*P*A - A'*P*B*((B'*P*B + R)\(B'*P*A));
        Pnmin = RicattiIter(A,B,Q,R,Pnmin,n);
%         n
    end
    
end

