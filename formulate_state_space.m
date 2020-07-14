% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
% *                                                               *
% *                 Program by Ben Kaye (c) 2020                  *
% *                         EUROP Project                         *
% *                                                               *
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

g = 9.81;
m = 27e-3; % should use model params.. oh well
Ts = 1/10;
cmd_2_newtons_conversion_quadratic_coefficient  =  1.3385e-10;
cmd_2_newtons_conversion_linear_coefficient     =  6.4870e-6;
nrotor_vehicle_thrust_max = 0.1597;


% LTI continuous-time system
A = [   zeros(3),   eye(3),     zeros(3);
        zeros(3),   zeros(3),   [0, g, 0; -g, 0, 0; 0, 0, 0];
        zeros(3),   zeros(3),   zeros(3) ];
    
B = [   zeros(3,4);
        [ [ 0; 0; 1/m ], zeros(3) ];
        [ zeros(3,1), eye(3)] ];
     
% C = [ diag(ones(1,3)), zeros(3,6) ];
C = eye(9);

D = 0;

plant = ss(A,B,C,D);
plant.InputName = 'u';
plant.OutputName = 'y';


[ Ac, Bc, Cc, Dc ] = linmod('crazyflie_controller_ak');

controller = ss(Ac,Bc,Cc,Dc);
controller.InputName = 'e';
controller.OutputName = 'u';

sum = sumblk('e = r-y', 9);

system = connect(plant, controller, sum, 'r', 'y');

d_system = c2d(system, Ts);

TF = minreal(zpk(d_system));


Nc = 4;
Np = 10;

% [ Ae, Be, Ce, PhiT_Phi, PhiT_F, PhiT_R ] = mpc_gains( Ad, Bd, Cd, Nc, Np);





% x_0 = [0.1; 0.2];
% deltaU = (PhiT_Phi + r_w * eye(Nc))\(PhiT_R - PhiT_F*x_0);


