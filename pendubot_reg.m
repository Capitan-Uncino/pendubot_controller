function u=pendubot_reg(alpha1r,x)
%function u=pendubot_reg(alpha21,x)
%
%INPUTS:
% alpha1r : reference value for the first link angle [rad]
% x       : measured state vector x =   [
%										x1 = alpha_1 - pi
%										x2 = alpha_2 - pi
%										x3 = dalpha_1/dt
%										x4 = dalpha_2/dt
%                                       ]
%
%OUTPUT:
% u       : control signal (motor torque) to be applied [Nm]



%u=3*((alpha1r-pi)-x(1))-0.3*x(3);
% x_e = [alpha1r; pi; 0; 0];
x_e = [alpha1r-pi; 0; 0; 0];
% x_r = [pi; pi; 0; 0];
u_e = -1.7205;
K = [0.7076, -5.8825, -0.2986, -0.7025];

% u = u_e - K*(x+x_r-x_e);
u = u_e - K*(x-x_e);
% u = 0;