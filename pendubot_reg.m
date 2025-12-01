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



u=3*((alpha1r-pi)-x(1))-0.3*x(3);