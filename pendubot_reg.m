function u=pendubot_reg(alpha1r,x)
% %function u=pendubot_reg(alpha21,x)
% %
% %INPUTS:
% % alpha1r : reference value for the first link angle [rad]
% % x       : measured state vector x =   [
% %										x1 = alpha_1 - pi
% %										x2 = alpha_2 - pi
% %										x3 = dalpha_1/dt
% %										x4 = dalpha_2/dt
% %                                       ]
% %
% %OUTPUT:
% % u       : control signal (motor torque) to be applied [Nm]
% 
% 
% 
% %u=3*((alpha1r-pi)-x(1))-0.3*x(3);
% x_e = [alpha1r; pi; 0; 0];
% % x_e = [alpha1r-pi; 0; 0; 0];
% x_r = [pi; pi; 0; 0];
% u_e = -1.7205;
% K = [-0.2165, -16.2534, -1.2892, -2.1864];
% 
% u = u_e - K*(x'+x_r-x_e);
% % u = u_e - K*(x-x_e);
% % u = 0;

%========================================================
% SWING UP
% xi = x'+ [pi; pi; 0; 0];
xi = [x(1)+pi; x(2)+pi; x(3); x(4)];

p1=0.0148;%kg m^2
p2=0.0051;
p3=0.0046;
p4=0.1003;%kg m
p5=0.0303;
g=9.81;

k = 3.9621/8;%V (Nm)^â?»1

sym x1 x2 x3 x4 v1 k1 k2

% xs = [x1; x2; x3; x4];

% New state variables
z1 = x1-alpha1r;
eta1 = x2;
z2 = x3;
eta2 = x4;

w = [z1; z2; eta1; eta2];

% Equilibrium point
x_e = [alpha1r; pi; 0; 0];

%M(alpha)
m11 = p1 + p2 + 2*p3*cos(eta1-z1);  m12 = p2 + p3*cos(eta1-z1);
m21 = p2 + p3*cos(eta1-z1);         m22 = p2;

%Vm(alpha, dalpha/dt)*dalpha/dt
Vmda = p3*sin(eta1-z1)*[-eta2*z2 + z2^2 - eta2^2 ; z2^2];

%G(alpha)
G = [p4*g*sin(z1); p5*g*sin(eta1)];

% Solve for equilibrium input
eqns = [z2 == 0; z1 == 0; eta2 == 0; -1/m22 * (Vmda(2)+G(2)) - (m12/m22 * v1) == 0];
v1_e = double(vpa(solve(subs(eqns, w, x_e ), v1), 5));

[kp, kd] = solve(k1 > 0, k2 > 0, v1_e == k1*(alpha1r-x1)-k2*x3, x, x_e);

V1 = kp*(alpha1r-x1)-kd*x3;

% Partial Feedback Linearization
%M(alpha)
m11 = p1 + p2 + 2*p3*cos(xi(2)-xi(1));  m12 = p2 + p3*cos(xi(2)-xi(1));
m21 = p2 + p3*cos(xi(2)-xi(1));         m22 = p2;

%Vm(alpha, dalpha/dt)*dalpha/dt
Vmda = p3*sin(xi(2)-xi(1))*[-xi(4)*xi(3) + xi(3)^2 - xi(4)^2 ; xi(3)^2];

%G(alpha)
G = [p4*g*sin(xi(1)); p5*g*sin(xi(2))];

m11_hat = m11 - (m12*m21/m22);
Vmda1_hat = Vmda(1) - (m12*Vmda(2)/m22);
G1_hat = G(1) - (m12*G(2)/m22);

tau = m11_hat*V1 + Vmda1_hat + G1_hat;
u = tau/k;

