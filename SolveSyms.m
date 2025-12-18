clear all
close all
clc

%pendubot's equations

syms x1 x2 x3 x4 v1 k1 k2  %alpha1, alpha2, dalpha1/dt, dalpha2/dt 

x = [x1; x2; x3; x4];
n = length(x);
x_ini = [pi/24; 0; 0; 0];
x_ref = [4/3*pi; pi; 0; 0];

p1=0.0148;%kg m^2
p2=0.0051;
p3=0.0046;
p4=0.1003;%kg m
p5=0.0303;
g=9.81;

k = 3.9621/8;%V (Nm)^â?»1

% %M(alpha)
% m11 = p1 + p2 + 2*p3*cos(x2-x1);  m12 = p2 + p3*cos(x2-x1);
% m21 = p2 + p3*cos(x2-x1);         m22 = p2;
% M=[m11 m12; m21 m22];
% 
% %Vm(alpha, dalpha/dt)
% Vm = p3*sin(x2-x1)*[x3-x4 -x4; x3 0];
% 
% %G(alpha)
% G = [p4*g*sin(x1); p5*g*sin(x2)];
% 
% %equations
% fx = [x3; x4; M\(-Vm*[x3; x4]-G)];
% gx = [0; 0; M\[k; 0]];
% 
% % Linearization for equilibrium point x_e and input u_e
% x_e = [4/3*pi; pi; 0; 0]; 
% 
% eqns = [x3 == 0; x4 == 0; M\(-Vm*[x3; x4]-G) + M\[k; 0]*u == 0];
% u_e = double(vpa(solve(subs(eqns, x, x_e ), u), 5));
% 
% A = double(vpa(subs(jacobian(fx + gx*u, x), [x; u], [x_e; u_e]), 3));
% B = double(vpa(subs(jacobian(fx + gx*u, u), x, x_e), 3));
% 
% C = eye(n);
% D = zeros(n,1);

xs = [x1; x2; x3; x4];

% New state variables
z1 = x1-pi;
eta1 = x2;
z2 = x3;
eta2 = x4;

w = [z1; z2; eta1; eta2];

% Equilibrium point
x_e = [pi; pi; 0; 0];

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

% [kp, kd] = solve(k1 >= 0, k2 >= 0, v1_e == k1*(pi-x1)-k2*x3, x_e);

V1 = kp*(pi-x1)-kd*x3;