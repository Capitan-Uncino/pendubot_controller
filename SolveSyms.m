clear all
close all
clc

%pendubot's equations

syms x1 x2 x3 x4 u   %alpha1, alpha2, dalpha1/dt, dalpha2/dt 

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

%M(alpha)
m11 = p1 + p2 + 2*p3*cos(x2-x1);  m12 = p2 + p3*cos(x2-x1);
m21 = p2 + p3*cos(x2-x1);         m22 = p2;
M=[m11 m12; m21 m22];

%Vm(alpha, dalpha/dt)
Vm = p3*sin(x2-x1)*[x3-x4 -x4; x3 0];

%G(alpha)
G = [p4*g*sin(x1); p5*g*sin(x2)];

%equations
fx = [x3; x4; M\(-Vm*[x3; x4]-G)];
gx = [0; 0; M\[k; 0]];
% f34 = M\(-Vm*[x3; x4]-G);
% g34 = M\[k; 0];
% eqn = fx + gx*u == 0;
% eqns = [x3 == 0; x4 == 0; f34 + g34*u == 0];
% S = solve(eqns, [x1 x2 x3 x4 u])

% Linearization for equilibrium point x_e and input u_e
x_e = [0; 0; 0; 0];
% x_dot = fx + gx*u == 0;
% eqn = fx + gx*u == 0;
eqns = [x3 == 0; x4 == 0; M\(-Vm*[x3; x4]-G) + M\[k; 0]*u == 0];
u_e = vpa(solve(subs(eqns, x, x_e ), u), 5);
u_e = double(u_e);
% u_e = 5;

% A1 = jacobian(fx + gx*u, x)
% B1 = jacobian(fx + gx*u, u)
% A2 = vpa(subs(A, x, x_e),2)
% % A = vpa(subs(A, [x; u], [x_ref; u_e]),2)
% B2 = vpa(subs(B, x, x_ref),2)
% A = vpa(subs(A, u, S),2)
% B = vpa(subs(B, u, S),2)

A = vpa(subs(jacobian(fx + gx*u, x), [x; u], [x_e; u_e]), 3);
A = double(A);
B = vpa(subs(jacobian(fx + gx*u, u), x, x_e), 3);
B = double(B);

C = eye(n);
D = zeros(n,1);