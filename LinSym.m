function xpo = LinSym(u,xi)
% function [A, B, C, D] = LinSym(xi,u)

%pendubot's equations

x=[xi(1); xi(3); xi(2); xi(4)]; % state vector alpha1, alpha2, dalpha1/dt, dalpha2/dt

p1=0.0148;%kg m^2
p2=0.0051;
p3=0.0046;
p4=0.1003;%kg m
p5=0.0303;
g=9.81;

k = 3.9621/8;%V (Nm)^â?»1

%M(alpha)
m11 = p1 + p2 + 2*p3*cos(x(2)-x(1));  m12 = p2 + p3*cos(x(2)-x(1));
m21 = p2 + p3*cos(x(2)-x(1));         m22 = p2;
M=[m11 m12; m21 m22];

%Vm(alpha, dalpha/dt)*dalpha/dt
Vm = p3*sin(x(2)-x(1))*[x(3)-x(4),  -x(4);  x(3), 0];

%G(alpha)
G = [p4*g*sin(x(1)); p5*g*sin(x(2))];

%equations 
fx = [x(3); x(4); M\(-Vm*[x(3); x(4)] - G)];
gx = [0; 0; M\[k; 0]];

x_dot = fx + gx*u == 0;
u_e = vpa(solve(subs(x_dot, x, x_e), u),5)

A = vpa(subs(jacobian(fx + gx*u, x), u, u_e),2)
B = vpa(subs(jacobian(fx + gx*u, u), u, u_e),2)

% new state vector
xpo = A*(x-x_e) + B*(u-u_e)