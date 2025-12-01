function xpo=pendubot_eqt(u,xi)

%||=================================================||
%||                                                 ||
%||   LINMA 2671                                    ||
%||                                                 ||
%||   Pr J. Hendrickx, Eng. A. Taylor & F. Wielant  ||
%||                                                 ||
%||=================================================||

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
Vmda = p3*sin(x(2)-x(1))*[-x(4)*x(3) + x(3)^2 - x(4)^2 ; x(3)^2];

%G(alpha)
G = [p4*g*sin(x(1)); p5*g*sin(x(2))];

% new state vector
xp(1) = x(3);
xp(2) = x(4);
xp(3:4) = inv(M)*([k*u ; 0] - Vmda - G);

xp=xp';

%reverse order
xpo=[xp(1); xp(3); xp(2); xp(4)];
