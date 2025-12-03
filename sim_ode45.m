
%% ORIGINAL CODE
clear all
close all
clc

% Pendubot equations
syms x1 x2 x3 x4 u   % alpha1, alpha2, dalpha1/dt, dalpha2/dt 

x = [x1; x2; x3; x4];
n = length(x);
x_ini = [pi/24; 0; 0; 0];

p1=0.0148;%kg m^2
p2=0.0051;
p3=0.0046;
p4=0.1003;%kg m
p5=0.0303;
g=9.81;

k = 3.9621/8;%V (Nm)^−1

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

% Linearization point
x_e = [0; 0; 0; 0];

% Solve for equilibrium input
eqns = [x3 == 0; x4 == 0; M\(-Vm*[x3; x4]-G) + M\[k; 0]*u == 0];
u_e = double(vpa(solve(subs(eqns, x, x_e ), u), 5));

% Linearized system
A = double(vpa(subs(jacobian(fx + gx*u, x), [x; u], [x_e; u_e]), 3));
B = double(vpa(subs(jacobian(fx + gx*u, u), x, x_e), 3));

C = eye(n);
D = zeros(n,1);

eig_A = eig(A)

omega_n = abs(eig_A);

omega_fast = max(omega_n)


f_fast = omega_fast / (2*pi);  
f_s    = 10 * f_fast            % sampling frequency (Hz)
Ts     = 1 / f_s                % sampling period (s)

sys_c = ss(A, B, C, D);

sys_d = c2d(sys_c, Ts, 'zoh');

[Ad, Bd, Cd, Dd] = ssdata(sys_d)

%% ---------------------------------------------------------------
%% SIMULATION OF LINEAR SYSTEM x_dot = A x + B u
%% ---------------------------------------------------------------

% Simulation parameters
Tend = 0.3;
dt = 0.001;
tspan = 0:dt:Tend;

% Define input: constant equilibrium input
u_fun = @(t) -u_e;

% Linear system dynamics for ode45
lin_dyn = @(t, x) A*(x - x_e) + B*u_fun(t);

% Initial condition (offset from equilibrium)
x0 = x_ini-x_e;   % use your original initial state

% Simulate
[t, X] = ode45(lin_dyn, tspan, x0);

X = X + repmat(x_e', size(X,1), 1);


%% ---------------------------------------------------------------
%% SIMULATION OF NONLINEAR SYSTEM
%% ---------------------------------------------------------------






% Simulation parameters
Tend = 0.3;
dt = 0.001;
tspan = 0:dt:Tend;

% Define input: constant equilibrium input
u_fun = @(t) 0;



fx_fun = matlabFunction(fx, 'Vars', {x1, x2, x3, x4});
gx_fun = matlabFunction(gx, 'Vars', {x1, x2, x3, x4});


% Linear system dynamics for ode45
dyn = @(t, x, u) fx_fun(x(1), x(2), x(3), x(4)) + gx_fun(x(1), x(2), x(3), x(4)) * u;

% Initial condition (offset from equilibrium)
x0 = x_ini;   % use your original initial state

% Simulate
[t_nonlinear, X_nonlinear] = ode45(@(t,x) dyn(t,x,u_fun(t)), tspan, x0);




%% ---------------------------------------------------------------
%% PLOTTING
%% ---------------------------------------------------------------

figure;
subplot(2,2,1)
plot(t, X(:,1), 'LineWidth', 1.6); hold on
plot(t, X(:,2), 'LineWidth', 1.6);
legend('\alpha_1','\alpha_2')
ylabel('Angle [rad]')
title('Pendubot Angles (Linearized Model)')

subplot(2,2,2)
plot(t, X(:,3), 'LineWidth', 1.6); hold on
plot(t, X(:,4), 'LineWidth', 1.6);
legend('d\alpha_1','d\alpha_2')
ylabel('Angular Velocity [rad/s]')
xlabel('Time [s]')
title('Pendubot Angular Velocities (Linearized Model)')

subplot(2,2,3)
plot(t, X_nonlinear(:,1), 'LineWidth', 1.6); hold on
plot(t, X_nonlinear(:,2), 'LineWidth', 1.6);
legend('\alpha_1','\alpha_2')
ylabel('Angle [rad]')
title('Pendubot Angles')

subplot(2,2,4)
plot(t, X_nonlinear(:,3), 'LineWidth', 1.6); hold on
plot(t, X_nonlinear(:,4), 'LineWidth', 1.6);
legend('d\alpha_1','d\alpha_2')
ylabel('Angular Velocity [rad/s]')
xlabel('Time [s]')
title('Pendubot Angular Velocities')
