
%% ORIGINAL CODE
clear all
close all
clc

% Pendubot equations
syms x1 x2 x3 x4 u   % alpha1, alpha2, dalpha1/dt, dalpha2/dt 

x = [x1; x2; x3; x4];
n = length(x);
x_ini = [pi+pi/6; pi; 0; 0];
Tend = 10;

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
x_e = [pi; pi; 0; 0];

% Solve for equilibrium input
eqns = [x3 == 0; x4 == 0; M\(-Vm*[x3; x4]-G) + M\[k; 0]*u == 0];
u_e = double(vpa(solve(subs(eqns, x, x_e ), u), 5));

% Linearized system
A = double(vpa(subs(jacobian(fx + gx*u, x), [x; u], [x_e; u_e]), 3));
B = double(vpa(subs(jacobian(fx + gx*u, u), x, x_e), 3));

C = eye(n);
D = zeros(n,1);



%% ---------------------------------------------------------------
%% SIMULATION OF LINEAR SYSTEM x_dot = A x + B u
%% ---------------------------------------------------------------

% Simulation parameters
dt = 0.001;
tspan = 0:dt:Tend;

% Define input: constant equilibrium input
u_fun = @(t,x) -u_e;

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
dt = 0.001;
tspan = 0:dt:Tend;

% Define input: constant equilibrium input
u_fun = @(t,x) 0;



fx_fun = matlabFunction(fx, 'Vars', {x1, x2, x3, x4});
gx_fun = matlabFunction(gx, 'Vars', {x1, x2, x3, x4});


% Linear system dynamics for ode45
dyn = @(t, x, u) fx_fun(x(1), x(2), x(3), x(4)) + gx_fun(x(1), x(2), x(3), x(4)) * u;

% Initial condition (offset from equilibrium)
x0 = x_ini;   % use your original initial state

% Simulate
[t_nonlinear, X_nonlinear] = ode45(@(t,x) dyn(t,x,u_fun(t)), tspan, x0);


%% ---------------------------------------------------------------
%% SIMULATION OF DISCRETIZED SYSTEM
%% ---------------------------------------------------------------


eig_A = eig(A);

omega_n = abs(eig_A);

omega_fast = max(omega_n);


f_fast = omega_fast / (2*pi)  
f_s    = 10 * f_fast
Ts     =  floor(1 / f_s * 100) / 100

sys_c = ss(A, B, C, D);

sys_d = c2d(sys_c, Ts, 'zoh');

[Ad, Bd, Cd, Dd] = ssdata(sys_d);


N = floor(Tend/Ts) + 1;  
t_discrete = (0:N-1)*Ts;


X_discrete = zeros(N, size(Ad,1));
U_discrete = zeros(N, 1);

Q = diag([1,1,1,1])
R = 1

K = dlqr(Ad,Bd,Q,R)


n = rows(Ad);    % = 4
m = columns(Bd);

% ---------- Weights ----------
s = tf('s');

% Continuous weights (choose as you like)
Ws_c = (s/2 + 1)/(s + 0.01);    % sensitivity weight
Wu_c = 0.1;                     % control weight
Wt_c = (s + 200)/(s/200 + 1);   % complementary sensitivity weight

% Discretize
Ws = c2d(Ws_c, Ts);
Wu = c2d(Wu_c, Ts);
Wt = c2d(Wt_c, Ts);

% ---------- Generalized Plant (Mixed Sensitivity) ----------

% e = r - y  (y = x)
% Build e-system
E = ss(-eye(n), eye(n), eye(n), zeros(n,n), Ts);
% inputs: [u ; r]  → outputs: e

% You can treat r as exogenous (size 4)

% Block partitions (standard form)
P11 = [ Ws*E ;
        Wu*ss([],[],[],eye(m),Ts) ;
        Wt*Gd ];

P12 = [ -Ws*Gd ;
         Wu ;
         Wt*Gd ];

P21 = E;
P22 = -Gd;

% Full generalized plant
P = [ P11  P12 ;
      P21  P22 ];

% ---------- H-infinity synthesis ----------
nmeas = n;      % number of measured outputs = 4
ncont = m;      % number of control inputs

[K, CL, gamma] = hinfsyn(P, nmeas, ncont, Ts);

disp("H-infinity gamma = "), disp(gamma);


x = x_ini- x_e;

  u_fun = @(x,t) -K*x -u_e;  

for k = 1:N
    t_curr = (k-1)*Ts;       
    X_discrete(k,:) = x;         
    
    u = u_fun(x, t_curr);    
    U_hist(k,:) = u;         
    
    x = Ad*x + Bd*u;          
end

X_discrete = X_discrete + repmat(x_e', size(X_discrete,1), 1);

%% ---------------------------------------------------------------
%% PLOTTING
%% ---------------------------------------------------------------


visualize_systems = [3]

  scale_angles_low = 2; 
  scale_angles_high = 4;
  scale_velocities_low = -2;
  scale_velocities_high = 2;


if any(ismember(visualize_systems, 1))



  figure(1);
  subplot(2,1,1)
  plot(t_nonlinear, X_nonlinear(:,1), 'LineWidth', 1.6); hold on
  plot(t_nonlinear, X_nonlinear(:,2), 'LineWidth', 1.6);
  legend('\alpha_0','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles')
  ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(t_nonlinear, X_nonlinear(:,3), 'LineWidth', 1.6); hold on
  plot(t_nonlinear, X_nonlinear(:,4), 'LineWidth', 1.6);
  legend('d\alpha_0','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities')
  ylim([scale_velocities_low scale_velocities_high])   % <-- SET SCALE

end

if any(ismember(visualize_systems, 2))

  figure(2);
  subplot(2,1,1)
  plot(t, X(:,1), 'LineWidth', 1.6); hold on
  plot(t, X(:,2), 'LineWidth', 1.6);
  legend('\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles (Linearized Model)')
  ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(t, X(:,3), 'LineWidth', 1.6); hold on
  plot(t, X(:,4), 'LineWidth', 1.6);
  legend('d\alpha_1','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities (Linearized Model)')
  ylim([scale_velocities_low scale_velocities_high])    % <-- SET SCALE

end

if any(ismember(visualize_systems, 3))


  figure(3);
  subplot(2,1,1)
  plot(t_discrete, ones(length(X_discrete)).*pi, 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,1), 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,2), 'LineWidth', 1.6);
  legend('\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles (Discretized Model)')
  ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(t_discrete, X_discrete(:,3), 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,4), 'LineWidth', 1.6);
  legend('d\alpha_1','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities (Discretized Model)')
  ylim([scale_velocities_low scale_velocities_high])   % <-- SET SCALE

end
