close all
clc

%Plot after running pendubot

% Pendubot equations
syms x1 x2 x3 x4 u   % alpha1, alpha2, dalpha1/dt, dalpha2/dt 

x = [x1; x2; x3; x4];
n = length(x);
x_ini = [pi; pi; 0; 0];
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
x_e = [4*pi/3; pi; 0; 0];

% Solve for equilibrium input
eqns = [x3 == 0; x4 == 0; M\(-Vm*[x3; x4]-G) + M\[k; 0]*u == 0];
u_e = double(vpa(solve(subs(eqns, x, x_e ), u), 5));

% Linearized system
A = double(vpa(subs(jacobian(fx + gx*u, x), [x; u], [x_e; u_e]), 3));
B = double(vpa(subs(jacobian(fx + gx*u, u), x, x_e), 3));

C = eye(n);
D = zeros(n,1);

% READ TABLE
% opts = detectImportOptions("Test_1.txt");
% opts.DecimalSeparator = ",";
% opts.Delimiter = " ";
% opts = setvartype(opts, "double");
% T = readtable("Test_1.txt", opts);
% fid = fopen("Test_1.txt","r");
% file = textscan(fid, "%f %f %f %f %f %f %f %f","HeaderLines",1,'Delimiter',' ','MultipleDelimsAsOne',1)
% fclose(fid);
% === READ TABLE ===

% 1) leggi tutto come testo
fid = fopen("test_lqr.txt.txt","r","n","UTF-8");
txt = fread(fid,"*char")';
fclose(fid);

% 2) sostituisci virgola -> punto
txt = strrep(txt, ",", ".");

% 3) file temporaneo corretto
fid = fopen("test_lqr.txt.txt","w");
fwrite(fid, txt);
fclose(fid);

% 4) leggi i dati numerici (salta intestazioni)
fid = fopen("test_lqr.txt.txt","r");

data = textscan(fid, ...
    "%f %f %f %f %f %f %f %f", ...
    "HeaderLines", 3, ...
    "Delimiter", {'\t',' '}, ...
    "MultipleDelimsAsOne", true, ...
    "CollectOutput", true);

fclose(fid);

M = data{1};

% % 1) leggi tutto il file come testo
% fid = fopen("Test_1.txt","r","n","UTF-8");
% txt = fread(fid,"*char")';
% fclose(fid);
% 
% % 2) sostituisci virgola con punto
% txt = strrep(txt, ",", ".");
% 
% % 3) dividi in righe e salta intestazioni
% lines = strsplit(txt, {'\r','\n'});
% lines = lines(4:end);  % le prime 3 righe sono intestazioni
% 
% % 4) inizializza matrice dati
% M = [];
% 
% for k = 1:length(lines)
%     line_clean = strtrim(lines{k});
%     if isempty(line_clean)
%         continue
%     end
%     % sostituisce tab e spazi multipli con uno spazio singolo
%     line_clean = regexprep(line_clean,'[\t ]+',' ');
%     % split riga in celle di testo
%     elems = strsplit(line_clean,' ');
%     % converte in numeri
%     nums = str2double(elems);
%     % aggiunge alla matrice se tutti sono numerici
%     if all(~isnan(nums))
%         M = [M; nums];
%     end
% end
% 
% % 5) converte in table con nomi colonne
% nomi = {'t','U','Us','Alpha1_r','Alpha1','dAlpha1','Alpha2','Curr_Monitor'};
% T = array2table(M,"VariableNames",nomi);

% % 1) leggi tutto il file come testo
% fid = fopen("test_lqr.txt.txt","r","n","UTF-8");
% if fid == -1
%     error('Impossibile aprire il file');
% end
% txt = fread(fid,"*char")';
% fclose(fid);
% 
% % 2) sostituisci virgola con punto
% txt = strrep(txt, ",", ".");
% 
% % 3) dividi il testo in righe
% lines = strsplit(txt, {'\r','\n'});
% 
% % 4) prendi la riga di intestazione (seconda riga) e prepara i nomi colonne
% header_line = strtrim(lines{2});
% header_line = regexprep(header_line,'[\t ]+',' '); % TAB/spazi multipli → 1 spazio
% col_names = strsplit(header_line,' ');
% 
% % 5) leggi le righe dei dati (da riga 3 in poi)
% M = [];
% for k = 3:length(lines)
%     line_clean = strtrim(lines{k});
%     if isempty(line_clean)
%         continue
%     end
%     % TAB/spazi multipli → 1 spazio
%     line_clean = regexprep(line_clean,'[\t ]+',' ');
%     elems = strsplit(line_clean,' ');
%     nums = str2double(elems);
%     if all(~isnan(nums))
%         M = [M; nums];
%     end
% end
% 
% % % 6) crea la table con nomi colonne corretti
% % T = array2table(M, "VariableNames", col_names);
% % 
% % % 7) mostra prime righe per controllo
% % T(1:10,:)


%% ---------------------------------------------------------------
%% SIMULATION OF LINEAR SYSTEM x_dot = A x + B u
%% ---------------------------------------------------------------

% Simulation parameters
dt = 0.001;
tspan = 0:dt:Tend;

% Define input: constant equilibrium input
u_fun = @(t,x) - u_e;

% Linear system dynamics for ode45
lin_dyn = @(t, x) A*(x - x_e) + B*u_fun(t);

% Initial condition (offset from equilibrium)
x0 = x_ini-x_e;   % use your original initial state

% Simulate
[t, X] = ode45(lin_dyn, tspan, x0);

X = X + repmat(x_e', size(X,1), 1);

% %% ---------------------------------------------------------------
% %% SIMULATION OF NONLINEAR SYSTEM
% %% ---------------------------------------------------------------
% 
% % Simulation parameters
% dt = 0.001;
% tspan = 0:dt:Tend;
% 
% % Define input: constant equilibrium input
% u_fun = @(t,x) 0;
% 
% 
% 
% fx_fun = matlabFunction(fx, 'Vars', {x1, x2, x3, x4});
% gx_fun = matlabFunction(gx, 'Vars', {x1, x2, x3, x4});
% 
% 
% % Linear system dynamics for ode45
% dyn = @(t, x, u) fx_fun(x(1), x(2), x(3), x(4)) + gx_fun(x(1), x(2), x(3), x(4)) * u;
% 
% % Initial condition (offset from equilibrium)
% x0 = x_ini;   % use your original initial state
% 
% % Simulate
% [t_nonlinear, X_nonlinear] = ode45(@(t,x) dyn(t,x,u_fun(t)), tspan, x0);

% DISCRETIZED SYSTEM
Ts = 0.004;

sys_c = ss(A, B, C, D);

sys_d = c2d(sys_c, Ts, 'zoh');

[Ad, Bd, Cd, Dd] = ssdata(sys_d);

N = floor(Tend/Ts) + 1;  
t_discrete = (0:N-1)*Ts;

X_discrete = zeros(N, size(Ad,1));
U_discrete = zeros(N, 1);

max_ex = 1/(1^2);
max_ev = 1/(1^2);
max_ue = 1;
Q = diag([max_ex,max_ex,max_ev,max_ev]);
% Q = diag([1,1,10,10]);
R = 1/(max_ue^2);

K = dlqr(Ad,Bd,Q,R)


x_tilde = x_ini- x_e;


%x_ref_tilde = x_ref - x_e

u_fun = @(x_tilde,t) -K*(x_tilde);  

for k = 1:N
    t_curr = (k-1)*Ts;       

    X_discrete(k,:) = x_tilde +x_e;         
    
    u = u_fun(x_tilde, t_curr);    
    U_discrete(k) = u +u_e;         
    
    x_tilde = Ad*x_tilde + Bd*u;          

end

%X_discrete = X_discrete + repmat(x_e', size(X_discrete,1), 1);

maxinput = max(abs(U_discrete));

figure(1);
  subplot(2,1,1)
  plot(tout(1:length(alpha)), alpha(:,1), 'LineWidth', 1.6); 
  hold on
  plot(tout(1:length(alpha)), alpha(:,3), 'LineWidth', 1.6);
  legend('\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles')
  % ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(tout(1:length(alpha)), alpha(:,2), 'LineWidth', 1.6); 
  hold on
  plot(tout(1:length(alpha)), alpha(:,4), 'LineWidth', 1.6);
  legend('d\alpha_0','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities')

 figure(2)
 subplot(2,1,1)
  plot(t, X(:,1), 'LineWidth', 1.6); hold on
  plot(t, X(:,2), 'LineWidth', 1.6);
  legend('\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles (Linearized Model)')
  % ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(t, X(:,3), 'LineWidth', 1.6); hold on
  plot(t, X(:,4), 'LineWidth', 1.6);
  legend('d\alpha_1','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities (Linearized Model)')

 % figure(3);
 %  subplot(2,1,1)
 %  plot(t_nonlinear, X_nonlinear(:,1), 'LineWidth', 1.6); hold on
 %  plot(t_nonlinear, X_nonlinear(:,2), 'LineWidth', 1.6);
 %  legend('\alpha_0','\alpha_2')
 %  ylabel('Angle [rad]')
 %  title('Pendubot Angles')
 %  % ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE
 % 
 %  subplot(2,1,2)
 %  plot(t_nonlinear, X_nonlinear(:,3), 'LineWidth', 1.6); hold on
 %  plot(t_nonlinear, X_nonlinear(:,4), 'LineWidth', 1.6);
 %  legend('d\alpha_0','d\alpha_2')
 %  ylabel('Angular Velocity [rad/s]')
 %  xlabel('Time [s]')
 %  title('Pendubot Angular Velocities')

 figure(3)
 figure(3);
  subplot(2,1,1)
  plot(t_discrete, ones(length(X_discrete)).*pi, 'LineWidth', 1.6); hold on 
  % plot(t_discrete, U_discrete, 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,1), 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,2), 'LineWidth', 1.6);
  legend('u','\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles (Discretized Model)')
  % ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(t_discrete, X_discrete(:,3), 'LineWidth', 1.6); hold on
  plot(t_discrete, X_discrete(:,4), 'LineWidth', 1.6);
  legend('d\alpha_1','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities (Discretized Model)')
  % ylim([scale_velocities_low scale_velocities_high])   % <-- SET SCALE

  figure(4);
  subplot(2,1,1)
  plot(M(:,1), M(:,4), 'LineWidth', 1.6); 
  hold on
  plot(M(:,1), M(:,6), 'LineWidth', 1.6);
  legend('\alpha_1','\alpha_2')
  ylabel('Angle [rad]')
  title('Pendubot Angles')
  % ylim([scale_angles_low scale_angles_high])    % <-- SET SCALE

  subplot(2,1,2)
  plot(M(:,1), M(:,5), 'LineWidth', 1.6); 
  hold on
  plot(M(:,1), M(:,7), 'LineWidth', 1.6);
  legend('d\alpha_0','d\alpha_2')
  ylabel('Angular Velocity [rad/s]')
  xlabel('Time [s]')
  title('Pendubot Angular Velocities')