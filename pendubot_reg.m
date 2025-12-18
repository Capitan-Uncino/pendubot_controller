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
% % SWING UP
% % xi = x'+ [pi; pi; 0; 0];
% xi = [x(1)+pi; x(2)+pi; x(3); x(4)];
% 
% p1=0.0148;%kg m^2
% p2=0.0051;
% p3=0.0046;
% p4=0.1003;%kg m
% p5=0.0303;
% g=9.81;
% 
% k = 3.9621/8;%V (Nm)^â?»1
% 
% sym x1 x2 x3 x4 v1 k1 k2
% 
% % xs = [x1; x2; x3; x4];
% 
% % New state variables
% z1 = x1-alpha1r;
% eta1 = x2;
% z2 = x3;
% eta2 = x4;
% 
% w = [z1; z2; eta1; eta2];
% 
% % Equilibrium point
% x_e = [alpha1r; pi; 0; 0];
% 
% %M(alpha)
% m11 = p1 + p2 + 2*p3*cos(eta1-z1);  m12 = p2 + p3*cos(eta1-z1);
% m21 = p2 + p3*cos(eta1-z1);         m22 = p2;
% 
% %Vm(alpha, dalpha/dt)*dalpha/dt
% Vmda = p3*sin(eta1-z1)*[-eta2*z2 + z2^2 - eta2^2 ; z2^2];
% 
% %G(alpha)
% G = [p4*g*sin(z1); p5*g*sin(eta1)];
% 
% % Solve for equilibrium input
% eqns = [z2 == 0; z1 == 0; eta2 == 0; -1/m22 * (Vmda(2)+G(2)) - (m12/m22 * v1) == 0];
% v1_e = double(vpa(solve(subs(eqns, w, x_e ), v1), 5));
% 
% [kp, kd] = solve(k1 > 0, k2 > 0, v1_e == k1*(alpha1r-x1)-k2*x3, x, x_e);
% 
% V1 = kp*(alpha1r-x1)-kd*x3;
% 
% % Partial Feedback Linearization
% %M(alpha)
% m11 = p1 + p2 + 2*p3*cos(xi(2)-xi(1));  m12 = p2 + p3*cos(xi(2)-xi(1));
% m21 = p2 + p3*cos(xi(2)-xi(1));         m22 = p2;
% 
% %Vm(alpha, dalpha/dt)*dalpha/dt
% Vmda = p3*sin(xi(2)-xi(1))*[-xi(4)*xi(3) + xi(3)^2 - xi(4)^2 ; xi(3)^2];
% 
% %G(alpha)
% G = [p4*g*sin(xi(1)); p5*g*sin(xi(2))];
% 
% m11_hat = m11 - (m12*m21/m22);
% Vmda1_hat = Vmda(1) - (m12*Vmda(2)/m22);
% G1_hat = G(1) - (m12*G(2)/m22);
% 
% tau = m11_hat*V1 + Vmda1_hat + G1_hat;
% u = tau/k;

%========================================================
% STATE RECONSTRUCTION
xi = [x(1)+pi; x(2)+pi; x(3); x(4)];
alpha1 = xi(1);
alpha2 = xi(2);
dalpha1 = xi(3);
dalpha2 = xi(4);

% Physical parameters
p1 = 0.0148;   % kg m^2
p2 = 0.0051;
p3 = 0.0046;
p4 = 0.1003;   % kg m
p5 = 0.0303;
g  = 9.81;

k = 3.9621/8;  % motor constant (Nm/V)
u_max   = 8;

% Mass matrix
m11 = p1 + p2 + 2*p3*cos(alpha2-alpha1);
m12 = p2 + p3*cos(alpha2-alpha1);
m21 = m12;
m22 = p2;

% Coriolis / centrifugal terms
Vmda = p3*sin(alpha2-alpha1) * [-dalpha2*dalpha1 + dalpha1^2 - dalpha2^2; dalpha1^2];

% Gravity terms
G = [p4*g*sin(alpha1); p5*g*sin(alpha2)];

% -------------------------------------------------------
% Partial Feedback Linearization
% Reduced dynamics
m11_hat   = m11 - (m12*m21)/m22;
Vmda1_hat = Vmda(1) - (m12/m22)*Vmda(2);
G1_hat    = G(1)    - (m12/m22)*G(2);

% % alpha_e1 = alpha1r-xi(1);
% % alpha_e2 = alpha1r-xi(2);
% delta = 5;
% lim_inf = alpha1r - delta/2;
% lim_sup = alpha1r + delta/2;
% kp = 25;     % <-- sostituisci con i tuoi
% kd = 15;
% K = [-0.2165, -16.2534, -1.2892, -2.1864];
% % K = [-2, -15, -3, -5];
% 
% 
% if xi(1) > lim_inf && xi(1) < alpha1r
%     V1 = -kp*(xi(1) - delta/2) - kd*xi(3);
%     if xi(2) > lim_inf/2 && xi(2) < alpha1r
%         x_e = [alpha1r; pi; 0; 0];
%         u_e = -1.7205;
%         u = u_e - K*(xi-x_e);
%         u_max = 8;   % Volt (esempio)
%         u = max(min(u, u_max), -u_max);
%         return
%     elseif xi(2) > alpha1r && xi(2) < lim_sup
%         x_e = [alpha1r; pi; 0; 0];
%         u_e = -1.7205;
%         u = u_e - K*(xi-x_e);
%         u_max = 8;   % Volt (esempio)
%         u = max(min(u, u_max), -u_max);
%         return
%     end
% 
% elseif xi(1) > alpha1r && xi(1) < lim_sup
%     V1 = -kp*(xi(1) + delta/2) - kd*xi(3);
%     if xi(2) > lim_inf/2 && xi(2) < alpha1r
%         x_e = [alpha1r; pi; 0; 0];
%         u_e = -1.7205;
%         u = u_e - K*(xi-x_e);
%         u_max = 8;   % Volt (esempio)
%         u = max(min(u, u_max), -u_max);
%         return
%     elseif xi(2) > alpha1r && xi(2) < lim_sup
%     x_e = [alpha1r; pi; 0; 0];
%     u_e = -1.7205;
%     u = (u_e - K*(xi-x_e));
%     u_max = 8;   % Volt (esempio)
%     u = max(min(u, u_max), -u_max);
%     return
%     end
% else 
%     % Virtual control (PD on first joint)
%     V1 = -kp*(xi(1) - alpha1r) - kd*xi(3);
% end
% 
% % if x(2) > lim_inf && x(2) < alpha1r
% %         x_e = [alpha1r; pi; 0; 0];
% %         u_e = -1.7205;
% %         K = [-0.2165, -16.2534, -1.2892, -2.1864];
% %         u = u_e - K*(xi-x_e);
% %         u_max = 8;   % Volt (esempio)
% %         u = max(min(u, u_max), -u_max);
% %         return
% %     end
% 
% % Torque
% tau = m11_hat*V1 + Vmda1_hat + G1_hat;
% 
% % U = [p4*g*sin(alpha1r);
% %      p5*g*sin(alpha1r)];
% % tau_ref = U;
% 
% % Motor input
% u = tau / k;
% 
% u_max = 8;   % Volt (esempio)
% u = max(min(u, u_max), -u_max);

% ---------- LQR REGION ----------
alpha_tol  = 0.15;     % rad
dalpha_tol = 0.4;

use_LQR = abs(alpha1-alpha1r) < alpha_tol && ...
          abs(alpha2-pi)      < alpha_tol && ...
          abs(dalpha1)        < dalpha_tol && ...
          abs(dalpha2)        < dalpha_tol;

if use_LQR
    % Equilibrium
    x_e = [alpha1r; pi; 0; 0];
    u_e = -1.7205;

    % LQR gain
    K = [-0.2165, -16.2534, -1.2892, -2.1864];

    tau = u_e - K*(xi - x_e);

else
    % ---------- PFL + PD ----------
    kp = 25;
    kd = 15;

    % Virtual control
    V1 = -kp*(alpha1 - alpha1r) - kd*dalpha1;

    % Torque
    tau = m11_hat*V1 + Vmda1_hat + G1_hat;
end

%% ===== MOTOR INPUT =====
u = tau / k;
u = max(min(u, u_max), -u_max);

end