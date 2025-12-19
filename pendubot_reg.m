function u=pendubot_reg(alpha1r,x)
% C1================================
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


%u=3*((alpha1r-pi)-x(1))-0.3*x(3);
x_e = [alpha1r; pi; 0; 0];

x_r = [pi; pi; 0; 0];
u_e = -1.7205;
%K = [0.0690   -9.8058   -0.7260   -1.2034];   %possibile c2 1 1 10 10 100
%K = [-0.0023  -10.9877   -0.8320   -1.3835];   %possibile c1
%K = [0.1042  -19.0514   -1.4656   -2.6460];   %possibile c1
%K = [0.0888  -13.0422   -0.9831   -1.7121];   %possibile c1 1 1 10 10 20
K = [0.0644   -9.2012   -0.6783   -1.1075]; % 1 1 10 10 200
%K = [0.2165, -16.2534, -1.2892, -2.1864];

u = u_e - K*(x'+x_r-x_e);
% u = u_e - K*(x-x_e);
% u = 0;

% %========================================================
% % % SWING UP
% %========================================================
% % STATE RECONSTRUCTION
% xi = [x(1)+pi; x(2)+pi; x(3); x(4)];
% alpha1 = xi(1);
% alpha2 = xi(2);
% dalpha1 = xi(3);
% dalpha2 = xi(4);
% 
% % Physical parameters
% p1 = 0.0148;   % kg m^2
% p2 = 0.0051;
% p3 = 0.0046;
% p4 = 0.1003;   % kg m
% p5 = 0.0303;
% g  = 9.81;
% 
% k = 3.9621/8;  % motor constant (Nm/V)
% u_max   = 8;
% alpha2r = alpha1r;
% 
% %M(alpha)
% m11 = p1 + p2 + 2*p3*cos(alpha2-alpha1);    m12 = p2 + p3*cos(alpha2-alpha1);
% m21 = p2 + p3*cos(alpha2-alpha1);           m22 = p2;
% 
% %Vm(alpha, dalpha/dt)*dalpha
% Vmda = p3*sin(alpha2-alpha1) * [-dalpha2*dalpha1 + dalpha1^2 - dalpha2^2; dalpha1^2];
% 
% %G(alpha)
% G = [p4*g*sin(alpha1); p5*g*sin(alpha2)];
% 
% 
% 
% % % 2=====================================
% % % -------------------------------------------------------
% % % Partial Feedback Linearization
% % % Reduced dynamics
% % m11_hat   = m11 - (m12*m21)/m22;
% % Vmda1_hat = Vmda(1) - (m12/m22)*Vmda(2);
% % G1_hat    = G(1)    - (m12/m22)*G(2);
% % 
% % % Thresholds
% % delta = 5*(2*pi/60);
% % lim_inf1 = alpha1r - delta/2;
% % lim_sup1 = alpha1r + delta/2;
% % lim_inf2 = alpha2r - delta/2;
% % lim_sup2 = alpha2r + delta/2;
% % 
% % kp = 25;     
% % kd = 15;
% % 
% % % lqr parameters
% % K = [-0.2165, -16.2534, -1.2892, -2.1864];
% % % K = [-2, -15, -3, -5];
% % x_e = [alpha1r; alpha2r; 0; 0];
% % u_e = -1.7205;
% % 
% % 
% % if xi(1) > lim_inf1 && xi(1) < alpha1r
% %     V1 = -kp*(xi(1) - delta/2) - kd*xi(3);
% %     if xi(2) > lim_inf2 && xi(2) < alpha2r
% %         u = u_e - K*(xi-x_e);
% %         u = max(min(u, u_max), -u_max);
% %         return
% %     elseif xi(2) > alpha2r && xi(2) < lim_sup2
% %         u = u_e - K*(xi-x_e);
% %         u = max(min(u, u_max), -u_max);
% %         return
% %     end
% % 
% % elseif xi(1) > alpha1r && xi(1) < lim_sup1
% %     V1 = -kp*(xi(1) + delta/2) - kd*xi(3);
% %     if xi(2) > alpha2r && xi(2) < lim_sup2
% %         u = (u_e - K*(xi-x_e));
% %         u = max(min(u, u_max), -u_max);
% %         return
% %     elseif xi(2) > lim_inf2 && xi(2) < alpha2r
% %         u = u_e - K*(xi-x_e);
% %         u = max(min(u, u_max), -u_max);
% %         return
% %     end
% % else 
% %     % Virtual control (PD on first joint)
% %     V1 = -kp*(xi(1) - alpha1r) - kd*xi(3);
% % end
% % 
% % % if x(2) > lim_inf && x(2) < alpha1r
% % %         u = u_e - K*(xi-x_e);
% % %         u = max(min(u, u_max), -u_max);
% % %         return
% % %     end
% % 
% % % Torque
% % tau = m11_hat*V1 + Vmda1_hat + G1_hat;
% % 
% % % U = [p4*g*sin(alpha1r);
% % %      p5*g*sin(alpha1r)];
% % % tau_ref = U;
% % 
% % % Motor input
% % u = tau / k;
% % u = max(min(u, u_max), -u_max);
% 
% % 3========================================
% % % ---------- LQR REGION ----------
% % alpha_tol  = 0.15;     % rad
% % dalpha_tol = 0.4;
% % 
% % use_LQR = abs(alpha1-alpha1r) < alpha_tol && ...
% %           abs(alpha2-pi)      < alpha_tol && ...
% %           abs(dalpha1)        < dalpha_tol && ...
% %           abs(dalpha2)        < dalpha_tol;
% % 
% % if use_LQR
% %     % Equilibrium
% %     x_e = [alpha1r; pi; 0; 0];
% %     u_e = -1.7205;
% % 
% %     % LQR gain
% %     K = [-0.2165, -16.2534, -1.2892, -2.1864];
% % 
% %     tau = u_e - K*(xi - x_e);
% % 
% % else
% %     % ---------- PFL + PD ----------
% %     kp = 25;
% %     kd = 15;
% % 
% %     % Virtual control
% %     V1 = -kp*(alpha1 - alpha1r) - kd*dalpha1;
% % 
% %     % Torque
% %     tau = m11_hat*V1 + Vmda1_hat + G1_hat;
% % end
% % 
% % %% ===== MOTOR INPUT =====
% % u = tau / k;
% % u = max(min(u, u_max), -u_max);
% % 
% % end
% 
% % 4===========================
% %% =================================================
% %% ========== 1) SWING-UP ENERGETICO ================
% %% =================================================
% 
% % Energia del secondo link (target = upright)
% E  = 0.5*p2*dalpha2^2 + p5*g*(1 - cos(alpha2));
% E0 = 2*p5*g;   % energia desiderata (pi rad)
% 
% kE = 3;        % gain energetico
% 
% u_swing = kE * (E - E0) * sign(dalpha2*cos(alpha2));
% 
% %% =================================================
% %% ========== 2) PFL + PD ==========================
% %% =================================================
% 
% kp = 25;
% kd = 15;
% 
% V1 = -kp*(alpha1 - alpha1r) - kd*dalpha1;
% tau_pfl = m11_hat*V1 + Vmda1_hat + G1_hat;
% 
% %% =================================================
% %% ========== 3) LQR LOCALE ========================
% %% =================================================
% 
% x_e = [alpha1r; pi; 0; 0];
% u_e = -1.7205;
% K   = [-0.2165, -16.2534, -1.2892, -2.1864];
% 
% tau_lqr = u_e - K*(xi - x_e);
% 
% %% =================================================
% %% ========== BLENDING CONTINUO ====================
% %% =================================================
% 
% % distanza dall’equilibrio
% d = norm([alpha1-alpha1r;
%           alpha2-pi;
%           dalpha1;
%           dalpha2]);
% 
% % pesi smooth
% w_lqr   = exp(-4*d^2);
% w_pfl   = (1-w_lqr) * exp(-0.5*d^2);
% w_swing = 1 - w_lqr - w_pfl;
% 
% tau = w_swing*u_swing + w_pfl*tau_pfl + w_lqr*tau_lqr;
% 
% %% ===== MOTOR INPUT =====
% u = tau / k;
% u = max(min(u, u_max), -u_max);

end
