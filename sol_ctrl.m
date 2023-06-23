function [ctrl_t, ctrl_alpha, ctrl_h, y] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag)
% Helper function to solver ODEs for control-dependent ambiguity aversion
% flag: only one solution does not have zeros
[ctrl_t, y] = ode45(@(t, y) ODE_ctrl(t, y, A, B, D, Q, R, xi, flag), [0, T], [G/2; nu; mu; G; nu; mu], ...
    odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'Stats', 'on'));

% reverse time
y = flipud(y);
% figure
% plot(ctrl_t, y(:,1), '-', 'LineWidth',1);
% hold on
% plot(ctrl_t, y(:,2), '--', 'LineWidth',1);
% plot(ctrl_t, y(:,3), '-.', 'LineWidth',1);
% plot(ctrl_t, y(:,4), '-', 'LineWidth',1);
% plot(ctrl_t, y(:,5), '--', 'LineWidth',1);
% plot(ctrl_t, y(:,6), '-.', 'LineWidth',1);
% leg1 = legend('L', 'H', 'J', 'M', 'N', '\Gamma');
% set(leg1, 'Location', 'best');
% xlabel('Time');
% title('Control dependent ODE solutions');
% hold off

delta = y(:, 4) - y(:, 5) - y(:, 6);
E = 2*y(:, 1) - y(:, 2) - y(:, 3);

kappa = D^2*y(:, 4) + R;
beta = B*delta;
gamma = D^2*(delta.*E - E.^2)/xi;

if flag == 1
    ctrl_alpha = (-beta + sqrt(beta.^2 - 4*kappa.*gamma))/2./kappa;
elseif flag == -1
    ctrl_alpha = (-beta - sqrt(beta.^2 - 4*kappa.*gamma))/2./kappa;
end
    
ctrl_h = D*E./ctrl_alpha/xi;

% figure
% plot(ctrl_t, ctrl_alpha, '-', 'LineWidth',1);
% hold on
% plot(ctrl_t, ctrl_h, '--', 'LineWidth',1);
% leg2 = legend('\alpha', 'h');
% set(leg2, 'Location', 'best');
% xlabel('Time');
% title('Control dependent controls');
% hold off

% validate positive definite condition when Q=0
% step_len = T/(length(ctrl_h) - 1);
% t_grid = 0:step_len:T;
% P = G*exp(2*A*(T-t_grid));
% min(R - ctrl_h.^2*xi + D^2*P')
end

