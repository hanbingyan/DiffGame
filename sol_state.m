function [st_t, st_alpha, st_h] = sol_state(A, B, C, D, Q, R, xi, T, G, nu, mu)

[st_t, y] = ode45(@(t, y) ODE_state(t, y, A, B, C, D, Q, R, xi), [0, T], [G/2; nu; mu; G; nu; mu], ...
    odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Stats', 'on'));

% reverse time
y = flipud(y);
% figure
% plot(t, y(:,1), '-', 'LineWidth',1);
% hold on
% plot(t, y(:,2), '--', 'LineWidth',1);
% plot(t, y(:,3), '-.', 'LineWidth',1);
% plot(t, y(:,4), '-', 'LineWidth',1);
% plot(t, y(:,5), '--', 'LineWidth',1);
% plot(t, y(:,6), '-.', 'LineWidth',1);
% leg1 = legend('L', 'H', 'F', 'M', 'N', '\Gamma');
% set(leg1, 'Location', 'best');
% xlabel('Time');
% title('State dependent ODE solutions');
% hold off

delta = y(:, 4) - y(:, 5) - y(:, 6);
E = 2*y(:, 1) - y(:, 2) - y(:, 3);
st_alpha = -(D*C*E.*delta/xi + B*delta + D*C*y(:, 4))./(R + y(:, 4)*D^2 + delta.*E*D^2/xi);
st_h = (C + D*st_alpha).*E/xi;

% figure
% plot(t, st_alpha, '-', 'LineWidth',1);
% hold on
% plot(t, st_h, '--', 'LineWidth',1);
% leg2 = legend('\alpha', 'h');
% set(leg2, 'Location', 'best');
% xlabel('Time');
% title('State dependent controls');
% hold off
end
