A = 0.02;
B = 0.06/0.16;
C = 0;
D = 1.0;
Q = 0.0;
R = 0.0;

T = 0.1;
G = -0.01;
nu = -0.01;
mu = -1.0;

% only one solution does not have zeros.
% flag = 1;
% 
% 
% len = 10001;
% ctrl_std = zeros(1, 151);
% ctrl_mean = zeros(1, 151);
% st_std = zeros(1, len);
% st_mean = zeros(1, len);
% close_std = zeros(1, len);
% close_mean = zeros(1, len);

xi = -1000000000000.0; 
[st_t, st_alpha, st_h, y] = sol_state(A, B, xi, T, G, nu, mu);
st_drift = A + (B+st_h).*st_alpha;
st_int = trapz(st_t, st_drift);
st_mean = exp(st_int);
st_quad = trapz(st_t, st_alpha.*st_alpha);
st_std = st_mean*sqrt(exp(st_quad) - 1);




function [st_t, st_alpha, st_h, y] = sol_state(A, B, xi, T, G, nu, mu)

[st_t, y] = ode45(@(t, y) ODE_st(t, y, A, B, xi), [0, T], [G/2; nu; mu; G; nu; mu], ...
    odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'Stats', 'on'));

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
st_alpha = -B*delta./(y(:, 4) + delta.*E/xi);
st_h = st_alpha.*E/xi;

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



function dydt = ODE_st(t, y, A, B, xi)
delta = y(4) - y(5) - y(6);
E = 2*y(1) - y(2) - y(3);
alpha = -B*delta/(y(4) + delta*E/xi);
h = alpha*E/xi;
dydt = [(2*A + 2*(B + h)*alpha + alpha^2)*y(1) - h^2*xi/2; ...
    0.0; ...
    (A + (B + h)*alpha)*y(3); ...
    (2*A + h*alpha + B*alpha)*y(4) - h^2*xi; ...
    0.0; ...
    A*y(6)];
end


%%%%%%%%%%%%%%%% linear relation %%%%%%%%%%%%%%%%%%%%%%

% for k = 1:151
%     mu = -0.02 + 0.02*k; %mu1
%     xi = 5.0 + 2.0*mu;
%     [ctrl_t, ctrl_alpha, ctrl_h] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag);
%     ctrl_drift = A + (B+ctrl_h).*ctrl_alpha;
%     ctrl_int = trapz(ctrl_t, ctrl_drift);
%     ctrl_mean(k) = exp(ctrl_int);
%     ctrl_quad = trapz(ctrl_t, ctrl_alpha.*ctrl_alpha);
%     ctrl_std(k) = ctrl_mean(k)*sqrt(exp(ctrl_quad) - 1);
% end

% for k = 1:len
%     mu = -0.02 + 0.02*k; %mu1
%     xi = 2.0*mu;
% %     xi = 5.0 + 2.0*mu;
% 
%     [st_t, st_alpha, st_h] = sol_state(A, B, C, D, Q, R, xi, T, G, nu, mu);
%     st_drift = A + (B+st_h).*st_alpha;
%     st_int = trapz(st_t, st_drift);
%     st_mean(k) = exp(st_int);
%     st_quad = trapz(st_t, st_alpha.*st_alpha);
%     st_std(k) = st_mean(k)*sqrt(exp(st_quad) - 1);
% 
%     [x, alpha, h] = closedloop(A, B, xi, T, mu);
%     drift = A + (B+h).*alpha;
%     intg = trapz(x, drift);
%     close_mean(k) = exp(intg);
%     close_quad = trapz(x, alpha.*alpha);
%     close_std(k) = close_mean(k)*sqrt(exp(close_quad) - 1);
% end
% 
% fig = figure;
% plot(st_std, st_mean, '-', 'LineWidth', 1, 'Color', [31, 119, 180]./255);
% hold on
% text(ctrl_std(151)-0.35, ctrl_mean(151), '$\mu=3, \xi=11 \rightarrow$', 'Interpreter', 'latex');
% plot(close_std, close_mean, '-','LineWidth', 1, 'Color', [255, 127, 14]./255);
% plot(ctrl_std, ctrl_mean, '-', 'LineWidth', 1, 'Color', [44, 160, 44]./255);
% 
% 
% text(close_std(1001)-0.37, close_mean(1001), '$\mu=20, \xi=45 \rightarrow$', 'Interpreter', 'latex');
% text(st_std(501), st_mean(501), '$\leftarrow \mu=10, \xi=25$', 'Interpreter', 'latex');
% 
% text(close_std(10001)-0.41, close_mean(10001), '$\mu=200, \xi=405 \rightarrow$', 'Interpreter', 'latex');
% text(st_std(10001)-0.41, st_mean(10001), '$\mu=200, \xi=405 \rightarrow$', 'Interpreter', 'latex');
% 
% leg = legend('Open', 'Closed', 'CDAA');
% set(leg, 'Interpreter', 'latex');
% set(leg, 'Location', 'southwest');
% 
% xlabel('std');
% ylabel('mean');
% hold off
% 
% 
% exportgraphics(gca, 'linear.pdf', 'Resolution', 300)

%export_fig test.pdf -m2.5

% set(gca, 'FontName', 'Times');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [8.1 8.1]);
% set(gcf, 'PaperPosition', [0 0 8 8]);
% print(gcf, '-dpdf', '-r150', 'linear.pdf');

% set(fig, 'FontName', 'Times');
% set(fig, 'PaperUnits', 'inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig, '-dpdf', '-r150', 'linear.pdf');
