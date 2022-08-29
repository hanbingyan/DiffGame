A = 0.02;
B = 0.06/0.16;
C = 0;
D = 1.0;
Q = 0.0;
R = 0.0;

T = 1.0;
G = 1.0;
nu = 1.0;


% only one solution does not have zeros.
flag = 1;


len = 10001;
ctrl_std = zeros(1, 151);
ctrl_mean = zeros(1, 151);
st_std = zeros(1, len);
st_mean = zeros(1, len);
close_std = zeros(1, len);
close_mean = zeros(1, len);


%%%%%%%%%%%%%%%% linear relation %%%%%%%%%%%%%%%%%%%%%%

for k = 1:151
    mu = -0.02 + 0.02*k; %mu1
    xi = 5.0 + 2.0*mu;
    [ctrl_t, ctrl_alpha, ctrl_h] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag);
    ctrl_drift = A + (B+ctrl_h).*ctrl_alpha;
    ctrl_int = trapz(ctrl_t, ctrl_drift);
    ctrl_mean(k) = exp(ctrl_int);
    ctrl_quad = trapz(ctrl_t, ctrl_alpha.*ctrl_alpha);
    ctrl_std(k) = ctrl_mean(k)*sqrt(exp(ctrl_quad) - 1);
end

for k = 1:len
    mu = -0.02 + 0.02*k; %mu1
%     xi = 2.0*mu;
    xi = 5.0 + 2.0*mu;

    [st_t, st_alpha, st_h] = sol_state(A, B, C, D, Q, R, xi, T, G, nu, mu);
    st_drift = A + (B+st_h).*st_alpha;
    st_int = trapz(st_t, st_drift);
    st_mean(k) = exp(st_int);
    st_quad = trapz(st_t, st_alpha.*st_alpha);
    st_std(k) = st_mean(k)*sqrt(exp(st_quad) - 1);

    [x, alpha, h] = closedloop(A, B, xi, T, mu);
    drift = A + (B+h).*alpha;
    intg = trapz(x, drift);
    close_mean(k) = exp(intg);
    close_quad = trapz(x, alpha.*alpha);
    close_std(k) = close_mean(k)*sqrt(exp(close_quad) - 1);
end

fig = figure;
plot(st_std, st_mean, '-', 'LineWidth', 1, 'Color', [31, 119, 180]./255);
hold on
text(ctrl_std(151)-0.36, ctrl_mean(151), '$\mu_1=3, \xi=11 \rightarrow$', 'Interpreter', 'latex');
plot(close_std, close_mean, '-','LineWidth', 1, 'Color', [255, 127, 14]./255);
plot(ctrl_std, ctrl_mean, '-', 'LineWidth', 1, 'Color', [44, 160, 44]./255);


text(close_std(1001)-0.38, close_mean(1001), '$\mu_1=20, \xi=45 \rightarrow$', 'Interpreter', 'latex');
text(st_std(501), st_mean(501), '$\leftarrow \mu_1=10, \xi=25$', 'Interpreter', 'latex');

text(close_std(10001)-0.43, close_mean(10001), '$\mu_1=200, \xi=405 \rightarrow$', 'Interpreter', 'latex');
text(st_std(10001)-0.42, st_mean(10001), '$\mu_1=200, \xi=405 \rightarrow$', 'Interpreter', 'latex');

leg = legend('Open', 'Closed', 'CDAA');
set(leg, 'Interpreter', 'latex');
set(leg, 'Location', 'southwest');

xlabel('standard deviation');
ylabel('mean');
hold off


exportgraphics(gca, 'linear.pdf', 'Resolution', 300)

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
