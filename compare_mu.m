%%%% Varying Risk Aversion %%%%%%%%%
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


%%%%%%%%%%%%%%%% varying mu %%%%%%%%%%%%%%%%%%%%%%

for k = 1:151
    mu = -0.02 + 0.02*k; %mu1
    xi = 10.0; 
    [ctrl_t, ctrl_alpha, ctrl_h] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag);
    ctrl_drift = A + (B+ctrl_h).*ctrl_alpha;
    ctrl_int = trapz(ctrl_t, ctrl_drift);
    ctrl_mean(k) = exp(ctrl_int);
    ctrl_quad = trapz(ctrl_t, ctrl_alpha.*ctrl_alpha);
    ctrl_std(k) = ctrl_mean(k)*sqrt(exp(ctrl_quad) - 1);
end

for k = 1:len
    mu = -0.02 + 0.02*k; %mu1
    xi = 10.0; 
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

plot(close_std, close_mean, '-','LineWidth', 1, 'Color', [255, 127, 14]./255);
plot(ctrl_std, ctrl_mean, '-', 'LineWidth', 1, 'Color', [44, 160, 44]./255);

text(ctrl_std(151)-0.2, ctrl_mean(151), '$\mu_1=3 \rightarrow$', 'Interpreter', 'latex');
text(st_std(51)-0.2, st_mean(51), '$\mu_1=1 \rightarrow$', 'Interpreter', 'latex');
text(close_std(51), close_mean(51), '$\leftarrow \mu_1=1$', 'Interpreter', 'latex');
text(st_std(151), st_mean(151), '$\leftarrow \mu_1=3$', 'Interpreter', 'latex');
text(close_std(151)-0.2, close_mean(151), '$\mu_1=3 \rightarrow$', 'Interpreter', 'latex');

leg = legend('Open', 'Closed', 'CDAA');
set(leg, 'Interpreter', 'latex');
set(leg, 'Location', 'southeast');

xlabel('standard deviation');
ylabel('mean');
hold off


exportgraphics(gca, 'mu.pdf', 'Resolution', 300);

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
