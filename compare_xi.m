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
ctrl_std = zeros(1, 1971);
ctrl_mean = zeros(1, 1971);
st_std = zeros(1, len);
st_mean = zeros(1, len);
close_std = zeros(1, len);
close_mean = zeros(1, len);


%%%%%%%%%%%%%%%% varying xi %%%%%%%%%%%%%%%%%%%%%%

for k = 1:1971
    xi = 2.9 + 0.1*k; %mu1
    mu = 2.0; % 1.0 + 3.0*mu;
    [ctrl_t, ctrl_alpha, ctrl_h] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag);
    ctrl_drift = A + (B+ctrl_h).*ctrl_alpha;
    ctrl_int = trapz(ctrl_t, ctrl_drift);
    ctrl_mean(k) = exp(ctrl_int);
    ctrl_quad = trapz(ctrl_t, ctrl_alpha.*ctrl_alpha);
    ctrl_std(k) = ctrl_mean(k)*sqrt(exp(ctrl_quad) - 1);
end

for k = 1:len
    xi = -0.02 + 0.02*k; %mu1
    mu = 2.0;

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
% fig.Position = [100 100 1000 800];
plot(st_std, st_mean, '-', 'LineWidth', 1, 'Color', [31, 119, 180]./255);
hold on

plot(close_std, close_mean, '-','LineWidth', 1, 'Color', [255, 127, 14]./255);
plot(ctrl_std, ctrl_mean, '-', 'LineWidth', 1, 'Color', [44, 160, 44]./255);

text(ctrl_std(1), ctrl_mean(1), '$\leftarrow \xi=3$', 'Interpreter', 'latex');
text(ctrl_std(1971), ctrl_mean(1971), '$\leftarrow \xi=200$', 'Interpreter', 'latex');
text(st_std(10001)-0.17, st_mean(10001), '$\xi=200 \rightarrow$', 'Interpreter', 'latex');
text(close_std(10001)-0.17, close_mean(10001), '$\xi=200 \rightarrow$', 'Interpreter', 'latex');

leg = legend('Open', 'Closed', 'CDAA');
set(leg, 'Interpreter', 'latex');
set(leg, 'Location', 'southwest');

xlabel('standard deviation');
ylabel('mean');
xlim([0 1.2]);
ylim([0.7 1.35]);
hold off

exportgraphics(gca, 'xi.pdf', 'Resolution', 300);

