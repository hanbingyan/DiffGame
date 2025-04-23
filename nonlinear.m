% % Quadratic relation
A = 0.02;
B = 0.06/0.16;
C = 0;
D = 1.0;
Q = 0.0;
R = 0.0;

T = 1.0;
G = 1.0;
nu = 1.0;

len = 111;
st_std = zeros(1, len);
st_mean = zeros(1, len);

%%%%%%%%%%%%%%%% Open nonlinear relation %%%%%%%%%%%%%%%%%%%%%%

for k = 1:len
    mu = 0.05*k; %mu1
    xi = 100*mu^2;

    [st_t, st_alpha, st_h] = sol_state(A, B, C, D, Q, R, xi, T, G, nu, mu);
    st_drift = A + (B+st_h).*st_alpha;
    st_int = trapz(st_t, st_drift);
    st_mean(k) = exp(st_int);
    st_quad = trapz(st_t, st_alpha.*st_alpha);
    st_std(k) = st_mean(k)*sqrt(exp(st_quad) - 1);

end



len = 1001;
close_std = zeros(1, len);
close_mean = zeros(1, len);

%%%%%%%%%%%%%%%% Closed nonlinear relation %%%%%%%%%%%%%%%%%%%%%%

for k = 1:len
    mu = 0.05*k; %mu1
    xi = 10*mu^2;

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
leg = legend('Open', 'Closed');
set(leg, 'Interpreter', 'latex');
set(leg, 'Location', 'northwest');

xlabel('standard deviation');
ylabel('mean');
hold off

exportgraphics(gca, 'nonlinear.pdf', 'Resolution', 300)
