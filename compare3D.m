A = 0.05;
B = 1.0;
C = 0;
D = 1.0;
Q = 0.0;
R = 0.0;
xi = 5.0; % \xi_1 should be quite large
T = 1.0;
G = 1.0;
nu = 1.0;
mu = 1.0;

% only one solution does not have zeros.
flag = 1;

xi_len = 21;
mu_len = 11;
ctrl_std = zeros(xi_len, mu_len);
ctrl_mean = zeros(xi_len, mu_len);
st_std = zeros(xi_len, mu_len);
st_mean = zeros(xi_len, mu_len);
close_std = zeros(xi_len, mu_len);
close_mean = zeros(xi_len, mu_len);

[MU, XI] = meshgrid(1.0:0.1:2.0, 4.0:0.1:6.0);

for l = 1:xi_len 
    for k = 1:mu_len
        mu = MU(l, k);
        xi = XI(l, k);
        [ctrl_t, ctrl_alpha, ctrl_h] = sol_ctrl(A, B, D, Q, R, xi, T, G, nu, mu, flag);
        ctrl_drift = A + (B+ctrl_h).*ctrl_alpha;
        ctrl_int = trapz(ctrl_t, ctrl_drift);
        ctrl_mean(l, k) = exp(ctrl_int);
        ctrl_quad = trapz(ctrl_t, ctrl_alpha.*ctrl_alpha);
        ctrl_std(l, k) = ctrl_mean(l, k)*sqrt(exp(ctrl_quad) - 1);

        [st_t, st_alpha, st_h] = sol_state(A, B, C, D, Q, R, xi, T, G, nu, mu);
        st_drift = A + (B+st_h).*st_alpha;
        st_int = trapz(st_t, st_drift);
        st_mean(l, k) = exp(st_int);
        st_quad = trapz(st_t, st_alpha.*st_alpha);
        st_std(l, k) = st_mean(l, k)*sqrt(exp(st_quad) - 1);

        [x, alpha, h] = closedloop(A, B, xi, T, mu);
        drift = A + (B+h).*alpha;
        intg = trapz(x, drift);
        close_mean(l, k) = exp(intg);
        close_quad = trapz(x, alpha.*alpha);
        close_std(l, k) = close_mean(l, k)*sqrt(exp(close_quad) - 1);
    end
end

figure
f1 = surf(MU, XI, st_std);
xlabel('Mu')
ylabel('Xi')
zlabel('Mean')
hold on 
f2 = surf(MU, XI, close_std);
legend([f1, f2], {'Open', 'Closed'});
hold off

