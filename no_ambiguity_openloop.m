% % Check the mean and variance go to infinity when mu -> Inf
A = 0.02;
B = 0.06/0.16;
T = 1.0;
mu = 1e50;

[st_t, y] = ode45(@(t, y) openloop(t, y, A, B), [0, T], [1; 1; mu], ...
    odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Stats', 'on'));

y = flipud(y);
delta = y(:, 1) - y(:, 2) - y(:, 3);
st_alpha = -B*delta./y(:, 1);

drift = A + B.*st_alpha;
intg = trapz(st_t, drift);
open_mean = exp(intg);
open_quad = trapz(st_t, st_alpha.*st_alpha);
open_std = open_mean*sqrt(exp(open_quad) - 1);

plot(st_alpha);

function dydt = openloop(t, y, A, B)
% The ODE system for state-dependent ambiguity aversion
delta = y(1) - y(2) - y(3);
alpha = -B*delta/y(1);

dydt = [(2*A + B*alpha)*y(1); ...
    (2*A + B*alpha)*y(2); ...
    A*y(3)];
end