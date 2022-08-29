function dydt = ODE_ctrl(t, y, A, B, D, Q, R, xi, flag)
% ODE system for control-dependent ambiguity aversion
% l is the no. of assets and can only be 1
% xi is chosen as Phi = xi*l*u^2
% flag controls positive or negative roots
l = 1;
E = 2*y(1) - y(2) - y(3);
delta = y(4) - y(5) - y(6);
kappa = D^2*y(4) + R;
beta = B*delta;
gamma = D^2*(delta*E - E^2)/xi/l;
if flag == 1
    alpha = (-beta + sqrt(beta^2 - 4*kappa*gamma))/2/kappa;
elseif flag == -1
     alpha = (-beta - sqrt(beta^2 - 4*kappa*gamma))/2/kappa;
end
dydt = [(2*A + 2*B*alpha + 2*D^2*E/xi + D^2*alpha^2)*y(1) + 0.5*Q + 0.5*R*alpha^2 - 0.5*D^2*E^2/xi; ...
    (2*A + 2*B*alpha + 2*D^2*E/xi)*y(2); ...
    (A + B*alpha + D^2*E/xi)*y(3); ...
    (2*A + B*alpha + D^2*E/xi)*y(4) + Q; ...
    (2*A + B*alpha + D^2*E/xi)*y(5); ...
    A*y(6)];
end