% % Check the mean and variance go to infinity when mu -> Inf
% % Note that xi in F&S is 1/xi in diff game work. 
A = 0.02;
B = 0.06/0.16;
T = 1.0;
mu = 10;
r = A;
theta = B;
xi = 0.0;
mu1 = mu;
f = @(u,y) [(2*r +(-2*(y(1)-y(2)-y(3))*y(1)*theta.^2+(y(1)-y(2)-y(3)).^2*theta.^2)/(xi*(y(1)-y(2)-y(3)).^2+y(1)).^2)*y(1) - xi*(y(1)-y(2)-y(3)).^4*theta.^2/(xi*(y(1)-y(2)-y(3)).^2+y(1)).^2;
    2*(r - y(1)*(y(1)-y(2)-y(3))*theta.^2/(xi*(y(1)-y(2)-y(3)).^2+y(1)).^2)*y(2);
    (r - y(1)*(y(1)-y(2)-y(3))*theta.^2/(xi*(y(1)-y(2)-y(3)).^2+y(1)).^2)*y(3)];
sol = ode45(f, [0 T], [1; 1; mu1]);
x = 0:0.01: T-0.01;
y = deval(sol, x);
% Invert the time from u=T-t to t;
y = fliplr(y);   
Delta = y(1,:)-y(2,:)-y(3,:);
% Investment strategy alpha and drift h
alpha = - Delta./(xi*Delta.^2 + y(1,:))*theta;
h = - Delta.^2*theta*xi./(xi*Delta.^2 + y(1,:));

drift = A + (B+h).*alpha;
intg = trapz(x, drift);
close_mean = exp(intg);
close_quad = trapz(x, alpha.*alpha);
close_std = close_mean*sqrt(exp(close_quad) - 1);

% figure
% plot(close_std);