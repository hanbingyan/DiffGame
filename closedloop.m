function [x, alpha, h, y] = closedloop(A, B, xi, T, mu)
% Closed-loop control in Han, Pun, and Wong (Finance and Stochastics, 2021)
% % Note that xi in F&S is 1/xi in diff game work. 

r = A;
theta = B;
xi = 1/xi;
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
%%%% integral(s) = int^T_s 2*(r+...) du
% integral = zeros(1, T/0.01-1);
% integrand = 2*(r+(theta+h).*alpha)+alpha.^2;
% k = 0;
% for s = 0:0.01:T-0.02
%     k = k+1;
%     integral(k) = trapz(s:0.01:T-0.01, integrand(k:end)); 
% end
% % Intexp(s) = int^T_s exp(..)/xi*h^2 dv
% Intexp = zeros(1, T/0.01-1);
% k = 0;
% for s = 0:0.01:T-0.02
%     inttmp = zeros(1, int16((T-s)/0.01));
%     for v = s+0.01:0.01:T-0.01
%         inttmp(int16((v-s)/0.01)) = trapz(s:0.01:v, integrand(int16(s/0.01)+1:int16(v/0.01)+1));
%     end
%     Intexp(int16(s/0.01)+1) = trapz(s:0.01:T-0.01, inttmp.*h(int16(s/0.01)+1:end).^2/xi);
% end
% figure
% plot(0:0.01:T-0.02,exp(integral) - Intexp,'k');
% title('$P(s;t)$','Interpreter','latex');
% xlabel('Time');
% % export_fig Pst.eps -m9 -transparent
% figure
% plot(x,y(1,:),'-k','LineWidth',1);
% hold on
% plot(x,y(2,:),'--k','LineWidth',1);
% plot(x,y(3,:),':k','LineWidth',1);
% plot(x,Delta,'-.k','LineWidth',1);
% leg1 = legend('$M$', '$N$', '$\Gamma$', '$\Delta$');
% set(leg1, 'Interpreter', 'latex');
% set(leg1, 'Location', 'best');
% title('ODE solutions');
% xlabel('Time');
% hold off
% % export_fig ODEsolFull.eps -m9 -transparent
% figure
% plot(x, alpha,'-k','LineWidth',1);
% hold on
% plot(x, h, '-.k','LineWidth',1);
% leg2 = legend('$\alpha$','$h$');
% set(leg2, 'Interpreter', 'latex');
% set(leg2, 'Location', 'best');
% title('Equilibrium strategy');
% xlabel('Time');
% hold off
% export_fig Strategy.eps -m9 -transparent
% Mlow: lower bound for M
% Mlow = @(t) exp(2*r*(T-t))*(1-theta.^2/(2*r*xi))+theta.^2/(2*r*xi);
% time = 0:0.01:T;
% Mlowy = feval(Mlow, time);
% c = min(Mlowy)-0.01;
% % DeltaUp: upper bound for Delta
% DeltaUp = @(t) exp((2*r + theta.^2/4/xi/c + 3*sqrt(3)/8/sqrt(xi*c)*theta.^2)*(T-t)) - exp(2*r*(T-t)) - mu1*exp(r*(T-t));
% DeltaUpy = feval(DeltaUp, time);
% figure
% plot(time, DeltaUpy,'--k','LineWidth',1);
% hold on
% plot(x,Delta,'-.k','LineWidth',1);
% plot(time, Mlowy, ':k','LineWidth',1);
% plot(x,y(1,:),'-k','LineWidth',1);
% leg3 = legend('$\Delta$ upper bound', '$\Delta$','$M$ lower bound','$M$');
% set(leg3, 'Interpreter', 'latex');
% set(leg3, 'Location', 'best');
% xlabel('Time');
% hold off
% export_fig BoundFull.eps -m9 -transparent
end
