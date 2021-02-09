%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate the deterministic and stochastic equation of the angle between
%%% tracer gradient and compressive strain direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pameters

% Physical parameters
r = 1.2;
% alpha2 = 0;
alpha2 = 1e-1;
alpha = sqrt(alpha2);

% Simulation parameters
T =  min ( [ 2/sqrt(abs(r^2-1)) 10/r ] ) * 10; % temps max
dt = min (  [ min( [ 2/sqrt(abs(r^2-1)) 1/r ] ) ...
             (2*pi/alpha)^2                       ] ) /100;
% T = 3e2 * 1;
% dt = 1e-1 * 1;
nb_realiz = 1e3; % nb de particules / trajectoires
% T = 3e0 * 1;
% % T = 3e1 * 1;
% % % T = 1e1 * 1;
% dt = 1e-2 * 1;
% % dt = 1e-1 * 1;
% % % T = 1e1 * alpha2;
% % % dt = 1e-1 * alpha2;
theta_0 = 0.4; % angle initial
rho_0 = 1; % nombre d'onde / écart entre particule initial
nb_hist = 30; % discretization en abcisse pour les histogrammes
n_pdf = 100; % discretization en abcisse pour les plots ksdensity
max_plot_rho= 1/4; % max/min pour les abcisses pour les plots ksdensity
% max_plot_rho= T/4; % max/min pour les abcisses pour les plots ksdensity
%max_plot_rho=5 *T/10;

% Pre-treatement
time = 0:dt:T;
N =length(time); % nd pas de temps
alpha = sqrt(alpha2);

% Sampling of noise
dBt = sqrt(dt) * randn(N,nb_realiz,2);
% dBt = sqrt(dt) * randn(N,2);

%% Time integration
theta= nan(N,nb_realiz);
ln_rho= nan(N,nb_realiz);
theta(1,:)=theta_0;
ln_rho(1,:)=log(rho_0);

XI_theta = linspace( - pi, pi,n_pdf); % discretization en abcisse pour les plots ksdensity

XI_rho = linspace( - max_plot_rho, max_plot_rho,n_pdf);% discretization en abcisse pour les plots ksdensity

XI_sin = linspace( - 1, 1,n_pdf);% discretization en abcisse pour les plots ksdensity

F_theta= nan(N,n_pdf);
F_rho= nan(N,n_pdf);
F_sin= nan(N,n_pdf);
F_theta(1,:)=ksdensity(theta(1,:),XI_theta);
F_rho(1,:)=ksdensity((ln_rho(+1,:)),XI_rho);
F_sin(1,:)=ksdensity(-sin(theta(+1,:)),XI_sin);

for k_t=1:(N-1)
% for t_local = time(2:end)
    theta(k_t+1,:) = theta(k_t,:) ...
        + ( r - cos( theta(k_t,:) ) ) * dt ...
        + alpha * dBt(k_t,:,1) ;
%     ln_rho(k_t+1,:) = ln_rho(k_t,:) ...
%         + ln_rho(k_t,:)/2 .* ( ( alpha2/4 - sin( theta(k_t,:) ) ) * dt ...
%                          + alpha * dBt(k_t,:,2) );
    ln_rho(k_t+1,:) = ln_rho(k_t,:) ...
        +  ( ( alpha2/12 - sin( theta(k_t,:) ) ) * dt ...
                      + alpha/sqrt(12) * dBt(k_t,:,2) );
    ln_rho(k_t+1, ln_rho(k_t+1,:)<0 ) = 0;
    F_theta(k_t+1,:)=ksdensity(theta(k_t+1,:),XI_theta);
    F_rho(k_t+1,:)=ksdensity((ln_rho(k_t+1,:))/time(k_t+1),XI_rho);
    F_sin(k_t+1,:)=ksdensity(-sin(theta(k_t+1,:)),XI_sin);
end
theta = mod(theta+pi,2*pi)-pi;

%% Plots
% close
% figure(1)
% subplot(1,2,1)
% if abs(r) <1
%     plot(time,-acos(r)*ones(N,1),'k')
% end
% if abs(r)==1
%     plot(time,(1-r)*pi/2*ones(N,1),'k')
% end
% hold on
% % plot(time,atan(theta/2))
% plot(time,theta)
% hold off
% ax=axis;
% delta_ax = ax(4)-ax(3);
% ax(3) = ax(3) - 0.1 * delta_ax;
% ax(4) = ax(4) + 0.1 * delta_ax;
% axis(ax);
% subplot(1,2,2)
% plot(time,ln_rho)

%% Plots pdf
figure(2)
subplot(1,3,1)
ksdensity(theta(end,:))
% hist(theta(end,:),nb_hist)
title('$\theta$',...
    'Interpreter','latex')
% hist(tan(theta(end,:)/2),nb_hist)
% title('$\tan(\theta/2)$',...
%     'Interpreter','latex')
ax=axis;
ax(1:2) =[-pi pi ];
axis(ax);
subplot(1,3,2)
ksdensity((ln_rho(end,:))/time(end))
% hist(log(ln_rho(end,:))/time(end),nb_hist)
title('$\ln(\rho)/t$',...
    'Interpreter','latex')
ax=axis;
ax(1:2) = max_plot_rho * [-1 1 ];
axis(ax);
subplot(1,3,3)
ksdensity(-sin(theta(end,:)))
axis xy;
title('$-\sin(\theta)$',...
    'Interpreter','latex')


%% Plots pdf - time
figure(3)
subplot(1,3,1)
imagesc(time,XI_theta,F_theta');
axis xy;
title('$\theta$',...
    'Interpreter','latex')
cax = caxis; cax(2)=cax(2)/100;caxis(cax);
subplot(1,3,2)
imagesc(time,XI_rho,F_rho');
axis xy;
title('$\ln(rho)/t$',...
    'Interpreter','latex')
subplot(1,3,3)
imagesc(time,XI_sin,F_sin');
axis xy;
title('$-\sin(\theta)$',...
    'Interpreter','latex')


%% Plots log-pdf - time
figure(4)
subplot(1,3,1)
imagesc(time,XI_theta,log(F_theta)');
axis xy;
title('$\theta$',...
    'Interpreter','latex')
subplot(1,3,2)
imagesc(time,XI_rho,log(F_rho)');
axis xy;
title('$\ln(\rho)/t$',...
    'Interpreter','latex')
subplot(1,3,3)
imagesc(time,XI_sin,log(F_sin)');
axis xy;
title('$-\sin(\theta)$',...
    'Interpreter','latex')
