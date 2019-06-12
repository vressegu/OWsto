%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate the deterministic and stochastic equation of the angle between
%%% tracer gradient and compressive strain direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init
%% Pameters

% Physical parameters
r = 0.9;
alpha2 = 0;
alpha2 = 1e-1;

% Simulation parameters
T =  min ( [ 2/sqrt(abs(r^2-1)) abs(10/r) ] ) * 10;
T_statio = T;
T = 10*T;
dt = min (  [ min( [ 2/sqrt(abs(r^2-1)) abs(1/r) ] ) ...
             (2*pi)^2/alpha2                       ] ) /100;
% dt=dt/10;
n_statio = ceil(T_statio/dt);
% T = min ( [ 2/sqrt(abs(r^2-1)) 10/r ] ) *5;
% T = 3e2 * 1;
% dt = 1e-1 * 1;
% T = 1e1 * alpha2;
% dt = 1e-1 * alpha2;
theta_0 = 0.4;
rho_0 = 1;

% Pre-treatement
time = 0:dt:T;
N =length(time);
alpha = sqrt(alpha2);

% Sampling of noise
dBt = sqrt(dt) * randn(N,2);

%% Time integration
theta= nan(N,1);
ln_rho= nan(N,1);
theta(1)=theta_0;
ln_rho(1)=log(rho_0);
ln_rho_res(1)=log(rho_0);
for k_t=1:(N-1)
% for t_local = time(2:end)
    theta(k_t+1) = theta(k_t) ...
        + ( r - cos( theta(k_t) ) ) * dt ...
        + alpha * dBt(k_t,1) ;
    ln_rho(k_t+1) = ln_rho(k_t) ...
        +  ( ( alpha2/12 - sin( theta(k_t) ) ) * dt ...
             + alpha/sqrt(12) * dBt(k_t,2) );
    ln_rho_res(k_t+1) = ln_rho_res(k_t) ...
        +  - sin( theta(k_t))  * dt ;
end
theta = mod(theta+pi,2*pi)-pi;

%% Remove transient
time(1:n_statio)=[];
ln_rho(1:n_statio)=[];
ln_rho_res(1:n_statio)=[];
theta(1:n_statio)=[];
% theta=sin(time);
N=N-n_statio;
hN=floor(N/2);
freq = 1/dt * 1/N * (0:(hN-1));

% theta = unwrap(theta);

%% Plots
close
figure(1)
subplot(2,2,1)
if abs(r) <1
    plot(time,-acos(r)*ones(N,1),'k')
end
if abs(r)==1
    plot(time,(1-r)*pi/2*ones(N,1),'k')
end
hold on
% plot(time,tan(theta/2))
% title('$\tan(\theta/2)$',...
%     'Interpreter','latex')
plot(time,theta)
title('$\theta$',...
    'Interpreter','latex')
hold off
ax=axis;
delta_ax = ax(4)-ax(3);
ax(3) = ax(3) - 0.1 * delta_ax;
ax(4) = ax(4) + 0.1 * delta_ax;
axis(ax);
subplot(2,2,2)
plot(time,(ln_rho_res))
% plot(time,(ln_rho))
% title('$\ln(\rho)$',...
title('$\ln(\rho_{res})$',...
    'Interpreter','latex')

subplot(2,2,3)
gamma=abs(fft(theta)).^2;
loglog(freq,gamma(1:hN))
% plot(freq,gamma(1:hN))
title('Spectrum $\theta$',...
    'Interpreter','latex')
% hold off
% ax=axis;
% delta_ax = ax(4)-ax(3);
% ax(3) = ax(3) - 0.1 * delta_ax;
% ax(4) = ax(4) + 0.1 * delta_ax;
% axis(ax);
subplot(2,2,4)
gamma=abs(fft(ln_rho_res)).^2;
% gamma=abs(fft(ln_rho)).^2;
loglog(freq,gamma(1:hN))
% plot(freq,gamma(1:hN))
% title('Spectrum $\ln(\rho)$',...
title('Spectrum $\ln(\rho_{res})$',...
    'Interpreter','latex')

