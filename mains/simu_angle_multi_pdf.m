%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate the deterministic and stochastic equation of the angle between
%%% tracer gradient and compressive strain direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pameters

% Physical parameters
r = -1;


for r = -2:0.1:2
    
    alpha2 = 1e0;
    
    % Simulation parameters
    nb_realiz = 3e3;
    T = 1e2 * 1;
    dt = 1e-1 * 1;
    % T = 1e1 * alpha2;
    % dt = 1e-1 * alpha2;
    theta_0 = 0.4;
    rho_0 = 1;
    nb_hist = 30;
    
    % Pre-treatement
    time = 0:dt:T;
    N =length(time);
    alpha = sqrt(alpha2);
    
    % Sampling of noise
    dBt = sqrt(dt) * randn(N,nb_realiz,2);
    % dBt = sqrt(dt) * randn(N,2);
    
    %% Time integration
    theta= nan(N,nb_realiz);
    rho= nan(N,nb_realiz);
    theta(1,:)=theta_0;
    rho(1,:)=rho_0;
    for k_t=1:(N-1)
        % for t_local = time(2:end)
        theta(k_t+1,:) = theta(k_t,:) ...
            + ( r - cos( theta(k_t,:) ) ) * dt ...
            + alpha * dBt(k_t,:,1) ;
        rho(k_t+1,:) = rho(k_t,:) ...
            + rho(k_t,:)/2 .* ( ( alpha2/4 - sin( theta(k_t,:) ) ) * dt ...
            + alpha * dBt(k_t,:,2) );
    end
    theta = mod(theta+pi,2*pi)-pi;

end


%% Plots pdf
figure(2)
subplot(1,2,1)
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
subplot(1,2,2)
ksdensity(log(rho(end,:)))
% hist(log(rho(end,:)),nb_hist)
title('$\ln(\rho)$',...
    'Interpreter','latex')
ax=axis;
ax(1:2) = 50* [1 1 ];
axis(ax);