%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate the deterministic and stochastic equation of the angle between
%%% tracer gradient and compressive strain direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init;
%% Pameters

% Physical parameters
r_v = -2:0.1:2;
% r_v = -2:0.1:2;
alpha2_v = 10.^([-inf -3:0.5:2]);

s = [length(r_v) length(alpha2_v) ];

mean_theta_v = nan(s);
std_theta_v = nan(s);
m3_theta_v = nan(s);
m4_theta_v = nan(s);

mean_sin_v = nan(s);
std_sin_v = nan(s);
m3_sin_v = nan(s);
m4_sin_v = nan(s);

% mean_rho_v = nan(s);
% std_rho_v = nan(s);
% m3_rho_v = nan(s);
% m4_rho_v = nan(s);


for i=1:length(r_v)
    r = r_v(i);
    for j=1:length(alpha2_v)
        alpha2 = alpha2_v(j);
%         r
%         alpha2
        
        % Simulation parameters
        T =  min ( [ 2/sqrt(abs(r^2-1)) abs(10/r) ] ) * 10;
        dt = min (  [ min( [ 2/sqrt(abs(r^2-1)) abs(1/r) ] ) ...
            (2*pi)^2/alpha2                       ] ) /100;
        % T = 3e2 * 1;
        % dt = 1e-1 * 1;
        nb_realiz = 1e3;
        % T = 3e0 * 1;
        % % T = 3e1 * 1;
        % % % T = 1e1 * 1;
        % dt = 1e-2 * 1;
        % % dt = 1e-1 * 1;
        % % % T = 1e1 * alpha2;
        % % % dt = 1e-1 * alpha2;
        theta_0 = 0.4;
        rho_0 = 1;
        nb_hist = 30;
        n_pdf = 100;
        max_plot_rho= T/4;
        %max_plot_rho=5 *T/10;
        
        % Pre-treatement
        time = 0:dt:T;
        N =length(time);
        alpha = sqrt(alpha2);
        
        %% Sampling of noise
%         dBt = sqrt(dt) * randn(N,nb_realiz,2);
%         % dBt = sqrt(dt) * randn(N,2);
        
        %% Time integration
%         theta= nan(N,nb_realiz);
%         ln_rho= nan(N,nb_realiz);
%         theta(1,:)=theta_0;
%         ln_rho(1,:)=log(rho_0);
%         
% %         XI_theta = linspace( - pi, pi,n_pdf);
% %         XI_rho = linspace( - max_plot_rho, max_plot_rho,n_pdf);
% %         XI_sin = linspace( - 1, 1,n_pdf);
% %         F_theta= nan(N,n_pdf);
% %         F_rho= nan(N,n_pdf);
% %         F_sin= nan(N,n_pdf);
% %         F_theta(1,:)=ksdensity(theta(1,:),XI_theta);
% %         F_rho(1,:)=ksdensity((ln_rho(+1,:)),XI_rho);
% %         F_sin(1,:)=ksdensity(-sin(theta(+1,:)),XI_sin);
%         
%         for k_t=1:(N-1)
%             % for t_local = time(2:end)
%             theta(k_t+1,:) = theta(k_t,:) ...
%                 + ( r - cos( theta(k_t,:) ) ) * dt ...
%                 + alpha * dBt(k_t,:,1) ;
%             %     ln_rho(k_t+1,:) = ln_rho(k_t,:) ...
%             %         + ln_rho(k_t,:)/2 .* ( ( alpha2/4 - sin( theta(k_t,:) ) ) * dt ...
%             %                          + alpha * dBt(k_t,:,2) );
%             ln_rho(k_t+1,:) = ln_rho(k_t,:) ...
%                 +  ( ( alpha2/12 - sin( theta(k_t,:) ) ) * dt ...
%                 + alpha/sqrt(12) * dBt(k_t,:,2) );
%             ln_rho(k_t+1, ln_rho(k_t+1,:)<0 ) = 0;
% %             F_theta(k_t+1,:)=ksdensity(theta(k_t+1,:),XI_theta);
% %             F_rho(k_t+1,:)=ksdensity((ln_rho(k_t+1,:)),XI_rho);
% %             F_sin(k_t+1,:)=ksdensity(-sin(theta(k_t+1,:)),XI_sin);
%         end
%         theta = mod(theta+pi,2*pi)-pi;
        
        %% Save
%         theta_end = theta(end,:);
%         rho_end = theta(end,:);
        
%         save(['save_r_alpha/all_distrib_theta_rho_alpha2_' ...
%             num2str(alpha2) '_r_' num2str(r) '.mat']);
        
        load(['save_r_alpha/last_time_distrib_theta_rho_alpha2_' ...
            num2str(alpha2) '_r_' num2str(r) '.mat'],...
            'r','alpha2','T','dt','theta_0','rho_0','theta_end','rho_end');
        
        %% Moments
        mean_theta = mean(theta_end,2);
        std_theta = std(theta_end,0,2);
        center_norm_theta = (theta_end - mean_theta)/std_theta;
        m3 = mean(center_norm_theta.^3,2);
        m4 = mean(center_norm_theta.^4,2);
        
        m4(m4<3)=3;
        m4 = log(m4-3);
        
        mean_theta_v(i,j) = mean_theta;
        std_theta_v(i,j) = std_theta;
        m3_theta_v(i,j) = m3;
        m4_theta_v(i,j) = m4;
        
        mean_sin = mean(-sin(theta_end),2);
        std_sin = std(-sin(theta_end),0,2);
        center_norm_sin = (-sin(theta_end) - mean_sin)/std_sin;
        m3 = mean(center_norm_sin.^3,2);
        m4 = mean(center_norm_sin.^4,2);
        
        m4(m4<3)=3;
        m4 = log(m4-3);
        
        mean_sin_v(i,j) = mean_sin;
        std_sin_v(i,j) = std_sin;
        m3_sin_v(i,j) = m3;
        m4_sin_v(i,j) = m4;
        
        
    end
end

%% PLots

alpha2_v(1)=10^(-3.5);

alpha2_v = log(alpha2_v)/log(10);


figure(5)
subplot(2,2,1)
imagesc(r_v,alpha2_v,mean_theta_v');
axis xy;colorbar;caxis([-pi pi]);
subplot(2,2,2)
imagesc(r_v,alpha2_v,std_theta_v');
axis xy;colorbar;caxis([-pi pi]);
subplot(2,2,3)
imagesc(r_v,alpha2_v,m3_theta_v');axis xy;colorbar;
subplot(2,2,4)
imagesc(r_v,alpha2_v,m4_theta_v');axis xy;colorbar;

figure(6)
subplot(2,2,1)
imagesc(r_v,alpha2_v,mean_sin_v');
axis xy;colorbar;caxis([-1 1]);
subplot(2,2,2)
imagesc(r_v,alpha2_v,std_sin_v');
axis xy;colorbar;caxis([-1 1]);
subplot(2,2,3)
imagesc(r_v,alpha2_v,m3_sin_v');axis xy;colorbar;
subplot(2,2,4)
imagesc(r_v,alpha2_v,m4_sin_v');axis xy;colorbar;



