function [alpha2,state] = fct_mix_sto(model,r2,sigma_w)
% Compute and plot okubo weiss cirterion
%



% close all
x=model.grid.x;
y=model.grid.y;
day_plot = model.plot.day_plot;
% warning('some time removed')
% w(:,:,:,33:end)=[];

%% Parameters of sigma dBt

warning('we assume a to be known');
load([model.folder.folder_simu '/stat_temp_smooth.mat'],'a');

% a=300

model.sigma.a0 = a;
model.sigma.k_m = 6e-5;
% model.sigma.beta = 3;
model.sigma.beta = 2;
model.sigma.k_M = 2*pi / (10e3);
% model.sigma.k_M = pi / max(model.grid.dX);

model.sigma.c_k_M = fct_c_k_M(...
    model.sigma.a0,model.sigma.k_m,model.sigma.k_M,model.sigma.beta);
model.sigma.c_k_M

%% Criterion alpha
alpha2 = 12 * model.sigma.c_k_M ./ sigma_w ;

log_alpha2 = log10(alpha2);
min(log_alpha2(:))
% x=model.grid.x;
% y=model.grid.y;
figure;imagesc(x,y,log_alpha2(:,:,:,day_plot)');
colorbar;axis xy;axis equal
title('log(alpha2)')
% caxis([0 1]*2);

%% Interpolation

s =size(alpha2);
resh_alpha2 = alpha2(:);
resh_r = sqrt(r2(:));

mean_growth_rate_interp = nan(size(resh_alpha2));
m4_interp = nan(size(resh_alpha2));
state = nan(size(resh_alpha2));

load('save_r_alpha/average_time.mat', ...
    'alpha2_v','r_v','mean_sin_v_modif','m4_sin_v'); 

% Use symmetry
l_r =length(r_v);
l_r_usefull = (l_r-1)/2 +1;
v_r_usefull = l_r_usefull:l_r;
r_v = r_v(v_r_usefull);
mean_sin_v_modif = mean_sin_v_modif(v_r_usefull,:);
m4_sin_v = m4_sin_v(v_r_usefull,:);

% Large r
i_lr = resh_r > max(r_v);
resh_r(i_lr) = max(r_v);

% Very weak noise case
i_wn = resh_alpha2 < alpha2_v(2);
mean_growth_rate_interp(i_wn) = ...
    interp1( r_v, mean_sin_v_modif(:,1), ...
    resh_r(i_wn) ,'spline');
m4_interp(i_wn) = 0;
state(i_wn) = fct_choose_states_weak_noise( resh_r(i_wn) );

% Remove case of zero noise
alpha2_v(1) = [];
mean_sin_v_modif(:,1) =[];
m4_sin_v(:,1) =[];

% Larger noise
% grid for interpolation
[r_grid,log10_alpha2_grid] = meshgrid(r_v,log10(alpha2_v));
% Growth rate interpolation
mean_growth_rate_interp(~i_wn) = ...
    interp2( r_grid, log10_alpha2_grid, mean_sin_v_modif', ...
    resh_r(~i_wn), log10(resh_alpha2(~i_wn)), 'spline');
% Kurtosis interpolation
m4_interp(~i_wn) = ...
    interp2( r_grid, log10_alpha2_grid, m4_sin_v', ...
    resh_r(~i_wn), log10(resh_alpha2(~i_wn)), 'spline');      
% Decide of the states
state(~i_wn) = fct_choose_states( ...
    mean_growth_rate_interp(~i_wn),m4_interp(~i_wn), ...
    log10(resh_alpha2(~i_wn)));

state = reshape(state,s);
mean_growth_rate_interp = reshape(mean_growth_rate_interp,s);
m4_interp = reshape(m4_interp,s); 

mean_growth_rate_interp = sigma_w .* mean_growth_rate_interp;

%% Plots
day_plot = model.plot.day_plot;

% x=model.grid.x;
% y=model.grid.y;
% figure;imagesc(x,y,(r2_OW(:,:,:,day_plot))');
% colorbar;axis xy;axis equal
% title('OW r2')
% caxis([0 1]*2);
% 
% drawnow
% eval( ['print -depsc ' model.folder.folder_simu '/' ...
%     'OW_r2_day' num2str(day_plot) '.eps']);

figure;imagesc(x,y,(r2(:,:,:,day_plot))');
colorbar;axis xy;axis equal
title('Lapeyre r2')
caxis([0 1]*2);

drawnow
% eval( ['print -depsc ' model.folder.folder_simu '/' ...
%     'Lapeyre_r2_day' num2str(day_plot) '.eps']);

% if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
%     cax = caxis;cax(1)=-cax(2);caxis(cax);
%     caxis([-1 1]*30);
% end
tol = 1e-1;
Lapeyre_state = 0 * (r2 < (1 - tol)) ...
    + 1 * ( ((1-tol) < r2) & (r2 < (1 + tol)) ) ...
    + 2 * (r2 > (1 +tol));

figure;imagesc(x,y,Lapeyre_state(:,:,:,day_plot)');axis xy;axis equal
caxis([0 2]);
title('Lapeyre states')

drawnow
% eval( ['print -depsc ' model.folder.folder_simu '/' ...
%     'Lapeyre_states_day' num2str(day_plot) '.eps']);

figure;imagesc(x,y,state(:,:,:,day_plot)');axis xy;axis equal
% caxis([0 2]);
title('My states')

figure;imagesc(x,y,mean_growth_rate_interp(:,:,:,day_plot)');axis xy;axis equal
% caxis([0 2]);
title('Growth rate')

figure;imagesc(x,y,m4_interp(:,:,:,day_plot)');axis xy;axis equal
% caxis([0 2]);
title('Excess kurtosis')

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'my_states_day' num2str(day_plot) '.eps']);

% figure;imagesc(x,y,(r2>1)');axis xy;axis equal

% keyboard

