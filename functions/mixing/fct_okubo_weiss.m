function [Q,sigma_w,r_OW] = fct_okubo_weiss(model,w)
% Compute and plot okubo weiss cirterion
%

%% Velocity gradient
grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% Mx My 2(grad) 2(w 2D) N_t

% Strain
sigma_s = grad_w(:,:,1,2,:) + grad_w(:,:,2,1,:);
sigma_n = grad_w(:,:,1,1,:) - grad_w(:,:,2,2,:);
vector_sigma_w = sigma_s + 1i * sigma_n;
vector_sigma_w = permute(vector_sigma_w, [ 1 2 3 5 4]); % Mx My 1 N_t
clear sigma_s sigma_n
% sigma_s = real(vector_sigma_w);
% sigma_n = imag(vector_sigma_w);

% sigma_w = sqrt(sigma_n.^2 + sigma_s.^2);
sigma_w = abs(vector_sigma_w);
sigma_w2 = sigma_w.^2;

% Rotation
vort = grad_w(:,:,1,2,:) - grad_w(:,:,2,1,:);
vort = permute(vort, [ 1 2 3 5 4]);

%% Criterion Q
Q = 1/4 * ( sigma_w2 - vort.^2);
% Q = 1/4 * ( sigma_n.^2 + sigma_s.^2 - vort.^2);
if nargout >2
    r_OW = vort./sigma_w;
end

%% Plots
x=model.grid.x;
y=model.grid.y;
figure;imagesc(x,y,...
    Q(:,:,:,model.plot.day_plot)');
axis xy;axis equal
if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
    caxis([-1 1]*2e-9);
end
figure;imagesc(x,y,(Q(:,:,:,model.plot.day_plot)<0)');
axis xy;axis equal

%% Criterion Okubo Weiss
r_OW = (vort)./sigma_w;
r2_OW=r_OW.^2;

%% Plots
day_plot = model.plot.day_plot;

x=model.grid.x;
y=model.grid.y;
figure;imagesc(x,y,(r2_OW(:,:,:,day_plot))');
colorbar;axis xy;axis equal
title('OW r2')
caxis([0 1]*2);
colormap('winter')

figure;imagesc(x,y,(1./r2_OW(:,:,:,day_plot))');
colorbar;axis xy;axis equal
title('OW 1/r2')
caxis([0 10]);
colormap('winter')

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'OW_r2_day' num2str(day_plot) '.eps']);


% % transport_bareer =  log(abs(1./Q)) ;
% transport_bareer =  1./Q ;
tol = 0.5;
transport_bareer =  ( ((1-tol) < r2_OW) & (r2_OW < (1 + tol)) ) ;

figure;imagesc(x,y,transport_bareer(:,:,:,day_plot)');axis xy;axis equal
% caxis([0 1e1]);
title('transport bareer')
colormap('winter')


tol = 1e-1;
state_OW = 0 * (r2_OW < (1 - tol)) ...
    + 1 * ( ((1-tol) < r2_OW) & (r2_OW < (1 + tol)) ) ...
    + 2 * (r2_OW > (1 +tol));

figure;imagesc(x,y,state_OW(:,:,:,day_plot)');axis xy;axis equal
caxis([0 2]);
title('OW states')
colormap('winter')

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'OW_states_day' num2str(day_plot) '.eps']);


