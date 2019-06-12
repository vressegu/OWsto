function [r,sigma_w] = fct_lapeyre(model,w)
% Compute and plot okubo weiss cirterion
%

% close all
x=model.grid.x;
y=model.grid.y;
% warning('some time removed')
% w(:,:,:,33:end)=[];

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

% phi = 1/2 * angle(vector_sigma_w);
% phi = unwrap(phi,[],4);


% x=model.grid.x;
% y=model.grid.y;
% figure;imagesc(x,y,(phi(:,:,:,model.plot.day_plot))');
% colorbar;axis xy;axis equal

% sigma_w = sqrt(sigma_n.^2 + sigma_s.^2);
sigma_w = abs(vector_sigma_w);
sigma_w2 = sigma_w.^2;
% mask =( sigma_w2 >  1e-2 * std(sigma_w2(:)) );
% % mask =( sigma_w2 > eps);

% Rotation
vort = grad_w(:,:,1,2,:) - grad_w(:,:,2,1,:);
vort = permute(vort, [ 1 2 3 5 4]);

%% Rotation speed of the strain axis

% D_phi_Dt = material_deriv_mat(phi,w,model.dt,model.grid.dX);
% 
% figure;imagesc(x,y,...
%     ( mask(:,:,:,model.plot.day_plot) .* D_phi_Dt(:,:,:,model.plot.day_plot) ...
%     )');
% colorbar;axis xy;axis equal

% D_sigma_s_Dt = ...
%     material_deriv_mat(sigma_s,w,model.dt,model.grid.dX);
% D_sigma_n_Dt = ...
%     material_deriv_mat(sigma_n,w,model.dt,model.grid.dX);
% D_phi_bis = sigma_s .* D_sigma_n_Dt - sigma_n .* D_sigma_s_Dt;
% D_phi_bis = D_phi_bis ./ (2*sigma_w2);

% figure;imagesc(x,y,...
%     ( mask(:,:,:,model.plot.day_plot).*D_phi_bis(:,:,:,model.plot.day_plot)...
%     )');
% colorbar;axis xy;axis equal
% caxis([-1 1]*0.5e-4);
% % caxis([-1 1]*1e-4);

if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
    warning('the time step is too large for the CFL');
%     error('the time step is too large for the CFL');
end

D_vector_sigma_w_Dt = ...
    material_deriv_mat(vector_sigma_w,w,model.dt,model.grid.dX);
D_phi_Dt = real( 1i * vector_sigma_w .* conj(D_vector_sigma_w_Dt) );
D_phi_Dt = D_phi_Dt ./ (2*sigma_w2);

figure;imagesc(x,y,...
    ( D_phi_Dt(:,:,:,model.plot.day_plot)...
    )');
colorbar;axis xy;axis equal
caxis([-1 1]*0.5e-4);
% caxis([-1 1]*1e-4);
title('rotation correction')

figure;imagesc(x,y,...
    (...
    (vort(:,:,:,model.plot.day_plot) + 2 * D_phi_Dt(:,:,:,model.plot.day_plot))...
    )');
colorbar;axis xy;axis equal
caxis([-1 1]*0.5e-4);
% caxis([-1 1]*1e-4);
title('effective rotation')

figure;imagesc(x,y,...
    (  ...
    (vort(:,:,:,model.plot.day_plot) )...
    )');
colorbar;axis xy;axis equal
caxis([-1 1]*0.5e-4);
% caxis([-1 1]*1e-4);
title('rotation')


% D_phi_Dt = D_phi;

%% Criterion r
r = (vort + 2 * D_phi_Dt)./sigma_w;
r2=r.^2;

%% Criterion Okubo Weiss
r_OW = (vort)./sigma_w;
r2_OW=r_OW.^2;

%% Plots
day_plot = model.plot.day_plot;

%% OW plots 

x=model.grid.x;
y=model.grid.y;
figure;imagesc(x,y,(r2_OW(:,:,:,day_plot))');
colorbar;axis xy;axis equal
title('OW r2')
caxis([0 1]*2);

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'OW_r2_day' num2str(day_plot) '.eps']);


tol = 1e-1;
state_OW = 0 * (r2_OW < (1 - tol)) ...
    + 1 * ( ((1-tol) < r2_OW) & (r2_OW < (1 + tol)) ) ...
    + 2 * (r2_OW > (1 +tol));

figure;imagesc(x,y,state_OW(:,:,:,day_plot)');axis xy;axis equal
caxis([0 2]);
title('OW states')

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'OW_states_day' num2str(day_plot) '.eps']);


%% Lapeyre plots

figure;imagesc(x,y,(r2(:,:,:,day_plot))');
colorbar;axis xy;axis equal
title('Lapeyre r2')
caxis([0 1]*2);

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'Lapeyre_r2_day' num2str(day_plot) '.eps']);

% if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
%     cax = caxis;cax(1)=-cax(2);caxis(cax);
%     caxis([-1 1]*30);
% end
tol = 1e-1;
state = 0 * (r2 < (1 - tol)) ...
    + 1 * ( ((1-tol) < r2) & (r2 < (1 + tol)) ) ...
    + 2 * (r2 > (1 +tol));

figure;imagesc(x,y,state(:,:,:,day_plot)');axis xy;axis equal
caxis([0 2]);
title('Lapeyre states')

drawnow
eval( ['print -depsc ' model.folder.folder_simu '/' ...
    'Lapeyre_states_day' num2str(day_plot) '.eps']);

% figure;imagesc(x,y,(r2>1)');axis xy;axis equal

% keyboard

