function Q = fct_okubo_weiss_old(model,w)
% Compute and plot okubo weiss cirterion
%

%% Velocity gradient
grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% w = permute( w,[ 1 2 4 3]);
% grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% grad_w = gradient_perso(model.grid, w);
sigma_n = grad_w(:,:,1,1) - grad_w(:,:,2,2);
sigma_s = grad_w(:,:,1,2) + grad_w(:,:,2,1);
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);

%% Criterion Q
Q = 1/4 * ( sigma_n.^2 + sigma_s.^2 - vort.^2);

%% Plots
x=model.grid.x;
y=model.grid.y;
figure;imagesc(x(3:end-2),y(3:end-2),Q(3:end-2,3:end-2)');axis xy;axis equal
if strcmp(model.type_data,'Gula') && ~model.filtering.smoothing
    caxis([-1 1]*2e-9);
end
figure;imagesc(x,y,(Q<0)');axis xy;axis equal
