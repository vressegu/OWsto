function fct_mezic3(model,X,t,dot_red)
% This function plots Mezic's mixing criterion as well as personnal
% carcterisation of evolution of gradient
%

day = num2str(floor(t*model.advection.dt_adv/(3600*24)));
loc_colorbar = 'southoutside';
% colormap_ = 'default';
colormap_ = 'jet';

% Grid
x = model.grid.x_ref;
y = model.grid.y_ref;
if model.mirror
    My=model.grid.MX(2)/2;
    y=y(1:My);
else
    My=model.grid.MX(2);    
end

%% Mixing criterion

Xplot=reshape(X,[model.grid.MX,2]);
% Mx My 2
nabla_phi_x = 1/(2*model.grid.dX(1)) ...
    *(Xplot(3:end,2:end-1,:)-Xplot(1:end-2,2:end-1,:));
nabla_phi_y = 1/(2*model.grid.dX(2)) ...
    *(Xplot(2:end-1,3:end,:)-Xplot(2:end-1,1:end-2,:));
% Mx-2 My-2 2 2

nabla_phi(:,:,:,1) = nabla_phi_x;
nabla_phi(:,:,:,2) = nabla_phi_y;

nabla_phi_x = nabla_phi(:,:,1,:);
nabla_phi_y = nabla_phi(:,:,2,:);

iiix = nabla_phi_x >= model.grid.MX(1)/4;
%         figure;imagesc(iiix(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiix(:,:,:,2)');axis xy; axis equal; drawnow;
iiix = iiix(:);
if any(iiix)
    sX = size(nabla_phi_x);
    nabla_phi_x = nabla_phi_x(:);
    nabla_phi_x(iiix) = nabla_phi_x(iiix) - model.grid.MX(1)/2;
    nabla_phi_x = reshape(nabla_phi_x,sX);
end

iiix = nabla_phi_x < -model.grid.MX(1)/4;
%         figure;imagesc(iiix(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiix(:,:,:,2)');axis xy; axis equal; drawnow;
iiix = iiix(:);
if any(iiix)
    sX = size(nabla_phi_x);
    nabla_phi_x = nabla_phi_x(:);
    nabla_phi_x(iiix) = nabla_phi_x(iiix) + model.grid.MX(1)/2;
    nabla_phi_x = reshape(nabla_phi_x,sX);
end

iiiy = nabla_phi_y >= model.grid.MX(2)/4;
%         figure;imagesc(iiiy(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiiy(:,:,:,2)');axis xy; axis equal; drawnow;
iiiy = iiiy(:);
if any(iiiy)
    sX = size(nabla_phi_y);
    nabla_phi_y = nabla_phi_y(:);
    nabla_phi_y(iiiy) = nabla_phi_y(iiiy) - model.grid.MX(2)/2;
    nabla_phi_y = reshape(nabla_phi_y,sX);
end

iiiy = nabla_phi_y < -model.grid.MX(2)/4;
%         figure;imagesc(iiiy(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiiy(:,:,:,2)');axis xy; axis equal; drawnow;
iiiy = iiiy(:);
if any(iiiy)
    sX = size(nabla_phi_y);
    nabla_phi_y = nabla_phi_y(:);
    nabla_phi_y(iiiy) = nabla_phi_y(iiiy) + model.grid.MX(2)/2;
    nabla_phi_y = reshape(nabla_phi_y,sX);
end

nabla_phi(:,:,1,:) = nabla_phi_x;
nabla_phi(:,:,2,:) = nabla_phi_y;
%         nabla_phi(:,:,:,1) = nabla_phi_x;
%         nabla_phi(:,:,:,2) = nabla_phi_y;

nabla_meso_v = nabla_phi;
nabla_meso_v(:,:,1,1) = nabla_meso_v(:,:,1,1)-1;
nabla_meso_v(:,:,2,2) = nabla_meso_v(:,:,2,2)-1;
time_t = model.advection.dt_adv*t;
nabla_meso_v = 1/ time_t * nabla_meso_v;
tr_meso_v = nabla_meso_v(:,:,1,1) + nabla_meso_v(:,:,2,2);
det_meso_v = - 1/ time_t * tr_meso_v;
criterion1 = det_meso_v < 0;
criterion2 = (det_meso_v >= 0) & (det_meso_v <= 4/time_t^2);
criterion3 = det_meso_v > 4/time_t^2;
criterion = 1.8*criterion1 + 2.5 * criterion2 + 3 * criterion3;


%% Plot

width=4;
% width=12;
height=4;
figure5=figure(5);
set(figure5,'Units','inches', ...
    'Position',[10 20 width height], ...
    'PaperPositionMode','auto');

alpha2 = 1/2*sum(sum(nabla_phi.^2,4),3)-1;
alpha2 = alpha2 .* (alpha2 >eps );
%         alpha2 = alpha2 .* (alpha2 >=0 );
beta_alpha = sqrt((alpha2+2)./alpha2);
id_beta = isinf(beta_alpha(:));
beta_alpha(id_beta)=0;
beta_alpha=reshape(beta_alpha,[model.grid.MX-2]);

imagesc(x,y,criterion(:,2:My-1-2)')
axis xy
% colormap(colormap_)
% colorbar('location',loc_colorbar)
axis equal
title(' 3 states',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Clim',[0 4], ...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('y',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('x',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
axis([x(1) x(end) y(1) y(end)])

hold on;
plot(dot_red(1,:),dot_red(2,:),'.r');
% plot3(dot_red_savex(1,:),dot_red_savey(1,:),(1:(size(dot_red_savey,2))));
%         keyboard;

drawnow


eval( ['print -depsc ' model.folder.folder_simu '/mezic_state/' ...
    num2str(day) '.eps']);
%         keyboard;
