function alpha2_m = fct_mezic(model,X,t)
% This function plots Mezic's mixing criterion as well as personnal
% carcterisation of evolution of gradient
%

loc_colorbar = 'southoutside';
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

width=12;
height=4;
figure5=figure(5);
set(figure5,'Units','inches', ...
    'Position',[10 20 width height], ...
    'PaperPositionMode','auto');

subplot(1,4,1)
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

meso_vort = time_t^2*(nabla_meso_v(:,:,2,1) ...
    - nabla_meso_v(:,:,1,2)).^2 /2 ;
plot_meso_vort=meso_vort;

subplot(1,4,2)
imagesc(x(2:end-1),y(2:end-1),((plot_meso_vort(2:end-1,2:My-3)))');
%         imagesc(x,y,(abs(plot_meso_vort(:,1:My-2)))');
axis xy
axis equal
colormap(colormap_)
colorbar('location',loc_colorbar)
%         title('$log \left (| \breve \omega | \right )$',...
title('$ \breve \omega ^2 \ t^2 /2 $',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
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
%         colorbar
%        title('log(mesovorticity) in mesoellictic areas')

%         plot_meso_vort=(1-criterion2).*meso_vort;
%         %        figure(6)
%         % %        imagesc(x,y,log(meso_vort(:,1:My-2))');
%         %        figure(7)
%         subplot(2,3,6)
%         imagesc(x,y,log(plot_meso_vort(:,1:My-2))');
%         axis equal
%         title('log(mesovorticity) in mesohyperbolic areas')

alpha2 = 1/2*sum(sum(nabla_phi.^2,4),3)-1;
alpha2 = alpha2 .* (alpha2 >eps );
%         alpha2 = alpha2 .* (alpha2 >=0 );
beta_alpha = sqrt((alpha2+2)./alpha2);
id_beta = isinf(beta_alpha(:));
beta_alpha(id_beta)=0;
beta_alpha=reshape(beta_alpha,[model.grid.MX-2]);
%        figure(7)
subplot(1,4,3)
%         subplot(2,2,3)
%         subplot(2,3,4)
imagesc(x(2:end-1),y(2:end-1),(alpha2((2:end-1),2:My-3))');
axis equal
axis xy
colormap(colormap_)
colorbar('location',loc_colorbar)
%         in mesohyperbolic areas
title('$\alpha^2$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
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
%         colorbar
%         title('$log(\alpha^2)$ in mesohyperbolic areas')
%        figure(8)
subplot(1,4,4)
%         subplot(2,2,4)
%         subplot(2,3,5)
s_beta_alpha = size(beta_alpha);
m_beta = (beta_alpha > max(beta_alpha(:))/2);
beta_alpha(m_beta)=0;
beta_alpha=reshape(beta_alpha,s_beta_alpha);
imagesc(x(2:end-1),y(2:end-1),(beta_alpha((2:end-1),2:My-3))');
set(gca,'CLim',[0 5]);
axis equal
axis xy
colormap(colormap_)
colorbar('location',loc_colorbar)
title('$\beta / \alpha$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
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
%         colorbar
%         title('$log(\beta / \alpha)$ in mesohyperbolic areas')

%         %         nabla_criterion(:,:,1) = abs ...
%         %             (criterion(2:end,1:end-1)-criterion(1:end-1,1:end-1));
%         %         nabla_criterion(:,:,2) = abs ...
%         %             (criterion(1:end-1,2:end)-criterion(1:end-1,1:end-1));% Mx-3 My-3 2
%         nabla_criterion(1:model.grid.MX(1)-3,:,1) = abs ...
%             (criterion(2:end,:)-criterion(1:end-1,:)) == 2;
%         nabla_criterion(model.grid.MX(1)-2,:,1) = false(1,model.grid.MX(2)-2);
%         nabla_criterion(2:end,:,1) = nabla_criterion(2:end,:,1) ...
%             | nabla_criterion(1:end-1,:,1);
%         nabla_criterion(:,1:model.grid.MX(2)-3,2) = abs ...
%             (criterion(:,2:end)-criterion(:,1:end-1)) == 2;
%         nabla_criterion(:,model.grid.MX(2)-2,2) = false(model.grid.MX(1)-2,1);
%         nabla_criterion(:,2:end,2) = nabla_criterion(:,2:end,2) ...
%             | nabla_criterion(:,1:end-1,2);
%         %         nabla_criterion(:,:,1,2) = abs ...
%         %             (criterion(2:end,:)-criterion(1:end-1,:)) == 2;
%         %         nabla_criterion(:,:,2) = abs ...
%         %             (criterion(:,2:end)-criterion(:,1:end-1));% Mx-3 My-3 2
%         %         criterion_mixing1 = nabla_criterion == 2;
%         criterion_mixing = nabla_criterion(:,:,1) | nabla_criterion(:,:,2);
%
%         %         N_ech=300;
%         % %         a_estim = 2 * criterion_mixing ...
%         % %              * prod(model.grid.dX)/time_t;
%         %         a_estim =  2 * alpha2 ./ (1 +alpha2)  ...
%         %              * prod(model.grid.dX)/time_t;
%         % %         a_estim =  2 * alpha2 .* criterion_mixing ...
%         % %              * prod(model.grid.dX)/time_t;
%         % %          a_estim(a_estim>0)
%         %         sigma_B_t = zeros([model.grid.MX 2 N_ech]);
%         %         sigma_B_t(2:end-1,2:end-1,:,:) = bsxfun(@times, sqrt(a_estim) * sqrt(time_t), ...
%         %             randn([model.grid.MX-2 2 N_ech ]));
%         %          figure(8);plot(sigma_B_t(sigma_B_t>0)/sqrt(prod(model.grid.dX)))
%         %         X_rand = bsxfun(@plus, X, reshape(sigma_B_t,[prod(model.grid.MX) 2 N_ech]));
%         %         X_rand(:,1,:)=mod(X_rand(:,1,:),model.grid.MX(1)*model.grid.dX(1));
%         %         X_rand(:,2,:)=mod(X_rand(:,2,:),model.grid.MX(2)*model.grid.dX(2));
%         % %         keyboard;
%         %         for k=1:N_ech
%         %             T_adv_rand(:,k) = interp2(X0{1},X0{2},T0',X_rand(:,1,k),X_rand(:,2,k));
%         %         end
%         %         T_adv = 1/N_ech * sum(T_adv_rand,2);
%         %         T_adv=reshape(T_adv,MX);
%         %         fft_T_adv=fft2(T_adv);
%         %
%         %         figure(6)
%         %         subplot(2,1,1)
%         %         fct_spectrum(model, fft_T0);
%         %         hold on;
%         %         fct_spectrum(model,fft_T_adv);
%         %         hold off
%         %         ax = axis;
%         %         ax(3)=0;
%         %         axis(ax);
%         %         subplot(2,1,2)
%         %         imagesc(x,y,10*T_adv(:,1:My)'+20)
%
%         %         figure(9)
%         subplot(2,3,2)
%         imagesc(x,y,criterion_mixing(:,1:My-3)')
%         axis equal
%         %         imagesc(x,y,log(a_estim(:,1:My-3)'))
%         title('mixing areas');

%                 keyboard;
drawnow
%
%
%         eval( ['print -depsc ' pwd '/plots/lagrangian_mixing_t='...
%             num2str(t) '.eps']);
%
%         drawnow


eval( ['print -depsc ' model.folder.folder_simu '/mezic_mixing/mezic_mixing_t='...
    num2str(t) '.eps']);
%         keyboard;

alpha2_m = mean(alpha2(:));
