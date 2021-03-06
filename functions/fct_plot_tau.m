function [on_tau2_global, on_tau_local] = fct_plot_tau(model,...
    on_tau_local_folding,on_tau_local_stretching,rate_switch,cax_alpha)
% function on_tau_local = fct_on_tau_2(model,nabla_phi,w_eul,w_Lag)
% Compute the inverse of the local time scale of differntial advection
%

% figure(25);imagesc(rate_switch);axis xy; axis equal;

% rate_switch = rate_switch >= pi ; 
% pas mal non plus et defendable un peu theoriquement mais tau un peu surestime
% rate_switch = rate_switch >= 4 ; % plus joli mais tau un peu surestime
rate_switch = rate_switch >=2 ; % pas mal

% figure(26);imagesc(rate_switch);axis xy; axis equal;

s=size(on_tau_local_folding);
on_tau_local = on_tau_local_folding(:);
on_tau_local(~rate_switch) = on_tau_local_stretching(~rate_switch);
on_tau_local =reshape(on_tau_local,s);
% on_tau_local = rate_switch .* on_tau_local_folding ...
%     + (1-rate_switch) .* on_tau_local_stretching;

% on_tau_local = min( on_tau_local_folding,on_tau_local_stretching);
on_tau2_global = 1/prod(model.grid.MX) * sum(on_tau_local(:).^2);
tau_global = 1/sqrt(on_tau2_global) /(3600*24)

% figure(24)
% imagesc(log(on_tau_local'.^2));axis xy;axis equal

%% Plot log

% Grid
x = model.grid.x_ref;
y = model.grid.y_ref;

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
map = 'default';
loc_colorbar = 'southoutside';
% map = model.folder.colormap;

width=9;
% width=12;
height=4;

% width = 13;
% % width = 3.3;
% %     height = 3.2;
ax = [x(end)-x(1) y(end)-y(1)] ;
aspect_ratio = ax(2)/ax(1);
% height = 0.265*aspect_ratio * width;
% % height = 0.26*aspect_ratio * width;
X0=[0 10];

figure(22);
close;
figure22=figure(22);
set(figure22,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');



subplot(1,3,3)
imagesc(x,y,log((on_tau_local.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $-2ln(\tau)$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)
drawnow
cax=caxis;
% cax(2)=max( [ max(on_tau_local_folding(:)) ...
%     max(on_tau_local_stretching(:)) ...
%     max(on_tau_local(:)) ]);
% caxis(cax);

% imagesc(x,y,log((on_tau_local.^2)'));
subplot(1,3,1)
imagesc(x,y,log((on_tau_local_folding.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $-2ln(\tau_f)$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
% colorbar
colorbar('location',loc_colorbar)
drawnow
caxis(cax);

subplot(1,3,2)
imagesc(x,y,log((on_tau_local_stretching.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $-2ln(\tau_s)$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)
drawnow
caxis(cax);

eval( ['print -depsc ' folder_simu '/log_on_tau_2_all.eps']);
% keyboard;


%% Plot

figure(23);
close;
figure23=figure(23);
set(figure23,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');


subplot(1,3,3)
imagesc(x,y,((on_tau_local.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $1/\tau^2$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)
drawnow
cax=caxis;
if nargin == 5
   cax(2)=cax_alpha;
   caxis(cax);
end

% imagesc(x,y,log((on_tau_local.^2)'));
subplot(1,3,1)
imagesc(x,y,((on_tau_local_folding.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $1/\tau_f^2$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)
drawnow
caxis(cax);

subplot(1,3,2)
imagesc(x,y,((on_tau_local_stretching.^2)'));
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('x(m)',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
title('\hspace{0.5cm} $1/\tau_s^2$ ',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
axis([x(1) x(end) y(1) y(end)])
colormap(map)
%colorbar
colorbar('location',loc_colorbar)
drawnow
caxis(cax);

eval( ['print -depsc ' folder_simu '/on_tau_2_all.eps']);
% keyboard;

end

function f_ortho = fct_ortho(f)
% Compute the orthogonal vector in each point of the space
%
f_ortho(:,:,1)= - f(:,:,2);
f_ortho(:,:,2)= + f(:,:,1);
end

function g = fct_normalize(g)
% Compute the orthogonal vector in each point of the space
%
ng = sqrt(sum(g.^2,3));
g = bsxfun( @times, 1./ng , g);
end