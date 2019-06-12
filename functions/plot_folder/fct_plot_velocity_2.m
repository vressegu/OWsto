function fct_plot_velocity(model,w)
% This function creates some plot online and save it
%

%% Get paramters

% Grid
x = model.grid.x;
y = model.grid.y;

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
% folder_simu = model.folder.folder_simu;
plot_moments = false;
map = model.folder.colormap;

%% One particle
X0=[0 0];

% if model.mirror
%     w=w(:,1:model.grid.MX(2)/2,:,:);
%     y=y(1:model.grid.MX(2)/2);
% end

width = 3.3;

%     height = 3.2;

ax = [x(end)-x(1) y(end)-y(1)] ;
aspect_ratio = ax(2)/ax(1);
height = aspect_ratio * width;

figure1=figure(1);
% figure1=figure;
set(figure1,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
imagesc(x,y,sqrt(sum(w.^2,3))');

% caxis([-1 1]*1e-3);
% caxis([0 1]);
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
title(['$\sqrt{KE}$'],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
colormap(map)
colorbar
drawnow
% eval( ['print -depsc ' folder_simu '/one_realization/'...
%     num2str(day) '.eps']);       

%% Spectrum

% Remove boundaries
mask_boundaries = ...
    fct_unity_approx6(model.grid.MX(1))' * ...
    fct_unity_approx6(model.grid.MX(2));
w = bsxfun(@times, mask_boundaries, w);

fft_w = fft2(w);

X0=[3.3 1];
% close(figure(4))
figure4=figure(4);
% figure4=figure;

widthtemp = 12;
heighttemp = 6;
set(figure4,'Units','inches', ...
    'Position',[X0(1) X0(2) widthtemp heighttemp], ...
    'PaperPositionMode','auto');
fct_spectrum2( model,fft_w(:,:,:,id_part),'b');
% fct_spectrum( model,fft_w(:,:,:,id_part),'b');
set(gca,'XGrid','on','XTickMode','manual');
width = 4.5;
% width = 4;
height = 3;
set(figure4,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')

set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('$\overline{\Gamma}_T$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
title('Velocity spectrum',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
drawnow
% eval( ['print -depsc ' folder_simu '/Spectrum/' day '.eps']);



end

