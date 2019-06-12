function fct_plot(model,fft_b_adv_part,day)
% This function creates some plot online and save it
%

%% Get paramters

% Grid
x = model.grid.x;
y = model.grid.y;
My = model.grid.MX(2);

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
folder_simu = model.folder.folder_simu;
plot_moments = model.advection.plot_moments;
map = model.folder.colormap;

%% One particle
X0=[0 0];
T_adv_part = real(ifft2( fft_b_adv_part(:,:,1,id_part) ));

if strcmp(model.type_data,'Back_grad') ...
        && isfield(model, 'add_back_grad')
    [X,Y]=ndgrid(x,y);
    T_adv_part = T_adv_part + init_Back_grad(model,X,Y);
end

if ( (eval(day) == 0) && ...
        strcmp(model.type_data,'Perturbed_vortices') )
    width = 3.2;
    height = 3.2;
    figure1=figure(1);
    set(figure1,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    contourf(x,y,T_adv_part');
    x= model.grid.dX(1)*(0:model.grid.MX(1)-1);
    y= model.grid.dX(2)*(0:model.grid.MX(2)-1);
    Lx = model.grid.dX(1) * model.grid.MX(1);
    sigma= 2 * Lx/15;
    center1x=x(1/4*model.grid.MX(1)+1);
    center1y=y(1/4*model.grid.MX(2)+1);
    nsig=40;
    dist = 1.5;
    rate = 0.3;
    sigma = sigma/nsig;
    center1x= 2e4;
    center1y=y(1/4*model.grid.MX(2)+1);
    coord1=[center1x center1y];
    size_square = 10e4;
    redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
    hold on;
    if strcmp(model.type_data,'Perturbed_vortices')
        plot(redline1(:,1),redline1(:,2),'r','LineWidth',3);
    elseif strcmp(model.type_data,'spot6')
        plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
    else
        error('wrong type of data?');
    end
    center1x= x(1/2*model.grid.MX(1)+1) - 2e4;
    center1y=y(3/4*model.grid.MX(2)+1);
    coord1=[center1x center1y];
    redline1 = [ [coord1(1)-size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)-size_square/2] ; ...
        [coord1(1)+size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)+size_square/2] ; ...
        [coord1(1)-size_square/2 coord1(2)-size_square/2] ];
    plot(redline1(:,1),redline1(:,2),'b','LineWidth',3);
    hold off;
    
else
    width = 3.3;
    height = 3.2;
    figure1=figure(1);
    set(figure1,'Units','inches', ...
        'Position',[X0(1) X0(2) width height], ...
        'PaperPositionMode','auto');
    imagesc(x,y,T_adv_part');
end
caxis([-1 1]*1e-3);
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
title({'One realization', ...
    ['\hspace{0.5cm} $t=' num2str(day) '$ day ']},...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
axis xy; axis equal
colormap(map)
colorbar
drawnow
eval( ['print -depsc ' folder_simu '/one_realization/'...
    num2str(day) '.eps']);       

%% Spectrum
X0=[3.3 1];
close(figure(4))
figure4=figure(4);

widthtemp = 12;
heighttemp = 6;
set(figure4,'Units','inches', ...
    'Position',[X0(1) X0(2) widthtemp heighttemp], ...
    'PaperPositionMode','auto');
fct_spectrum( model,fft_b_adv_part(:,:,:,id_part),'b');
set(gca,'XGrid','on','XTickMode','manual');
width = 4;
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
ylabel('$|\hat{b}(\kappa)|^2$',...
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
title({'Spectrum of' ...
      '\hspace{0.5cm} one realization'},...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
drawnow
eval( ['print -depsc ' folder_simu '/Spectrum/' day '.eps']);

%% Moments
if plot_moments
    T_adv_part = real(ifft2(fft_b_adv_part));
    tol = 1e-2;
    if strcmp(type_data,'Spectrum')
        tol = 1e-1;
    end
    mean_T = mean(T_adv_part,4);
    std_T = std(T_adv_part,0,4);
    odg_b = sqrt( mean(mean_T(:).^2) );
    
    % First and second order moments
    X0=[0 4.2];
    width = 3.65;
    height = 3;
    figure2=figure(2);
    set(figure2,'Units','inches', ...
        'Position',[X0(1) X0(2) 2*width height], ...
        'PaperPositionMode','auto');
    subplot(1,2,1)
    subimage(x,y,mean_T');axis xy;
    imagesc(x,y,mean_T');axis xy;
    axis equal
    caxis([-1 1]*1e-3);
    if model.folder.colormap_freeze
        colormap(map);
        colorbar;
        cbfreeze;
    else
        colormap(map);
        ax1 = gca;
        colorbar('peer',ax1);
        colormap(ax1,map);
    end
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
    title('\hspace{0.5cm} Mean',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    subplot(1,2,2)
    subimage(x,y,std_T');
    imagesc(x,y,std_T');axis xy;
    axis equal
    caxis([0 1.5e-4]);
    if strcmp(type_data,'Spectrum')
        caxis([0 1e-3]);
    end
    if model.folder.colormap_freeze
        colormap('default');
        colorbar;
    else
        ax2 = gca;
        colorbar('peer',ax2);
        colormap(ax2,'default');
    end
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
    title('\hspace{1cm} Standard deviation',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    drawnow;
    eval( ['print -depsc ' folder_simu '/1st_2nd_order_moments/' day '.eps']);
    
    % Third and fourth order moments
    X0=[0 8.2];
    figure3=figure(3);
    set(figure3,'Units','inches', ...
        'Position',[X0(1) X0(2) 2*width height], ...
        'PaperPositionMode','auto');
    % These moments are shown only where the variance is large enough
    jjj = ( std_T < tol *odg_b );
    s=size(std_T);
    std_T=std_T(:);
    % The standard deviation used to compute skewness and kurtosis is set 
    % to infinity when this standard deviation is too low
    % As such, skewness and kurtosis are equal to zero when this standard 
    % deviation is too low
    std_T(jjj(:))=inf;
    std_T=reshape(std_T,s);
    
    % Centering
    T_prime = bsxfun(@plus,  real(ifft2( fft_b_adv_part)) , - mean_T);
    % Normalization
    T_prime = bsxfun(@times, T_prime, 1./ std_T);
    m3 = T_prime.^3 ;
    m3 = mean(m3,4);
    m4 = T_prime.^4 ;
    m4 = mean(m4,4);
    
    sm4=size(m4);
    m4=m4(:);
    m4(m4<3)=3;
    m4=reshape(m4,sm4);
    
    subplot(1,2,1)
    subimage(x,y,m3');
    imagesc(x,y,m3');axis xy;
    axis equal
    caxis(2*[-1 1]);
    if strcmp(type_data,'Spectrum')
        caxis(0.5*[-1 1]);
    end
    if model.folder.colormap_freeze
        colormap(map);
        colorbar;
        cbfreeze;
    else
        colormap(map);
        ax1 = gca;
        colorbar('peer',ax1);
        colormap(ax1,map);
    end
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
    title('\hspace{0.5cm} Skewness',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    
    subplot(1,2,2)
    subimage(x,y,log(m4'-3));
    imagesc(x,y,log(m4'-3));axis xy;
    axis equal
    if model.folder.colormap_freeze
        colormap('default');
        colorbar;
    else
        ax2 = gca;
        colorbar('peer',ax2);
        colormap(ax2,'default');
    end
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
    title('\hspace{1.5cm} log(Kurtosis-3)',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'interpreter','latex',...
        'FontSize',12,...
        'FontName','Times')
    if ~ strcmp(type_data, 'Spectrum')
        caxis([-5 3]);
    end
    drawnow
    eval( ['print -depsc ' folder_simu '/3rd_4th_order_moments/' day '.eps']);
end


end

