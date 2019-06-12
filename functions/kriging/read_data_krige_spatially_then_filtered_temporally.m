init

%%

vectname_data = {'zeta'};
type = 'filtered';
% n_sub = 4*20;
n_sub = 4;

% Colormap
load('BuYlRd.mat');
map = BuYlRd; clear BuYlRd

x= 900e3/1200 * (0:(1200-1));
y= 1050e3/1400 * (0:(1400-1));

xx= x(500:end);
yy= y(1:1100);
model.grid.x = xx;
model.grid.y = yy;
model.grid.origin = [xx(1) yy(2)];
model.grid.dX = [xx(2)-xx(1) yy(2)-yy(1)];
model.grid.MX = [length(xx) length(yy)];
model.grid.BOX = [ xx(1) xx(end)  ; yy(1) yy(end) ];

[xx,yy]=ndgrid(xx,yy);
% xt=stk_dataframe([xx(:) yy(:)]);

xx_LR = x(500:n_sub:end);
yy_LR = y(1:n_sub:1100);
% yy_LR = yy_LR(1:end-50);
xx_LRref = xx_LR;
yy_LRref = yy_LR;
model.grid_LR.x = xx_LR;
model.grid_LR.y = yy_LR;
model.grid_LR.origin = [xx_LR(1) yy_LR(2)];
model.grid_LR.dX = [xx_LR(2)-xx_LR(1) yy_LR(2)-yy_LR(1)];
model.grid_LR.MX = [length(xx_LR) length(yy_LR)];
model.grid_LR.BOX = [ xx_LR(1) xx_LR(end)  ; yy_LR(1) yy_LR(end) ];

[xx_LR,yy_LR]=ndgrid(xx_LR,yy_LR);
xt=stk_dataframe([xx_LR(:) yy_LR(:)]);

%% Kriging

% Algorithnm parametrisation
model.error_obs_h=true;
model.kriging.param_estimated_by_MLE=false;
model.grid_trace_larger = false;

grid_trace = fct_grid_satellite(model);
grid_trace = grid_trace';
% figure;
% plot(grid_trace(:,1),grid_trace(:,2),'g.')

if model.error_obs_h
    model.measurement_error_h=1e-2;
else
    model.measurement_error_h=0;
end
% [model.kriging,model.obs.x,model.obs.h]=kriging_on_high(error_obs,param_estim);
if ~model.kriging.param_estimated_by_MLE
    model.kriging = krige_spatially2(model);
end

%     N_t = 1;
%     N_t = 738-380+2;
    t_ini = 380;
%     t_final = 380;
    t_final = 738;
    N_t = t_final-t_ini+2;

for k=1:length(vectname_data)
    
%     x= 900e3/1200 * (0:(1200-1));
%     y= 1050e3/1400 * (0:(1400-1));

    name_data=vectname_data{k};
    fprintf([name_data ' \n'])
    data = nan(1200,1400,1,N_t);
%     data = nan(1201,1401,1,N_t);
    for t2 = t_ini:2:t_final
        data_temp=ncread( [ ...
        '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
        '/CHABU_surface/chabu_surf.0' num2str(t2) '.nc'],name_data );
        data(:,:,1,t2-t_ini+1:(t2-t_ini+2))=data_temp(2:1201,2:1401,:,:);
%         data(:,:,1,t2-t_ini+1:(t2-t_ini+2))=data_temp(1:1201,1:1401,:,:);
    end
%     data=data(:,:,:,1);
    time=ncread([ ...
        '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
        '/CHABU_surface/chabu_surf.0380.nc'],'ocean_time');
    dt= time(2)-time(1);
    fprintf('loaded \n')
    MX = size(data);
    MX = MX(1:2);
    data = reshape(data,[ MX 1 N_t]);
    
    %% Remove cost
    sum(isnan(data(:)))/length(isnan(data(:)))
    data=data(500:end,1:1100,:,:);
%     data(:,640+1:end,:,:)=[];
    sum(isnan(data(:)))/length(isnan(data(:)))
    
    MX = size(data);
    x= x(500:end);
    y= y(1:1100);
    dX = [ x(2)-x(1) y(2)-y(1) ];
    [X,Y]=meshgrid(x,y);
%     [X,Y]=ndgrid(x,y);
    
    %% Filter
    data_temp = nan([model.grid_LR.MX 1 size(data,4)]);
    for t3 = 1:size(data,4)
        
        %         [min(grid_trace(:,1)) max(grid_trace(:,1))]
        %         [min(x) max(x)]
        %         [min(grid_trace(:,2)) max(grid_trace(:,2))]
        %         [min(y) max(y)]
        width=8;
        height=5;
        X0=[0 0];
        taille_police = 11;
        close;
        figure4=figure(10);
        set(figure4,'Units','inches', ...
            'Position',[X0(1) X0(2) width height], ...
            'PaperPositionMode','auto');
        subplot(1,2,1)
        imagesc(x,y,data(:,:,:,t3)');axis equal;axis xy;
        hold on;
        plot(grid_trace(:,1),grid_trace(:,2),'g.')
        colormap(map);colorbar;
        set(gca,...
            'Units','normalized',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
        ylabel('$y(m)$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('$x(m)$',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'interpreter','latex',...
            'FontName','Times')
        title('\hspace{1cm} High resolution',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'interpreter','latex',...
            'FontSize',12,...
            'FontName','Times')

%         data_save = data(:,:,:,t3);
        
        % Satellite measurement
        ssh_on_trace = ...
            interp2(X,Y,data(:,:,:,t3)', grid_trace(:,1), grid_trace(:,2),'spline');
        
        % Measurement error
        ssh_on_trace  = ssh_on_trace ...
            + model.measurement_error_h * sqrt(2) * randn(size(ssh_on_trace)) ;
        
        if model.kriging.param_estimated_by_MLE
            model.kriging,BOX ...
                = krige_spatially_ML(model,ssh_on_trace,grid_trace);
        end
        
        % Interpolation by kriging
        % Here, we compute the kriging prediction on each point of the grid
        xi = stk_dataframe (grid_trace);
        zi = stk_dataframe(ssh_on_trace);
        zp = stk_predict(model.kriging, xi, zi, xt);
        data_temp(:,:,:,t3) = reshape(double(zp.mean),model.grid_LR.MX);
        
        subplot(1,2,2)
        imagesc(xx_LRref,yy_LRref,data_temp(:,:,:,t3)');axis equal;axis xy;
        colormap(map);colorbar;
        set(gca,...
            'Units','normalized',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
        ylabel('$y(m)$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('$x(m)$',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'interpreter','latex',...
            'FontName','Times')
        title('\hspace{1cm} Interpolation',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'interpreter','latex',...
            'FontSize',12,...
            'FontName','Times')
        drawnow
        %         eval( kriging_plots
        %         if plot_high_kriging
        %             fct_plot_ssh_kriging(model.kriging,x,y,data(:,:,:,t3),grid_trace)
        %         end
        eval( ['print -depsc ' pwd '/images/kriging_plots/' num2str(t3) '.eps'])
        
        fprintf(['Kriged spatially for day' num2str(t3) '\n'])
    end
    data = data_temp;
    fprintf(['Kriged spatially for all days \n'])
    
    %% Remove cost
%     sum(isnan(data(:)))/length(isnan(data(:)))
%     data(:,512/n_sub+2+1:end,:,:)=[];
%     sum(isnan(data(:)))/length(isnan(data(:)))
    %%
    n_day=1; % 7 2
    [data,data_var,n_day,residu_data] = fct_filtered_several_day(data,n_day,dt);
    fprintf('filtered temporally\n')
    
    %% Save
    dX = dX * n_sub;
    save([ ...
        'images/krige_space_then_time/' name_data ...
        '_filtered.mat'],'data','model');
%     save([ ...
%         '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
%         '/images/krige_space_then_time/' name_data ...
%         '_filtered.mat'],'data','model');
    clear data
    fprintf('saved \n')
    
    
end
