init

%%
% model.type_data = 'Loop_current_kriging';
model.type_data = 'Gula_kriging';
param_estimated_by_MLE  = true;


model.grid.dX = 500e3/1024 * [1 1];
model.grid.MX = [1538 512];

n_x_min = 501;
n_y_max = 900;

error_obs_h = true;
init_model

vectname_data = {'zeta'};
type = 'filtered';
% n_sub = 4*20;
n_sub = 4;


switch model.type_data
    case 'Loop_current_kriging'
        model.folder.data = ...
            ['/Volumes/WD_Ressegui/These/from_tour/data_brest/' ....
            'data_molemaker/data/' ];
    case 'Gula_kriging'
        model.folder.data = ...
            ['/Volumes/WD_Ressegui/These/from_tour/data_brest/'...
            'data_gula/CHABU_surface/'];
end

% Colormap
load('BuYlRd.mat');
map = BuYlRd; clear BuYlRd


% dt = 24 * 3600;
% model.grid.x = model.grid.dX(1) * ( 0:(model.grid.MX(1)-1));
% model.grid.y = model.grid.dX(2) * ( 0:(model.grid.MX(2)-1));
% x=model.grid.x;
% y=model.grid.y;
x= 900e3/1200 * (0:(1200-1));
y= 1050e3/1400 * (0:(1400-1));

% xx= x(500:end);
xx= x(n_x_min:end);
yy= y(1:n_y_max);
model.grid.x = xx;
model.grid.y = yy;
model.grid.origin = [xx(1) yy(2)];
model.grid.dX = [xx(2)-xx(1) yy(2)-yy(1)];
model.grid.MX = [length(xx) length(yy)];
model.grid.BOX = [ xx(1) xx(end)  ; yy(1) yy(end) ];

[xx,yy]=ndgrid(xx,yy);
% xt=stk_dataframe([xx(:) yy(:)]);
xt_HR=stk_dataframe([xx(:) yy(:)]);

xx_LR = x(n_x_min:n_sub:end);
yy_LR = y(1:n_sub:n_y_max);
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
model.error_obs_h=error_obs_h;
% model.kriging.param_estimated_by_MLE=true;
model.kriging.param_estimated_by_MLE=param_estimated_by_MLE;
model.grid_trace_larger = false;

% Position of measurement : satellite trace
grid_trace = fct_grid_satellite(model);
grid_trace = grid_trace';
xi = stk_dataframe (grid_trace);

% grid_trace_index = fct_grid_satellite_index(model);


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
% t_final = 380;
t_final = 61*2+380;
%     t_final = 738;
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
    %     sum(isnan(data(:)))/length(isnan(data(:)))
    data=data(n_x_min:end,1:n_y_max,:,:);
    % %     data(:,640+1:end,:,:)=[];
    %     sum(isnan(data(:)))/length(isnan(data(:)))
    
    MX = size(data);
    x= x(n_x_min:end);
    y= y(1:n_y_max);
    dX = [ x(2)-x(1) y(2)-y(1) ];
    [X,Y]=meshgrid(x,y);
    %     [X,Y]=ndgrid(x,y);
    
    
    % fct_plot_f(model,data(:,:,:,1));
    
    %% Filter
    data_temp = nan([model.grid_LR.MX 1 size(data,4)]);
    data_temp_HR = nan([model.grid.MX 1 size(data,4)]);
%     for t3 = 1


        model_ref = model;
        
      for t3 = 1:size(data,4)
        t3
        model = model_ref;
        
        %% Satellite measurement
        %         ssh_on_trace = ...
        %             interp2(X,Y,data(:,:,:,t3)', grid_trace(:,1), grid_trace(:,2),'spline');
        %         zi = stk_dataframe(ssh_on_trace);
        
        %% Satellite measurement trace by trace
        
        grid_trace_index = fct_grid_satellite_index(model);
        
        
        ssh_on_trace = [];
        grid_trace =[];
        nb_trace = length(grid_trace_index);
        
        
        
        for k=1:nb_trace
            % Trace k
            trace_k = grid_trace_index{k};
            
            % Sort measurement points of the trace k
            [~,I] = sort(trace_k(3,:));
            trace_k=trace_k(:,I);
            
            % Measurements of the trace k
            ssh_on_trace_k = ...
                interp2(X,Y,data(:,:,:,t3)', ...
                trace_k(1,:), trace_k(2,:),'spline');
            
            % Measurement with error
            ssh_on_trace_k  = ssh_on_trace_k ...
                + model.measurement_error_h ...
                * sqrt(2) * randn(size(ssh_on_trace_k)) ;
            trace_k(4,:) = ssh_on_trace_k;
            
            % Gather all traces
            ssh_on_trace = [ ssh_on_trace; ssh_on_trace_k'];
            grid_trace = [ grid_trace; ...
                [trace_k(1,:)' trace_k(2,:)']];
            
            grid_trace_index{k} = trace_k;
        end
        fct_plot_f_under_trace(model,grid_trace_index);
        
        %% Kriging
        xi = stk_dataframe(grid_trace);
        zi = stk_dataframe(ssh_on_trace);
        
        if model.kriging.param_estimated_by_MLE
            
            %             kriging_model.lognoisevariance = ...
            %                 2*log(model.measurement_error_h * sqrt(2));
            %             kriging_model.lognoisevariance=log(nugget_val);
            
            model.kriging.model_order= 0; 
            % ordinnary krigin : unknown  constant mean
            model_ = model;
            nugget_val=(5.6e-2)^2; 
            % from GILLE & KELLY 1996 (Southern ocean)
            % (ecart type 10 fois+ grand que estime MLE : 0.63cm)
            model_.measurement_error_h = sqrt( nugget_val/2 );
            model.kriging ...
                = krige_spatially_ML(model_,ssh_on_trace,grid_trace);
            
%             model.kriging ...
%                 = krige_spatially_ML(model,ssh_on_trace,grid_trace);
        end
        
        
        %% Satellite residual measurement trace by trace
        grid_trace_index_residual = grid_trace_index;
        residual_ssh_on_trace = [];
        for k=1:nb_trace
            % Trace k
            trace_k = grid_trace_index_residual{k};
            
            % Mesurement
            ssh_on_trace_k = trace_k(4,:);
            
            % Kriging esimate
            dataframe_k = stk_dataframe ([trace_k(1,:)' trace_k(2,:)']);
            kriged_ssh_on_trace_k = ...
                stk_predict(model.kriging, xi, zi, dataframe_k);
            
            % Residuals
            trace_k(4,:) = ssh_on_trace_k ...
                - double(kriged_ssh_on_trace_k.mean)';
            
            % Test
            %             warning('test')
            %             trace_k(4,:) = ones(size(trace_k(4,:)));
            
            % Gather all traces
            residual_ssh_on_trace = [ residual_ssh_on_trace; trace_k(4,:)'];
            
            grid_trace_index_residual{k} = trace_k;
        end
        zi_residuals = stk_dataframe(residual_ssh_on_trace);
        
        fct_plot_f_under_trace(model,grid_trace_index_residual);
        
        %% Measurement error
        %         ssh_on_trace  = ssh_on_trace ...
        %             + model.measurement_error_h * sqrt(2) * randn(size(ssh_on_trace)) ;
        %
        %         if model.kriging.param_estimated_by_MLE
        %             model.kriging ...
        %                 = krige_spatially_ML(model,ssh_on_trace,grid_trace);
        %         end
        
        %% Interpolation by kriging and residuals fields
        %         zp = stk_predict(model.kriging, xi, zi, xt);
        %         data_temp(:,:,:,t3) = reshape(double(zp.mean),model.grid_LR.MX);
        
        
%         if exist([model.folder.folder_simu ...
%                 '/Kriged_SSH_on_HR.mat'] ,'file') == 2 && ...
%                 ~ model.kriging.param_estimated_by_MLE
%             load([model.folder.folder_simu '/Kriged_SSH_on_HR.mat'],'zp');
%         elseif exist([model.folder.folder_simu ...
%                 '/Kriged_SSH_on_HR_estim_MLE.mat'] ,'file') == 2 && ...
%                 model.kriging.param_estimated_by_MLE
%             load([model.folder.folder_simu '/Kriged_SSH_on_HR_estim_MLE.mat'],'zp');
%         else
%             zp = stk_predict(model.kriging, xi, zi, xt_HR);
%             if ~ model.kriging.param_estimated_by_MLE
%                 save([model.folder.folder_simu '/Kriged_SSH_on_HR'],...
%                     'model','zp');
%             else
%                 save([model.folder.folder_simu '/Kriged_SSH_on_HR_estim_MLE'],...
%                     'model','zp');
%             end
%         end
        
        
        zp = stk_predict(model.kriging, xi, zi, xt_HR);
        if ~ model.kriging.param_estimated_by_MLE
            save([model.folder.folder_simu '/Kriged_SSH_on_HR_time_' ...
                num2str(t3)],...
                'model','zp');
        else
            save([model.folder.folder_simu '/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(t3)],...
                'model','zp');
        end
        
        data_temp_HR(:,:,:,t3) = reshape(double(zp.mean),model.grid.MX);
        data_residuals(:,:,:,t3) = data_temp_HR(:,:,:,t3) - data(:,:,:,t3);
        
        %% Plots
        fct_plot_f(model,data(:,:,:,t3));
        fct_plot_f(model,data_temp_HR(:,:,:,t3));
        fct_plot_f(model,data_residuals(:,:,:,t3));
        %         fct_plot_f(model,abs(data_residuals(:,:,:,t3)));
        close all
        
        %% Velocity
        
        % Physical parameters
        model = fct_physical_param(model);
        switch model.type_data
            case 'Loop_current_kriging'
                model.physical_constant.f = ...
                    ncread(...
                    ['/Volumes/WD_Ressegui/These/from_tour/data_brest/' ....
                    'data_molemaker' ...
                    'data/GOM_500_grd.nc'],'f');
            case 'Gula_kriging'
                model.physical_constant.f = ...
                    ncread(...
                    ['/Volumes/WD_Ressegui/These/from_tour/data_brest/'...
                    'data_gula' ...
                    '/CHABU_surface/chabu_surf.0'],'f');                
        end
        
        
        % Coriolis
        nabla_ssh = gradient_mat_2( data(:,:,:,t3),model.grid.dX);
        w_smooth(:,:,1,:)= - nabla_ssh(:,:,2,:);
        w_smooth(:,:,2,:)= nabla_ssh(:,:,1,:);
        clear nabla_ssh;
%         w_smooth = model.physical_constant.g / model.physical_constant.f0 ...
        w_smooth = bsxfun(@times, ...
            model.physical_constant.g ./ model.physical_constant.f ...
                    , w_smooth);
        
        % Coriolis
        nabla_ssh = gradient_mat_2( data_temp_HR(:,:,:,t3),model.grid.dX);
        w_HR(:,:,1,:)= - nabla_ssh(:,:,2,:);
        w_HR(:,:,2,:)= nabla_ssh(:,:,1,:);
        clear nabla_ssh;
        w_HR = bsxfun(@times, ...
            model.physical_constant.g ./ model.physical_constant.f ...
                    , w_HR);
%         w_HR = model.physical_constant.g / model.physical_constant.f0 * w_HR;

        
        if ~ model.kriging.param_estimated_by_MLE
            save([model.folder.folder_simu '/multi_time/w_Kriged_SSH_on_HR_time_' ...
                num2str(t3)],...
                'model','w_HR','w_smooth');
        else
            save([model.folder.folder_simu '/multi_time/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(t3)],...
                'model','w_HR','w_smooth');
        end
        
%         keyboard;
        
        %% Plot true SSH field
        %         %         [min(grid_trace(:,1)) max(grid_trace(:,1))]
        %         %         [min(x) max(x)]
        %         %         [min(grid_trace(:,2)) max(grid_trace(:,2))]
        %         %         [min(y) max(y)]
        %         width=12;
        %         height=5;
        %         X0=[0 0];
        %         taille_police = 11;
        %         close;
        %         figure4=figure(10);
        %         set(figure4,'Units','inches', ...
        %             'Position',[X0(1) X0(2) width height], ...
        %             'PaperPositionMode','auto');
        %         subplot(1,3,1)
        %         imagesc(x,y,data(:,:,:,t3)');axis equal;axis xy;
        %         hold on;
        %         plot(grid_trace(:,1),grid_trace(:,2),'g.')
        %         colormap(map);colorbar;
        %         set(gca,...
        %             'Units','normalized',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         ylabel('$y(m)$',...
        %             'FontUnits','points',...
        %             'interpreter','latex',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         xlabel('$x(m)$',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'interpreter','latex',...
        %             'FontName','Times')
        %         title('\hspace{1cm} High resolution',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'interpreter','latex',...
        %             'FontSize',12,...
        %             'FontName','Times')
        %
        % %         data_save = data(:,:,:,t3);
        %
        %         %% Plot kriged ssh
        %
        %         subplot(1,3,2)
        %         imagesc(x,y,data_temp_HR(:,:,:,t3)');axis equal;axis xy;
        %         colormap(map);colorbar;
        %         set(gca,...
        %             'Units','normalized',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         ylabel('$y(m)$',...
        %             'FontUnits','points',...
        %             'interpreter','latex',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         xlabel('$x(m)$',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'interpreter','latex',...
        %             'FontName','Times')
        %         title('\hspace{1cm} Interpolation',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'interpreter','latex',...
        %             'FontSize',12,...
        %             'FontName','Times')
        %         drawnow
        %         %         eval( kriging_plots
        %         %         if plot_high_kriging
        %         %             fct_plot_ssh_kriging(model.kriging,x,y,data(:,:,:,t3),grid_trace)
        %         %         end
        % %         eval( ['print -depsc ' pwd '/images/kriging_plots/' num2str(t3) '.eps'])
        %
        %         %% Plot residual ssh
        %
        %         subplot(1,3,3)
        %         imagesc(x,y,data_residuals(:,:,:,t3)');axis equal;axis xy;
        %         colormap(map);colorbar;
        %         set(gca,...
        %             'Units','normalized',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         ylabel('$y(m)$',...
        %             'FontUnits','points',...
        %             'interpreter','latex',...
        %             'FontSize',taille_police,...
        %             'FontName','Times')
        %         xlabel('$x(m)$',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'FontSize',taille_police,...
        %             'interpreter','latex',...
        %             'FontName','Times')
        %         title('\hspace{1cm} Residual',...
        %             'FontUnits','points',...
        %             'FontWeight','normal',...
        %             'interpreter','latex',...
        %             'FontSize',12,...
        %             'FontName','Times')
        %         drawnow
        %         %         eval( kriging_plots
        %         %         if plot_high_kriging
        %         %             fct_plot_ssh_kriging(model.kriging,x,y,data(:,:,:,t3),grid_trace)
        %         %         end
        % %         eval( ['print -depsc ' pwd '/images/kriging_plots/' num2str(t3) '.eps'])
        %
        %         fprintf(['Kriged spatially for day ' num2str(t3) '\n'])
        
        %         keyboard;
        %% MLE on observation residuals
%         %         zp_on_trace = stk_predict(model.kriging, xi, zi, xi);
%         %         zi_residuals = zi - zp_on_trace.mean;
%         %         %         zp_mean_on_trace = reshape(double(zp_on_trace.mean),model.grid_LR.MX);
%         

%         model_temp = model;
%         model_temp = rmfield(model_temp,'kriging');
%         model_temp.kriging.param_estimated_by_MLE = true;
%         model_temp.kriging.model_order= - 1; 
%         % ordinnary krigin : constant mean =0
%         model_temp.kriging.covariance_type = 'matern';
%         model.kriging_residuals ...
%             = krige_spatially_ML(model_temp,zi_residuals,xi);
%         clear model_temp
%         
%         keyboard;
        
      end
    fprintf(['Kriged spatially for all days \n'])
    
end
