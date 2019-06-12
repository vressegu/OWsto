%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
init;

% type_data =  'Gula';
% type_data =  'Gula_kriging';
% type_data =  'Loop_current';
type_data =  'Loop_current_kriging';
smoothing = true;
subsampling = false;
Mx_LR = 8; % filter of 30 km -> resolution ~ 100km
multitime = true;

model.type_data=type_data;
init_model

%% Load data
switch  type_data
    case  'Back_grad'
        % Day of the simulation
        day = 93;
        day = num2str(day);
        
        % Type of initial condtions
        type_data = 'Back_grad';
        
        % Physical parameters
        model = fct_physical_param();
        
        % Gather parameters in the structure model
        model.type_data=type_data;
        model.folder.folder_simu = [ ...
            '/Volumes/WD_Ressegui/These/from_tour/simu_SQG/images/usual_SQG/' ...
            model.type_data '/files/'];
        
        load([ model.folder.folder_simu  day ]);
        model.add_back_grad = true;
        T_adv = fct_plot_loaded(model,fft_buoy_part,day);
        fft_w = SQG_large_UQ(model, fft_buoy_part);
        w = real(ifft2(fft_w));
    case 'Gula'
        if multitime
            load([model.folder.folder_simu '/w_HR_multi_time.mat'],'w');
            [~,model] = read_data_gula(30,'night',model);
        else
            day = 30;
            %         day = 1;
            time_day = 'night';
            [w,model] = read_data_gula(day,time_day);
            model.advection.meth_anti_alias = 'none';
            model = init_grid_k (model);
        end
    case 'Loop_current'
        
        model.grid.dX = 500e3/1024 * [1 1];
        model.grid.MX = [1536 512];
        % model.grid.MX = [1538 512];
        model.grid.x = model.grid.dX(1) * ( 0:(model.grid.MX(1)-1));
        model.grid.y = model.grid.dX(2) * ( 0:(model.grid.MX(2)-1));
        dt = 0.5 *24 *3600;
        N_t = 180;
        
        w = nan([model.grid.MX 2 N_t]);
        model.folder.data = ...
            ['/Volumes/WD_Ressegui/These/from_tour/data_brest/' ....
            'data_molemaker/data/' ];
        data=ncread( [ ...
            model.folder.data 'surface_summer.nc'],'u' );
        w(:,:,1,:) = data(2:(model.grid.MX(1)+1),2:(model.grid.MX(2)+1),:,:);
        data=ncread( [ ...
            model.folder.data 'surface_summer.nc'],'v' );
        w(:,:,2,:) = data(2:(model.grid.MX(1)+1),2:(model.grid.MX(2)+1),:,:);
%         N_t = size(data,3);
%         t_final = t_ini-1 + N_t;
        clear data
        
    case 'Gula_kriging'
%         % kriged SSH from very high reoslution model output + noise
%         load('images/krige_space_then_time/zeta_filtered.mat',...
%             'data','model');
%         day = 30;
%         ssh=data(:,1:end-1,:,day); clear data;
%         model.grid_HR = model.grid;
%         model.grid = model.grid_LR;
%         model = rmfield(model,'grid_LR');
%         model.grid.y(end) = [];
%         model.grid.MX(2) = model.grid.MX(2)-1;
%         
%         % Physical parameters
%         model = fct_physical_param(model);
%         
%         % Coriolis
%         nabla_ssh = gradient_mat_2( ssh,model.grid.dX);
%         w(:,:,1,:)= - nabla_ssh(:,:,2,:);
%         w(:,:,2,:)= nabla_ssh(:,:,1,:);
%         clear nabla_ssh;
%         w = model.physical_constant.g / model.physical_constant.f0 * w;
%         warning('coriolis not accurate');
        
        smoothing = false;
        
        t_ini = 380;
        % t_final = 380;
        t_final = 61*2+380;
        %     t_final = 738;
        N_t = t_final-t_ini+2;

        load([model.folder.folder_simu '/multi_time/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(1)]);
        s = size(w_HR);
        
        w_HR_tot=nan([s N_t]);
        w_smooth_tot=nan([s N_t]);
        model_tot=model
        
        for t3 = 1:N_t
            load([model.folder.folder_simu '/multi_time/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(t3)]);
            w_HR_tot(:,:,:,t3)=w_HR;
            w_smooth_tot(:,:,:,t3)=w_smooth;
            model_tot(t3)=model;
        end
        w_HR = w_HR_tot;
        w_smooth = w_smooth_tot;
        
        
    case 'Loop_current_kriging'
        smoothing = false;
        
%         t_ini = 380;
%         % t_final = 380;
%         t_final = 61*2+380;
%         %     t_final = 738;
%         N_t = t_final-t_ini+2;
        N_t = 180;

        load([model.folder.folder_simu '/multi_time/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(1)]);
        s = size(w_HR);
        
        w_HR_tot=nan([s N_t]);
        w_smooth_tot=nan([s N_t]);
        model_tot=model
        
        for t3 = 1:N_t
            load([model.folder.folder_simu '/multi_time/Kriged_SSH_on_HR_estim_MLE_time_' ...
                num2str(t3)]);
            w_HR_tot(:,:,:,t3)=w_HR;
            w_smooth_tot(:,:,:,t3)=w_smooth;
            model_tot(t3)=model;
        end
        w_HR = w_HR_tot;
        w_smooth = w_smooth_tot;
end
model.type_data = type_data;
model.filtering.smoothing = smoothing;
model.filtering.subsampling = subsampling;

init_model


% grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% % vort = vorticity_perso(model.grid, w);
% figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
% fct_plot_velocity(model,w);

%% Smoothing
if model.filtering.smoothing
    w_HR = w;
    w_smooth = nan(size(w));
    clear w
    
    name_file = [model.folder.folder_simu '/smooth_w_multitime.mat'];
        
    if strcmp(type_data,'Gula') && ~multitime ...
            && exist([model.folder.folder_simu ...
            '/smooth_w.mat'] ,'file') == 2
         load([model.folder.folder_simu '/smooth_w.mat'],'w');
    elseif strcmp(type_data,'Gula') && multitime ...
            && exist([model.folder.folder_simu ...
            '/smooth_w_multitime.mat'] ,'file') == 2
         load([model.folder.folder_simu '/smooth_w_multitime.mat']);    
    elseif exist(name_file,'file') == 2
            load(name_file);
    else
        Mx_HR = min(model.grid.MX);
        %     Mx_HR = model.grid.MX(1);
        
        if  Mx_LR < Mx_HR
            %% Filter
            sigma_filter = min(model.grid.dX)/2 * Mx_HR/Mx_LR ;
            %     sigma_filter = 1e4; % 10 km
            % %     sigma_filter = 0.36*1/180*pi*model.coriolis.earth_radius; % 0.36 degre
            model.advection.sigma_filter=sigma_filter;
            
            %     fft_filter = design_filter(model);
            %     figure;imagesc(fft_filter);
            %     filter_perso = real(ifft2(fft_filter));
            
            filter_perso = design_filter_3(model);
            %     figure;imagesc(model.grid.x,model.grid.y,filter_perso');axis xy;axis equal
                            %     figure;imagesc(w(:,:,1));
                
            s =size(w_HR);
            N_t = s(4);
            for k_t=1:N_t
                k_t
                w_smooth(:,:,1,k_t) = conv2(w_HR(:,:,1,k_t),filter_perso,'same');
                w_smooth(:,:,2,k_t) = conv2(w_HR(:,:,2,k_t),filter_perso,'same');
            end
            
            
            %     w_f(:,:,1) = filter2(w(:,:,1),filter_perso);
            %     figure;imagesc(w_f(:,:,1));
            
            %         fft_w = fft2(w);
            %         fft_w((Mx_LR/2+1):(end-1-Mx_LR/2+1),:,:,:)=0;
            %         fft_w(:,(Mx_LR/2+1):(end-1-Mx_LR/2+1),:,:)=0;
            %         w = real(ifft2(fft_w));
            % %         fft_buoy_part((Mx_LR/2+1):(end-1-Mx_LR/2+1),:,:,:)=0;
            % %         fft_buoy_part(:,(Mx_LR/2+1):(end-1-Mx_LR/2+1),:,:)=0;
            % %         T_adv = real(ifft2(fft_buoy_part));
            
        else
            warning('the resolution is already too coarse');
        end
        if strcmp(type_data,'Gula') 
            if multitime
                save([model.folder.folder_simu '/smooth_w_multitime.mat'], ...
                    'w_smooth','-v7.3');                
            else
                save([model.folder.folder_simu '/smooth_w.mat'],'w');
            end
        else
           save(name_file , ...
                    'w_smooth','-v7.3');                
        end
        
    end
end

%% Subsampling
% if model.filtering.subsampling
%     Mx_HR = min(model.grid.MX);
%     %     Mx_HR = model.grid.MX(1);
%     if  Mx_LR < Mx_HR
%         model.grid.MX = model.grid.MX * (Mx_LR/Mx_HR);
%         model.grid.dX = model.grid.dX / (Mx_LR/Mx_HR);
%         model.grid.x = model.grid.x(1:Mx_HR/Mx_LR:end);
%         model.grid.y = model.grid.y(1:Mx_HR/Mx_LR:end);
%         model.advection.meth_anti_alias = 'none';
%         model = init_grid_k (model);
%         
%         w = w(1:Mx_HR/Mx_LR:end,1:Mx_HR/Mx_LR:end,:,:);
%         fft_w = fft2(w);
%         %         T_adv = T_adv(1:Mx_HR/Mx_LR:end,1:Mx_HR/Mx_LR:end,:,:);
%         %         fft_buoy_part = fft2(T_adv);
%     else
%         warning('the resolution is already too coarse');
%     end
% end
% 
% % fft_w = SQG_large_UQ(model, fft_buoy_part);
% % w = real(ifft2(fft_w));
% grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
% vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% % vort = vorticity_perso(model.grid, w);
% figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
% fct_plot_velocity(model,w);

%% Residual
w_residual = w_HR - w_smooth;
% fct_plot_velocity(model,w_residual);

w = w_residual;

time_step = 0.5 * 24 * 3600; % half a day
s =size(w);
N_t = s(4);
MX = s(1:2);
M=prod(MX);
w = reshape(w,[M 2 N_t]);

cov_tot = zeros([2 2 N_t]);
for k_t_final=1:N_t
    w_final = w(:,:,k_t_final);
    % w_final = permute(w(:,:,end),[ 1 3 2 ]);
    
    cov_temp = nan([2 2 k_t_final]);
    for k_t=k_t_final:-1:1
        cov_temp(:,:,k_t_final-k_t+1) = 1/M * w(:,:,k_t)' * w_final;
    end
%     figure(1);plot(squeeze(cov_temp(1,1,:)));drawnow
% %     cov_temp = k_t_final * cov_temp;
% %     if k_t_final>1
% %         cov_temp(:,:,1:end-1)=
% %     end
    cov_tot(:,:,1:k_t_final) = ...
        cov_tot(:,:,1:k_t_final) ...
        + cov_temp;
end
weight = 1./permute(N_t:-1:1,[1 3 2]); 
cov_tot = bsxfun(@times, weight, cov_tot);


%% Plots


% load([model.folder.folder_simu '/stat_temp_smooth.mat'],'a','tau_a','time_day','cov_tot');
N_t = size(cov_tot,3);

% time = time_step * (0:(N_t-1));
time_day = 0.5 * (0:(N_t-1));
figure(1);plot(time_day,squeeze(cov_tot(1,1,:)));drawnow
hold on;plot(time_day,squeeze(cov_tot(2,2,:)));drawnow
hold on;plot(time_day,zeros(size(time_day)),'r');drawnow
title('Temporal covariance')


spectrum = abs(fft(cov_tot,[],3));
% spectrum = spectrum(:,:,1:(N_t/2));
% figure(2);loglog(squeeze(spectrum(1,1,:)));drawnow
% hold on;loglog(squeeze(spectrum(2,2,:)));drawnow
figure(2);plot((squeeze(spectrum(1,1,:))));drawnow
hold on;plot((squeeze(spectrum(2,2,:))));drawnow
hold off;
ax=axis;ax(3)=0;axis(ax);
title('Spectrum')

cov_subsample = cov_tot(:,:,1:6:end);

figure(3);plot(time_day(1:6:end),squeeze(cov_subsample(1,1,:)));drawnow
hold on;plot(time_day(1:6:end),squeeze(cov_subsample(2,2,:)));drawnow
hold on;plot(time_day,zeros(size(time_day)),'r');drawnow
hold off
title('Temporal covariance subsample')

spectrum_subsample = abs(fft(cov_subsample,[],3));
% spectrum_subsample = spectrum_subsample(:,:,1:(N_t/2));
% figure(2);loglog(squeeze(spectrum(1,1,:)));drawnow
% hold on;loglog(squeeze(spectrum(2,2,:)));drawnow
figure(4);plot((squeeze(spectrum_subsample(1,1,:))));drawnow
hold on;plot((squeeze(spectrum_subsample(2,2,:))));drawnow
ax=axis;ax(3)=0;axis(ax);
hold off;
title('Spectrum of subsample signal')


% figure(2);plot(time_day,log(squeeze(cov_tot(1,1,:))));drawnow


%     w_final = w(:,:,end);
%     % w_final = permute(w(:,:,end),[ 1 3 2 ]);
%     
%     cov = nan([2 2 N_t]);
%     for k_t=N_t:-1:1
%         cov(:,:,k_t) = 1/M * w(:,:,k_t)' * w_final;
%     end

a = sum(cov_tot(:,:,1:20),3)*time_step
% a = sum(cov_tot,3)*time_step

var_tot = cov_tot(:,:,1);
sigma_tot = chol(var_tot,'lower');
cov_norm = nan(size(cov_tot));
for t=1:N_t
    cov_norm(:,:,t) = sigma_tot \ ( cov_tot(:,:,t) / (sigma_tot') );
end
tau_a = sum(cov_norm(:,:,1:20),3)*time_step/(3600*24)

figure(5);plot(time_day,squeeze(cov_norm(1,1,:)));drawnow
hold on;plot(time_day,squeeze(cov_norm(2,2,:)));drawnow
title('Temporal covariance normalized')

%% Save

save([model.folder.folder_simu '/stat_temp_smooth.mat'],'a','tau_a','time_day','cov_tot');


% a = sum(cov_tot(:,:,1:20),3)*time_step

keyboard;