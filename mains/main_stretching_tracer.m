%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
init;

type_data =  'Gula';
% type_data =  'Gula_kriging';
smoothing = false;
subsampling = false;
Mx_LR = 8; % filter of 30 km
multitime = true;

model.type_data=type_data;
init_model
day_plot=30;

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
            model.dt = 0.5 * 24 * 3600;
        else
            day = 30;
            %         day = 1;
            time_day = 'night';
            [w,model] = read_data_gula(day,time_day);
            model.advection.meth_anti_alias = 'none';
            model = init_grid_k (model);
        end
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
        if smoothing
            w = w_smooth;
        else
            w = w_HR;
        end
        
        model.dt = 0.5 *24 * 3600;

end
model.type_data = type_data;
model.filtering.smoothing = smoothing;
model.filtering.subsampling = subsampling;
model.plot.day_plot=day_plot;

init_model

grad_w = gradient_mat_2(permute( w(:,:,:,model.plot.day_plot),[ 1 2 4 3]),model.grid.dX);
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% vort = vorticity_perso(model.grid, w);
figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
fct_plot_velocity(model,w(:,:,:,model.plot.day_plot));

%% Smoothing
if model.filtering.smoothing
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
    w(:,:,1) = conv2(w(:,:,1),filter_perso,'same');
    w(:,:,2) = conv2(w(:,:,2),filter_perso,'same');
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
end

%% Subsampling
if model.filtering.subsampling
    Mx_HR = min(model.grid.MX);
%     Mx_HR = model.grid.MX(1);
    if  Mx_LR < Mx_HR
        model.grid.MX = model.grid.MX * (Mx_LR/Mx_HR);
        model.grid.dX = model.grid.dX / (Mx_LR/Mx_HR);
        model.grid.x = model.grid.x(1:Mx_HR/Mx_LR:end);
        model.grid.y = model.grid.y(1:Mx_HR/Mx_LR:end);
        model.advection.meth_anti_alias = 'none';
        model = init_grid_k (model);
        
        w = w(1:Mx_HR/Mx_LR:end,1:Mx_HR/Mx_LR:end,:,:);
        fft_w = fft2(w);
%         T_adv = T_adv(1:Mx_HR/Mx_LR:end,1:Mx_HR/Mx_LR:end,:,:);
%         fft_buoy_part = fft2(T_adv);
    else
        warning('the resolution is already too coarse');
    end
end

% fft_w = SQG_large_UQ(model, fft_buoy_part);
% w = real(ifft2(fft_w));
grad_w = gradient_mat_2(permute( ...
    w(:,:,:,model.plot.day_plot),[ 1 2 4 3]),model.grid.dX);
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% vort = vorticity_perso(model.grid, w);
figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
fct_plot_velocity(model,w(:,:,:,model.plot.day_plot));

%% Mixing criterion

% close all
% warning('some time removed')
% w(:,:,:,33:end)=[];

Q = fct_okubo_weiss(model,w);
% Q = fct_okubo_weiss(model,w(:,:,:,model.plot.day_plot));

[r,sigma_w] = fct_lapeyre(model,w);

save([model.folder.folder_simu '/Lapeyre_OW_mixing_multitime.mat'],...
    'Q','r','sigma_w','model','w');



keyboard;