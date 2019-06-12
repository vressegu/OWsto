%%%%%%%%%%%%%%%%%%%%
%%% Main
%%%%%%%%%%%%%%%%%%%%
init;

type_data =  'Loop_current';
% type_data =  'Gula_Summer';
% type_data =  'Gula';
% type_data =  'Gula_kriging';
smoothing = true;
subsampling = false;
% Mx_LR = 8; % filter of 30 km -> resolution ~ 100km
Mx_LR = 64; % filter of 4 km -> resolution ~ 10km

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
    case 'Loop_current'
        day = 30;
        
%         init_model
        model.grid.dX = 500e3/1024 * [1 1];
        model.grid.MX = [1536 512];
        dt = 24 * 3600;
        model.grid.x = model.grid.dX(1) * ( 0:(model.grid.MX(1)-1));
        model.grid.y = model.grid.dX(2) * ( 0:(model.grid.MX(2)-1));
        x=model.grid.x;
        y=model.grid.y;
        w = nan([model.grid.MX 2 1]);
        load(['/Volumes/WD_Ressegui/These/from_tour/data_brest/' ...
            'data_molemaker/filtered_time/summer/u_filtered.mat'],...
            'data_smooth');
        w(:,:,1)=real(data_smooth(2:end,1:512,1,day)); clear data_smooth
        load(['/Volumes/WD_Ressegui/These/from_tour/data_brest/' ...
            'data_molemaker/filtered_time/summer/v_filtered.mat'],...
            'data_smooth');
        w(:,:,2)=real(data_smooth(2:(end-1),2:513,1,day)); clear data_smooth
%         '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_molemaker/filtered_time/summer/temp_filtered.mat'
        
%         model.advection.meth_anti_alias = 'none';
%         model = init_grid_k (model);
    case 'Gula'
        day = 30;
%         day = 1;
        time_day = 'night';
        [w,model] = read_data_gula(day,time_day);
        model.advection.meth_anti_alias = 'none';
        model = init_grid_k (model);
    case 'Gula_Summer'
        day = 150;
%         day = 1;
        time_day = 'night';
        [w,model] = read_data_gula(day,time_day);
        model.advection.meth_anti_alias = 'none';
        model = init_grid_k (model);
    case 'Gula_kriging'
        % kriged SSH from very high reoslution model output + noise
        load('images/krige_space_then_time/zeta_filtered.mat',...
            'data','model');
        day = 30;
        ssh=data(:,1:end-1,:,day); clear data;
        model.grid_HR = model.grid;
        model.grid = model.grid_LR;
        model = rmfield(model,'grid_LR');
        model.grid.y(end) = [];
        model.grid.MX(2) = model.grid.MX(2)-1;
        
        % Physical parameters
        model = fct_physical_param(model);
        
        % Coriolis
        nabla_ssh = gradient_mat_2( ssh,model.grid.dX);
        w(:,:,1,:)= - nabla_ssh(:,:,2,:);
        w(:,:,2,:)= nabla_ssh(:,:,1,:);
        clear nabla_ssh;
        w = model.physical_constant.g / model.physical_constant.f0 * w;
        warning('coriolis not accurate');
        
        smoothing = false;
end
model.type_data = type_data;
model.filtering.smoothing = smoothing;
model.filtering.subsampling = subsampling;

model.plot.day_plot=1;

init_model

grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% vort = vorticity_perso(model.grid, w);
figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
fct_plot_velocity(model,w);

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
grad_w = gradient_mat_2(permute( w,[ 1 2 4 3]),model.grid.dX);
vort = grad_w(:,:,1,2) - grad_w(:,:,2,1);
% vort = vorticity_perso(model.grid, w);
figure;imagesc(model.grid.x,model.grid.y,vort');axis xy;axis equal
fct_plot_velocity(model,w);

%% Mixing criterion

[Q,sigma_w,r_OW] = fct_okubo_weiss(model,w);

% r = fct_lapeyre(model,w);


%% My criterion

[alpha2,states] = fct_mix_sto(model,r_OW.^2,sigma_w);

keyboard;