function [salt,temp,model] = smooth_growth_rate_gula()
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

model.type_data = 'Gula';
model.filtering.smoothing= true;
init_model
vectname_data = {'u','v','salt','temp'};


% Mx_LR = 64; % filter of 4 km -> resolution ~ 10km
Mx_LR = 32; % filter of 8 km -> resolution ~ 20km
% Mx_LR = 8; % filter of 30 km -> resolution 50-100km


one_day = 24 * 3600;
one_week = 7 * one_day;

N_t = 3;

%% Plots tracer gradient growth rate

load([model.folder.folder_simu '/tracer_multi_time.mat'],...
    'growth_rate_temp','growth_rate_salt',...
    'model');

% N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

cax = [-1 1]*100;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
    drawnow;
end

% clear growth_rate_temp growth_rate_salt

%% Smooth growth rate

growth_rate_temp_smooth = nan([model.grid.MX 1 N_t+1]);
growth_rate_salt_smooth = nan( [model.grid.MX 1 N_t+1]);

Mx_HR = min(model.grid.MX);

if  Mx_LR < Mx_HR
    % Filter
    sigma_filter = min(model.grid.dX)/2 * Mx_HR/Mx_LR ;
    model.advection.sigma_filter=sigma_filter;
    filter_perso = design_filter_3(model);
    
    for k_t = 1:1:(N_t+1)
        data_temp=growth_rate_temp(:,:,1,k_t);
        
        % Convolution
        data_temp_smooth = conv2(data_temp,filter_perso,'same');
        growth_rate_temp_smooth(:,:,1,k_t)=data_temp_smooth;
    end
else
    warning('the resolution is already too coarse');
end


for k_t = 1:1:(N_t+1)
    data_temp=growth_rate_salt(:,:,1,k_t);
    
    if  Mx_LR < Mx_HR
        % Filter
        sigma_filter = min(model.grid.dX)/2 * Mx_HR/Mx_LR ;
        model.advection.sigma_filter=sigma_filter;
        filter_perso = design_filter_3(model);
        
        % Convolution
        data_temp_smooth = conv2(data_temp,filter_perso,'same');
    else
        warning('the resolution is already too coarse');
    end
    growth_rate_salt_smooth(:,:,1,k_t)=data_temp_smooth;
end


save([model.folder.folder_simu '/smooth_tracer_growth_rate.mat']);


%% Plots tracer gradient growth rate smooth



% N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

cax = [-1 1]*100;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_salt_smooth(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_temp_smooth(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
    drawnow;
end

keyboard  
    
end
