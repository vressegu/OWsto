function [salt,temp,model] = plot_read_data_gula_tracer_multitime_2()
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

model.type_data = 'Gula';
model.filtering.smoothing= true;
init_model
vectname_data = {'u','v','salt','temp'};

one_day = 24 * 3600;
one_week = 7 * one_day;

N_t = 3;


%% Plots tracer
% N_t = size(salt,4)-1;

load([model.folder.folder_simu '/tracer_multi_time.mat'],...
    'temp','salt',...
    'model');


marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t))
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t))
    drawnow;
end

% clear temp salt

%% Plots of Lagrangian source of tracer : are they really tracers?
% N_t = size(salt,4)-1;

load([model.folder.folder_simu '/tracer_multi_time.mat'],...
    'source_temp','source_salt',...
    'model');

marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

std_salt = salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,:);
std_salt = std(std_salt(:));

cax = [-1 1]*20;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,source_salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t)/(std_salt/one_week),...
        cax)
    %     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end


std_temp = temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,:);
std_temp = std(std_temp(:));

for k_t = 1:(N_t+1)
    fct_plot_f(model_,source_temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t)/(std_temp/one_week),...
        cax)
    drawnow;
end

clear source_temp source_salt

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

clear growth_rate_temp growth_rate_salt


%% Plots tracer gradient growth rate HR

load([model.folder.folder_simu '/tracer_multi_time.mat'],...
    'growth_rate_temp_HR','growth_rate_salt_HR',...
    'model');

% N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

cax = [-1 1]*100;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_salt_HR(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_temp_HR(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
    drawnow;
end

clear growth_rate_temp_HR growth_rate_salt_HR


end
