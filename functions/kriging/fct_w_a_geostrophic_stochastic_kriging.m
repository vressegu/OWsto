function model = fct_w_a_geostrophic_stochastic_kriging(recompute_kriging)

%% Algorithnm parametrisation
model.error_obs_h=true;
model.kriging.param_estimated_by_MLE=false;
model.plot_high_kriging=false;
% recompute_kriging=true;
% recompute_u=true;
model.model_order=-1;

% recompute_kriging=false;
model.obs_local=false;
model.plot_dot=false;

% model.plot_PDE=false;
% model.nb_mesh_refin=5;
% model.plot_h_w=false;
% % nb_contour=100;

dx_under_trace = 10e3; % 10 km
model.l_pixel_nugget=dx_under_trace;
% model.l_pixel_nugget=10e3;
% model.l_pixel_nugget=inf;
% % model.dt_a=7*24*3600; % 1 week
% model.dt_a=24*3600; % 1 day
% % model.dt_a=3600; % 1 hour
% % model.dt_a=60; % 1 minute
% % model.dt_a=0;

% Coriolis parameter
model.coriolis.f_model='f_plane';
model.coriolis.earth_radius = 6.4e6;
angle_grid = 45/360*2*pi; % rad
OMEGA = 2*pi/(24*60*60); % rad.s^-1
% OMEGA = 1/(24*60*60);
model.coriolis.f0 = 2* OMEGA * sin( angle_grid );
% coriolis.beta = 2* OMEGA * 2*pi/earth_radius * cos( angle_grid );

model.physical_constant.rho=1e3;
model.physical_constant.g=9.81;
%% Get usefull parameters
% load('data/simu2008_30.mat','dx','dy');
% model.dX=[dx dy];

%% Kriging
if model.obs_local || recompute_kriging
    [model.kriging,model.obs.x,model.obs.h]=kriging_on_high(model);
    if recompute_kriging
        model_save=model;
        save('model_global','model_save');
    end
else
    model_ref=model;
    load('model_global.mat');
    % model.l_pixel_nugget=dx_under_trace/10;
    model_ref.kriging=model_save.kriging;
    model_ref.obs.x=model_save.obs.x;
    model_ref.obs.h=model_save.obs.h;
    model=model_ref;
%     %     model.kriging=model_save.kriging;
%     %     model.obs.x=model_save.obs.x;
%     %     model.obs.h=model_save.obs.h;
%     model.l_pixel_nugget=inf;
end