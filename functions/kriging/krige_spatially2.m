function [kriging_model] = krige_spatially2(model)
% Krige data with the example kb03 of the toolbox stk
% to krige an SSH field



%% Algorithnm parametrisation
% model.error_obs_h=true;
% model.kriging.param_estimated_by_MLE=false;
% model.plot_high_kriging=false;
% % recompute_kriging=true;
% % recompute_u=true;
model.model_order=-1;
warning('model order is set to zero')

% % recompute_kriging=false;
% model.obs_local=false;
% model.plot_dot=false;

% error_obs_h=model.error_obs_h;
param_estimated_by_MLE=model.kriging.param_estimated_by_MLE;

covariance_type ='gaussienne'; % 'gaussienne' 'exponentielle'
% % model_order=2; % quadratic mean
% model_order=1; % intrinsec kriging : unknown  affine mean
model_order=model.model_order; % intrinsec kriging : unknown constant mean
% model_order= - 1; % ordinnary krigin : constant mean =0

% plot_high_kriging=model.plot_high_kriging;
unity_approx=false;

%% SSH of reference


DIM = 2;
BOX=model.grid.BOX;
BOX=BOX';
% BOX = L';

%% COMPUTE AND VISUALIZE THE FUNCTION ON A 80 x 80 REGULAR GRID

% if plot_high_kriging
%     
%     % Size of the regular grid
%     NT = 30^2;
%     % NT = 80^2;
%     
%     % The function stk_sampling_regulargrid() does the job of creating the grid
%     xt = stk_sampling_regulargrid (NT, DIM, BOX);
%     
%     % Compute the corresponding responses
%     zt = stk_feval (f, xt);
%     
%     CONTOUR_LINES = 40; % number of levels in contour plots
%     DOT_STYLE = {'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4};
%     % Since xt is a regular grid, we can do a contour plot
%     figure;
%     stk_subplot (2, 2, 1);  stk_plot2d (@contour, xt, f, CONTOUR_LINES);
%     % stk_subplot (2, 2, 1);  stk_plot2d (@contour, xt, f, CONTOUR_LINES);
%     axis (BOX(:));  stk_title ('function to be approximated');
%     axis equal
% end

%% CHOOSE A KRIGING (GAUSSIAN PROCESS) kriging_model

switch covariance_type
    case 'gaussienne'
        kriging_model = stk_model ('stk_gausscov_iso', DIM);
    case 'exponentielle'
        kriging_model = stk_model ('stk_expcov_iso', DIM);
end
kriging_model.param_estimated_by_MLE=param_estimated_by_MLE;

% As a default choice, a constant (but unknown) mean is used,
% i.e.,  kriging_model.order = 0.
kriging_model.order = model_order;
% kriging_model.order = 1;  %%% UNCOMMENT TO USE A LINEAR TREND %%%
% kriging_model.order = 2;  %%% UNCOMMENT TO USE A "FULL QUADRATIC" TREND %%%
kriging_model.BOX=BOX;

% %% EVALUATE THE FUNCTION ON A "MAXIMIN LHS" DESIGN
% 
% % % xi = stk_dataframe (fct_grid_satellite(L,0,200e3,1000e3)');
% % % xi = fct_grid_satellite(L,pi/8)';
% % xi = stk_dataframe (fct_grid_satellite(BOX)');
% xi = stk_dataframe (grid_trace');
% NI = size(xi,1);
% 
% zi = stk_dataframe (ssh_on_trace);
% % zi = stk_feval (f_obs, xi);
% % % zi = stk_feval (f_obs, xi);
% % % close;find(isnan(double(zi)))
% % % keyboard;
% 
% % if plot_high_kriging
% %     % Add the design points to the first plot
% %     hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
% % end

%% ESTIMATE THE PARAMETERS OF THE COVARIANCE FUNCTION
nugget=model.error_obs_h;
% nugget=(measurement_error>0);
% nugget=true;
if kriging_model.param_estimated_by_MLE
    % Compute an initial guess for the parameters of the Matern covariance (param0)
    % and a reasonable log-variance for a small "regularization noise"
    [param0, kriging_model.lognoisevariance] = stk_param_init (kriging_model, xi, zi, BOX,nugget);
    
    % kriging_model.lognoisevariance =log(1e-2);
    
    % % Alternative: user-defined initial guess for the parameters of
    % % the Matern covariance (see "help stk_materncov_aniso" for more information)
    % SIGMA2 = var (zi);
    % NU     = 2;
    % RHO1   = (BOX(2,1) - BOX(1,1)) / 10;
    % RHO2   = (BOX(2,2) - BOX(1,2)) / 10;
    % param0 = log ([SIGMA2; NU; 1/RHO1; 1/RHO2]);
    % kriging_model.lognoisevariance = 2 * log (1e-5);
    
    kriging_model.param = stk_param_estim (kriging_model, xi, zi, param0);
else % from GILLE & KELLY 1996 (Southern ocean)
    sigma2=1.20e-2; % 120 cm^2 (ref abderahim) (2 fois + que estime MLE:4.15e-2 )
    rossby_radius=85e3; % 83 km (proche de MLE : 52.327 km )
    %     rossby_radius=85e3; % 83 km (3 fois moins)
    %     nugget_val=(sqrt(2)*measurement_error)^2; % measurement error
    %     nugget_val=4e-3;
    % meme ref (instrument error, orbit error, and residual atmospheric and electromagnetic (EM) biais error
    nugget_val=(5.6e-2)^2; % meme ref (ecart type 10 fois+ grand que estime MLE : 0.63cm)
    kriging_model.param(1)=log(sigma2);
    kriging_model.param(2)=-log(rossby_radius);
    kriging_model.lognoisevariance=log(nugget_val);
end

%% CARRY OUT KRIGING PREDICITION AND VISUALIZE
% 
% if plot_high_kriging
%     % Here, we compute the kriging prediction on each point of the grid
%     zp = stk_predict (kriging_model, xi, zi, xt);
%     
%     % x=double(xt);
%     % z=double(zt(:,1));
%     % save('data_test','x','z');
%     
%     % Display the result using a contour plot, to be compared with the contour
%     % lines of the true function
%     stk_subplot (2, 2, 2);  stk_plot2d (@contour, xt, zp.mean, CONTOUR_LINES);
%     % stk_subplot (2, 2, 2);  stk_plot2d (@contour, xt, zp.mean, CONTOUR_LINES);
%     tsc = sprintf ('approximation from %d points', NI);  hold on;
%     plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
%     hold off;  axis (BOX(:));  stk_title (tsc);
%     axis equal
%     
%     %% VISUALIZE THE ACTUAL PREDICTION ERROR AND THE KRIGING STANDARD DEVIATION
%     
%     stk_subplot (2, 2, 3);  stk_plot2d (@pcolor, xt, log (abs (zp.mean - zt)));
%     %     stk_subplot (2, 2, 3);  stk_plot2d (@pcolor, xt, log (abs (zp.mean - zt)));
%     [cmin,cmax] = caxis;
%     %     map =colormap;
%     hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
%     hold off;  axis (BOX(:));  stk_title ('true approx error (log)');
%     axis equal
%     
%     
%     stk_subplot (2, 2, 4);  stk_plot2d (@pcolor, xt, 0.5 * log (zp.var));
%     %     stk_subplot (2, 2, 4);  stk_plot2d (@pcolor, xt, 0.5 * log (zp.var));
%     caxis([cmin,cmax]);
%     %     colormap(map);
%     hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
%     hold off;  axis (BOX(:));  stk_title ('kriging std (log)');
%     axis equal
%     
% end

%%
% kriging_model
% save('model_exp.mat','kriging_model');