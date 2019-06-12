function kriging_model = krige_spatially_ML(model,ssh_on_trace,grid_trace)
% Krige data with the example kb03 of the toolbox stk
% to krige an SSH field

% error_obs_h=model.error_obs_h;
param_estimated_by_MLE=model.kriging.param_estimated_by_MLE;

if isfield(model.kriging,'covariance_type')
    covariance_type = model.kriging.covariance_type;
else
    covariance_type ='gaussienne'; % 'gaussienne' 'exponentielle' 'matern'
end

model_order=model.kriging.model_order;
% % model_order=2; % quadratic mean
% model_order=1; % intrinsec kriging : unknown  affine mean
% model_order= 0; % ordinnary krigin : unknown  constant mean
% model_order= - 1; % ordinnary krigin : constant mean =0

% plot_high_kriging=model.plot_high_kriging;
% unity_approx=false;

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
    case 'matern'
        kriging_model = stk_model ('stk_materncov_iso', DIM);
end
kriging_model.param_estimated_by_MLE=param_estimated_by_MLE;

% As a default choice, a constant (but unknown) mean is used,
% i.e.,  kriging_model.order = 0.
kriging_model.order = model_order;
% kriging_model.order = 1;  %%% UNCOMMENT TO USE A LINEAR TREND %%%
% kriging_model.order = 2;  %%% UNCOMMENT TO USE A "FULL QUADRATIC" TREND %%%
kriging_model.BOX=BOX;

%% EVALUATE THE FUNCTION ON A "MAXIMIN LHS" DESIGN

% % % xi = stk_dataframe (fct_grid_satellite(L,0,200e3,1000e3)');
% % % xi = fct_grid_satellite(L,pi/8)';
% % xi = stk_dataframe (fct_grid_satellite(BOX)');
% xi = stk_dataframe (grid_trace');
xi = stk_dataframe (grid_trace);
NI = size(xi,1);

zi = stk_dataframe (ssh_on_trace);
% zi = stk_feval (f_obs, xi);
% % zi = stk_feval (f_obs, xi);
% % close;find(isnan(double(zi)))
% % keyboard;

% if plot_high_kriging
%     % Add the design points to the first plot
%     hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
% end

%% ESTIMATE THE PARAMETERS OF THE COVARIANCE FUNCTION
nugget=model.error_obs_h;
% nugget=true;
if kriging_model.param_estimated_by_MLE
    
    
    if nugget
%         if kriging_model.lognoisevariance < ...
%                 2*log(model.measurement_error_h * sqrt(2))
            kriging_model.lognoisevariance = ...
                2*log(model.measurement_error_h * sqrt(2));
%         end
    end
    
    
    % Compute an initial guess for the parameters of the Matern covariance (param0)
    % and a reasonable log-variance for a small "regularization noise"
    [param0, kriging_model.lognoisevariance] = ...
        stk_param_init (kriging_model, xi, zi, BOX);
%         stk_param_init (kriging_model, xi, zi, BOX,nugget);
    
    switch covariance_type
        case {'gaussienne','exponentielle'}
            ampli = exp(1/2*param0(1))
            rho = exp(-param0(2))
            nugget_value = exp(kriging_model.lognoisevariance)
        case 'matern'
            ampli = exp(1/2*param0(1))
            nu = exp(param0(2));
            beta = 2*nu-1
            rho = exp(-param0(3))
            nugget_value = exp(kriging_model.lognoisevariance)
            
%             beta = 1;
%             nu = (beta+1)/2;
%             param0(2) = log(nu);
%             warning('value of nu changed')
    end
    

%     if nugget
% %         if kriging_model.lognoisevariance < ...
% %                 2*log(model.measurement_error_h * sqrt(2))
%             kriging_model.lognoisevariance = ...
%                 2*log(model.measurement_error_h * sqrt(2));
% %         end
%     end
    


    kriging_model.param = stk_param_estim (kriging_model, xi, zi, param0);
    
    switch covariance_type
        case {'gaussienne','exponentielle'}
            ampli = exp(1/2*kriging_model.param(1))
            rho = exp(-kriging_model.param(2))
            nugget_value = exp(kriging_model.lognoisevariance)
        case 'matern'
            ampli = exp(1/2*kriging_model.param(1))
            nu = exp(kriging_model.param(2));
            beta = 2*nu-1
            rho = exp(-kriging_model.param(3))
            nugget_value = exp(kriging_model.lognoisevariance)
    end
    
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