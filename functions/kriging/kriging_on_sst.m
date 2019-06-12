function [kriging_model,xi,zi] = kriging_on_sst(model)
% Used of the example kb03 of the toolbox stk
% to krige an SSH field

% error_obs_sst=model.error_obs_sst;
param_estimated_by_MLE=true;
% param_estimated_by_MLE=model.kriging.param_estimated_by_MLE;

covariance_type ='gaussienne'; % 'gaussienne' 'exponentielle'
% % model_order=2; % quadratic mean
% model_order=1; % intrinsec kriging : unknown  affine mean
model_order=model.kriging.order; % intrinsec kriging : unknown constant mean
% model_order= - 1; % ordinnary krigin : constant mean =0

plot_sst_kriging=model.plot_sst_kriging;
unity_approx=false;
% if nargin <2
% model.kriging.param_estimated_by_MLE =false;
% % model.kriging.param_estimated_by_MLE =true;
% end

% if nargin < 1 || error_obs_sst
%     measurement_error=1e-2;
% else
measurement_error=0;
% end

%% SSH of reference


DIM = 2;
% % L = 6.371e6* [ 0.0021 2.1823 ; -0.4841 0.4841];
% % L = 1e6* [ 0.0021 2.1823 ; -0.4841 0.4841];
% % L = 1e6* [ 0.0021 1.0890 ; -0.9711 0.9711];
% L = 1e6* [ 0.0021 1.0890 ; 0.0019 1.9440];
%
% if model.obs_local
%     ref_L = repmat(L(:,1),[1 2]);
%     L = L(:,2)-L(:,1);
% %     L = [ 0.5 0.6 ; 0.45 0.5 ] .* repmat(L,[1 2]);
%     L = [ 0.4 0.7 ; 0.45 0.6 ] .* repmat(L,[1 2]);
% %     L = [ 0.4 0.58 ; 0.45 0.5 ] .* repmat(L,[1 2]);
% % %     L = [ 0.4 0.58 ; 0.3 0.5 ] .* repmat(L,[1 2]);
% % %     % L = [ 0.4 0.6 ; 0.4 0.6 ] .* repmat(L,[1 2]);
%     L = L + ref_L;
% end
%
% BOX = L';
BOX=model.kriging.BOX;
% marge=1e5;
% BOX(1,:)=BOX(1,:)-marge/2;
% BOX(2,:)=BOX(2,:)+marge/2;

f  = @(x) (sst_data(x,unity_approx,BOX));
if measurement_error==0
    f_obs = f;
else
    error(['not coded yet. the measurement noise need to be periodic' ...
        ' (cf function kriging_on_high)']);
end
% f_obs  = @(x) (sst_data(x,unity_approx,BOX) + measurement_error * sqrt(2) * randn) ;

%% COMPUTE AND VISUALIZE THE FUNCTION ON A 80 x 80 REGULAR GRID

if plot_sst_kriging
    % Size of the regular grid
    NT = 30^2;
    % NT = 80^2;
    
    % The function stk_sampling_regulargrid() does the job of creating the grid
    xt = stk_sampling_regulargrid (NT, DIM, BOX);
    
    % Compute the corresponding responses
    zt = stk_feval (f, xt);
    
    CONTOUR_LINES = 40; % number of levels in contour plots
    DOT_STYLE = {'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4};
    % Since xt is a regular grid, we can do a contour plot
    figure;
    stk_subplot (2, 2, 1);  stk_plot2d (@contour, xt, f, CONTOUR_LINES);
    % stk_subplot (2, 2, 1);  stk_plot2d (@contour, xt, f, CONTOUR_LINES);
    axis (BOX(:));  stk_title ('function to be approximated');
    axis equal
end

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

%% EVALUATE THE FUNCTION ON A "MAXIMIN LHS" DESIGN

% xi = stk_dataframe (fct_grid_satellite(L,0,200e3,1000e3)');
% xi = fct_grid_satellite(L,pi/8)';
% xi = stk_dataframe (fct_grid_satellite(BOX)');

% NI = 10^2;
NI = 29^2;
% NI = (floor(100/model.sst_n_sub))^2;
% xi = stk_sampling_regulargrid (NI, DIM, BOX);
% xi = stk_sampling_randomlhs (NI, DIM, BOX_period);
% NI = size(xi,1);

load('data/simu2008_30.mat','Mx','My','dx','dy');
% [My,Mx,~]=size(ssh);clear ssh
My=My-2;
BOX_period = [ dx*[0 Mx]  ; dy*[0 My] ]';
% BOX_period = [ dx*[1 Mx+1]  ; dy*[1 My+1] ]';
marge=1e5;
BOX_period(1,:)=BOX_period(1,:)-marge/2;
BOX_period(2,:)=BOX_period(2,:)+marge/2;
xi = stk_sampling_randomlhs (NI, DIM, BOX_period);
 
zi = stk_feval (f_obs, xi);
xi = [xi ; xi+[dx*(Mx) 0 ] ; xi-[dx*(Mx) 0 ] ]; 
zi = [zi ; zi; zi ]; 
clear Mx My 
% zi = stk_feval (f_obs, xi);

if plot_sst_kriging
    % Add the design points to the first plot
    hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
end

%% ESTIMATE THE PARAMETERS OF THE COVARIANCE FUNCTION
% nugget=(measurement_error>0);
nugget=true;
if kriging_model.param_estimated_by_MLE
    % Compute an initial guess for the parameters of the Matern covariance (param0)
    % and a reasonable log-variance for a small "regularization noise"
    %%
    %     sigma2=2^2;
    %     nugget_val=1^2;
    %     rho=1/180*pi*model.coriolis.earth_radius;% 1 degre
    %
    %     param0(1)=log(sigma2);
    %     param0(2)=-log(rho);
    %     kriging_model.lognoisevariance=log(nugget_val);
    %     zi
    %     keyboard;
    [param0, kriging_model.lognoisevariance] = stk_param_init (kriging_model, xi, zi, BOX,nugget);
    %%
    
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
    
else 
    error('not coded yet. Need reference on the subjest');
end


%% CARRY OUT KRIGING PREDICITION AND VISUALIZE

if plot_sst_kriging
    % Here, we compute the kriging prediction on each point of the grid
    zp = stk_predict (kriging_model, xi, zi, xt);
    
    % x=double(xt);
    % z=double(zt(:,1));
    % save('data_test','x','z');
    
    % Display the result using a contour plot, to be compared with the contour
    % lines of the true function
    stk_subplot (2, 2, 2);  stk_plot2d (@contour, xt, zp.mean, CONTOUR_LINES);
    % stk_subplot (2, 2, 2);  stk_plot2d (@contour, xt, zp.mean, CONTOUR_LINES);
    tsc = sprintf ('approximation from %d points', NI);  hold on;
    plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
    hold off;  axis (BOX(:));  stk_title (tsc);
    axis equal
    
    %% VISUALIZE THE ACTUAL PREDICTION ERROR AND THE KRIGING STANDARD DEVIATION
    
    stk_subplot (2, 2, 3);  stk_plot2d (@pcolor, xt, log (abs (zp.mean - zt)));
    %     stk_subplot (2, 2, 3);  stk_plot2d (@pcolor, xt, log (abs (zp.mean - zt)));
    [cmin,cmax] = caxis;
    %     map =colormap;
    hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
    hold off;  axis (BOX(:));  stk_title ('true approx error (log)');
    axis equal
    
    stk_subplot (2, 2, 4);  stk_plot2d (@pcolor, xt, 0.5 * log (zp.var));
    %     stk_subplot (2, 2, 4);  stk_plot2d (@pcolor, xt, 0.5 * log (zp.var));
    caxis([cmin,cmax]);
    %     colormap(map);
    hold on;  plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
    hold off;  axis (BOX(:));  stk_title ('kriging std (log)');
    axis equal
    
end

%%
% kriging_model
% save('model_exp.mat','kriging_model');
