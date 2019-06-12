function grad_h = fct_grad_h(XP,model)
% Compute (semi-analiticaly) the gradient of the SSH kriging posterior mean
%

[n_x,d]=size(XP);
n_x_m=size(model.obs.x,1);
threshold=8e7;
big_data = n_x_m*n_x>threshold;
if big_data
    grad_h = nan(d,n_x);
    n_block=floor(threshold/n_x_m);
    for k=1:ceil(n_x*n_x_m/threshold)
        grad_h(:,(k-1)*n_block+1:min(k*n_block,end)) =  ...
           fct_grad_h(XP((k-1)*n_block+1:min(k*n_block,end),:),model); % d x n_x
    end
    return
end

%% Compute the gradient of \lambda : the linear relation on observed and interpolate data in Kriging
grad_lambda = fct_a_V(XP,model); % n_x_m x d x n_x
% Compute a kind of Cholesky decomposition of the inhomogeneous par of a :
% chol ( \gamma(x^{obs}) ) ^{-1} )^t \beta^t J^t 

% Inverse the (triangular) Cholesky decomposition of the kriging covariance
opts_LS.UT = true;
grad_lambda = multi_linsolve( fct_chol_gamma(model), ...
     grad_lambda, opts_LS);% n_x_m x d x n_x
% = ( \gamma(x^{obs}) ) ^{-1} )^t \beta^t J^t 

grad_lambda = multitrans(grad_lambda);% d x n_x_m x n_x
% = J \beta ( \gamma(x^{obs}) ) ^{-1} )
grad_lambda(1,:,:) = - grad_lambda(1,:,:);
grad_lambda = grad_lambda([2 1],:,:);% d x n_x_m x n_x
% = \beta ( \gamma(x^{obs}) ) ^{-1} )

%% Compute the part of the gradient corresponding to the "innovation"
% Value of the prior mean at the measurements locations
h_prior_obs = h_prior(model.obs.x, model.kriging.order);
% (Measurements values) - (Value of the prior mean at the measurements locations)
% ( = "innovation" )
delta_h = model.obs.h - h_prior_obs; % n_x_m
% Corresponding part in gradient of the SSH kriging posterior mean
grad_h = squeeze(sum( bsxfun( @times, multitrans(grad_lambda), delta_h ) ,1) ); % d x n_x
clear grad_lambda delta_h

%% Compute the gradient of the prior SSH mean
switch model.kriging.order
    case {-1,0}
        grad_h_prior = 0;
    case 1
        error('not coded yet');
    otherwise
       error('not coded yet');
end

%% Gather the two parts
grad_h = grad_h + grad_h_prior; % d x n_x

%% Sub function
    function h_prior_res = h_prior(x,order)
        switch order
            case -1
                h_prior_res = 0;
            case 0
                error('not coded yet');
            case 1
                error('not coded yet');
            otherwise
                error('not coded yet');
        end
        
    end
end
