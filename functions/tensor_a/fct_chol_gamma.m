function chol_gamma_res = fct_chol_gamma(model_input)
% Compute the Cholesky decomposition of the height kriging covariance
%

persistent model chol_gamma

%% Check if the result is already computed
if  ( isempty(model) && isempty( chol_gamma))
% if ~ ( exist('model','var') && exist( 'chol_gamma','var'))
    update_pers = true;
else
    update_pers = ~ ( comp_struct(model_input,model) );
end

if update_pers
    model=model_input;clear model_input;
    
    gamma = stk_make_matcov (model.kriging, model.obs.x); % n_x_m x n_x_m
%     gamma = stk_make_matcov (model.kriging, model.obs.x, model.obs.x); % n_x_m x n_x_m
%     gamma = stk_make_matcov (model.kriging, model.obs.x); % n_x_m x n_x_m
    chol_gamma = chol(gamma); % n_x_m x n_x_m
    if size(double(model.obs.x),2) ~= 2
        error('wrong size');
    end
end
chol_gamma_res=chol_gamma;

end