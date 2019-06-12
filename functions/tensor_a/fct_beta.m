function beta_res = fct_beta(X_i_input, model_input)
% Compute \beta_{i,n} (x) = 2 (x_i-z^i_m)  
% = \nabla ( \gamma ( \| x - z_m \|^2 ) )
% where z^i_m are the coordinates of the observations points and
% \gamma ( \| x - z_m \|^2 ) is the prior kriging covariance
% 

persistent model X_i beta

%% Check if the result is already computed
if  ( isempty(model) && isempty(X_i) && isempty( beta))
% if ~ ( exist('model','var') && exist('X_i','var') && exist( 'beta','var'))
    update_pers = true;
else
    update_pers = ~ ( comp_struct(model_input,model) && ...
        isequal(X_i, X_i_input)  );
end

if update_pers
    % X_m  % d x n_x_m
    % X_i  % n_x x d
    model=model_input;clear model_input;
    X_i=X_i_input;    clear X_i_input;
    
%     X_i_input = permute(double(X_i_input), [ 2 3 1]); % d x 1 x n_x
%     beta = bsxfun( @minus, X_i_input, ...
%         permute(double(model.obs.x), [ 2 1])); % d x n_x_m x n_x
    beta = bsxfun( @minus, ...
        permute(double(X_i), [ 2 3 1]), ...
        permute(double(model.obs.x), [ 2 1])); % d x n_x_m x n_x
%     clear X_i_input;
    norm2_k_betat = sum( beta.^2 ,1); % 1 x n_x_m x n_x
    g_prime = gamma_cov ( norm2_k_betat , model,1); % 1 x n_x_m x n_x
    clear norm2_k_betat
    beta = 2 * bsxfun( @times, beta, g_prime); % d x n_x_m x n_x
    clear g_prime
end
beta_res=beta;