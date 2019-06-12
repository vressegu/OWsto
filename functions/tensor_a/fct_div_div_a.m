function div_div_a_res = fct_div_div_a(X_i_input,model_input)
% Compute 1/2 div( ( div( a(x,x) )^t )
%

%% Check if div_div_a is already computed
persistent model X_i div_div_a

if  ( isempty(model) && isempty(X_i) && isempty( div_div_a))
% if ~ ( exist('model','var') && exist('X_i','var') && exist( 'div_div_a','var'))
    update_pers = true;
else
    update_pers = ~ ( comp_struct(model_input,model) && ...
        isequal(X_i, X_i_input)  );
end

if ~ update_pers
    div_div_a_res = div_div_a;
    return
end
model=model_input;clear model_input;
X_i=X_i_input;clear X_i_input;

%% Preliminary calculations
N=size(double(model.obs.x),1);
chol_gamma = fct_chol_gamma(model);% n_x_m x n_x_m
beta = bsxfun( @minus, ...
    permute(double(X_i), [ 2 3 1]), ...
    permute(double(model.obs.x), [ 2 1])); % d x n_x_m x n_x
norm2_k_betat = sum( beta.^2 ,1); % 1 x n_x_m x n_x
    
%% 1st part of div_div_a
g_prime = permute(gamma_cov ( norm2_k_betat , model,1),[2 3 1]); % n_x_m x n_x

% Inverse the (triangular) Cholesky decomposition of the kriging covariance
opts_LS.LT = true;
div_div_a_1 = multi_linsolve(chol_gamma', g_prime, opts_LS);% n_x_m x n_x
clear g_prime
div_div_a_1 = (4 * N * (N-1) ) * sum(div_div_a_1.^2 , 1)' ; % n_x
% div_div_a_1 = (4 * N * (N-1) ) * squeeze(sum(div_div_a_1.^2 , 1)); % n_x

%% 2nd part
g_prime_prime = gamma_cov ( norm2_k_betat , model,2); % 1 x n_x_m x n_x
clear norm2_k_betat
beta = 2 * bsxfun( @times, beta, g_prime_prime); % d x n_x_m x n_x
clear g_prime_prime

% Inverse the (triangular) Cholesky decomposition of the kriging covariance
opts_LS.LT = true;
div_div_a_2 = multi_linsolve(chol_gamma', multitrans(beta), opts_LS);% n_x_m x d x n_x
clear chol_gamma;
div_div_a_2 = squeeze(sum(sum(div_div_a_2.^2 , 2) ,1)); % n_x

%% Gather two parts
div_div_a_res = div_div_a_1 + div_div_a_2;
clear div_div_a_1 div_div_a_2
g=9.81;
delta_t=model.dt_a;
div_div_a_res = - ( delta_t * (g/f_coriolis(model.coriolis)).^2 ) /2 ...
                                         .* div_div_a_res ; % n_x

%% Save as persistent variable
div_div_a = div_div_a_res;

end