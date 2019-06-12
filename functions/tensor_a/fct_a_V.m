function V_res = fct_a_V(X_i_input,model_input)
% Compute a kind of Cholesky decomposition of the inhomogeneous par of
% a(x,x)
% V(x) = chol ( \gamma(x^{obs}) ) ^{-1} )^t \beta^t J^t 

%% Check if V is already computed
persistent model X_i V

if  ( isempty(model) && isempty(X_i) && isempty( V))
% if ~ ( exist('model','var') && exist('X_i','var') && exist( 'V','var'))
    update_pers = true;
else
    update_pers = ~ ( comp_struct(model_input,model) && ...
        isequal(X_i, X_i_input)  );
end

if ~ update_pers
    V_res = V;
    return
end
model=model_input;clear model_input;
X_i=X_i_input;clear X_i_input;

%% Inhomogeneous part of V
beta = fct_beta(X_i, model); % d x n_x_m x n_x
chol_gamma = fct_chol_gamma(model);% n_x_m x n_x_m

% Compute the vector directly orthogonal
beta = multi_k_x(beta); % d x n_x_m x n_x

% Inverse the (triangular) Cholesky decomposition of the kriging covariance
opts_LS.LT = true;
V_res = multi_linsolve(chol_gamma', multitrans(beta), opts_LS);% n_x_m x d x n_x
clear beta chol_gamma;

%% Save as persistent variable
V = V_res;

end