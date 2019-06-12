function [coef_homogen_a, tra, lambda_min, lambda_max] = fct_vp_a(X_i,model)
% Compute singular value of the quadratic covariation of the geostrophic velocity
% a(x,x)

%% Compute a
a = fct_a(X_i,model);% n_x x d x d
% Homogeneous part of a
l_pixel=model.l_pixel_nugget;
nugget = exp(model.kriging.lognoisevariance);
coef_homogen_a = ( nugget / (l_pixel/2)^2 ...
    - 2 * gamma_cov (0, model,1) ) ; 
g=9.81;
delta_t=model.dt_a;
coef_homogen_a = delta_t * bsxfun(@times, ...
    (g/f_coriolis(model.coriolis)).^2, ...
    coef_homogen_a ); 
homogen_a = coef_homogen_a * permute( eye(2),[3 1 2]);% 1 x d x d

% - Inhomogeneous part of a
a= bsxfun(@plus, -a,  homogen_a) ; % nx x d x d
% Symmetric positive

%% Trace
tra = a(:,1,1)+a(:,2,2); % > 0
% (coef_homogen_a - tra ) is a minorant of the smallest eigenvalue

%% Discriminant
delta = (a(:,1,1)-a(:,2,2)).^2 + 4 * a(:,1,2).^2;
delta=sqrt(delta);

%% Smallest eigenvalue
lambda_min = coef_homogen_a - (tra+delta)/2;

%% Biggest eigenvalue
lambda_max = coef_homogen_a - (tra-delta)/2;


end