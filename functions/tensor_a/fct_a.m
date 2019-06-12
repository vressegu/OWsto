function a_res = fct_a(X_i_input,model_input)
% Compute the quadratic covariation of the geostrophic velocity
% a(x,x) = dt  (g/f)^2 *
%       ( nugget / (l_pixel/2)^2    - 2 \gamma'(0) Id
%           - 2 J beta  ( \gamma(x^{obs}) ) ^{-1} \beta^t J^t  )

%% Check if a is already computed
persistent model X_i a

if  ( isempty(model) && isempty(X_i) && isempty( a))
    % if ~ ( exist('model','var') && exist('X_i','var') && exist( 'a','var'))
    update_pers = true;
else
    update_pers = ~ ( comp_struct(model_input,model) && ...
        isequal(X_i, X_i_input)  );
end

if ~ update_pers
    a_res = a;
    return
end
model=model_input;clear model_input;
X_i=X_i_input;clear X_i_input;

%% Inhomogeneous part of a
% keyboard;
[n_x,d]=size(X_i);
n_x_m=size(model.obs.x,1);
threshold=8e7;
big_data = n_x_m*n_x>threshold;
if big_data
    inhomogen_a = nan(d,d,n_x);
    n_block=floor(threshold/n_x_m);
    for k=1:ceil(n_x*n_x_m/threshold)
        inhomogen_a_temp = fct_a_V(X_i((k-1)*n_block+1:min(k*n_block,end),:), model);% n_x_m x d x n_x
        inhomogen_a(:,:,(k-1)*n_block+1:min(k*n_block,end)) = ...
            - multiprod( multitrans(inhomogen_a_temp), inhomogen_a_temp);% d x d x n_x
    end
    clear inhomogen_a_temp
%     for k=1:n_x
%         inhomogen_a_temp = fct_a_V(X_i(k,:), model);% n_x_m x d x 1
%         inhomogen_a(:,:,k) = - inhomogen_a_temp' * inhomogen_a_temp;% d x d x 1
%     end

%     inhomogen_a_temp = fct_a_V(X_i, model);% n_x_m x d x n_x
%     inhomogen_a = nan(d,d,n_x);
%     for k=1:n_x
%         inhomogen_a(:,:,k) = - inhomogen_a_temp(:,:,k)' * inhomogen_a_temp(:,:,k);% d x d x 1
%     end
%     clear inhomogen_a_temp
else
    inhomogen_a = fct_a_V(X_i, model);% n_x_m x d x n_x
    inhomogen_a = - multiprod( multitrans(inhomogen_a), inhomogen_a);% d x d x n_x
end


%% Homogeneous part of a
l_pixel=model.l_pixel_nugget;
nugget = exp(model.kriging.lognoisevariance);
homogen_a = ( nugget / (l_pixel/2)^2 ...
    - 2 * gamma_cov (0, model,1) ) * eye(2); % d x d

%% Gather the two parts
a_res= bsxfun(@plus, homogen_a,  inhomogen_a) ; % d x d x n_x
clear inhomogen_a homogen_a;
g=model.physical_constant.g;
delta_t=model.dt_a;
a_res =  bsxfun(@times, ...
    delta_t * (g/f_coriolis(model.coriolis)).^2, ...
    a_res ); % d x d x n_x
a_res = permute(a_res,[3 1 2]);% nx x d x d

%% Save as persistent variable
a = a_res;

end