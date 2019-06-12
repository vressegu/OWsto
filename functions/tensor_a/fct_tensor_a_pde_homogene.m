function a_res = fct_tensor_a_pde_homogene(XP,model)
% Compute the quadratic covariation of the geostrophic velocity for PDE
% toolbox
% a(x,x) = dt  (g/f)^2 * 
%       ( nugget / (l_pixel/2)^2    - 2 \gamma'(0) Id 
%           - 2 J beta  ( \gamma(x^{obs}) ) ^{-1} \beta^t J^t  )

% a = fct_a(XP,model);% d x d x n_x
% a_res(1,:)=permute(a(1,1,:),[1 3 2]);
% a_res(2,:)=permute(a(1,2,:),[1 3 2]);
% a_res(3,:)=permute(a(2,2,:),[1 3 2]);

a = fct_a_prior(XP,model);% nx d x d
a_res(1,:)=a(:,1,1)';
a_res(2,:)=a(:,1,2)';
a_res(3,:)=a(:,2,2)';


