function T = tau_a_homogene_unit_approx(T, model)
% Compute \tau(T) = sum_{i,j} d_{x_i} ( a_{ij}(XP) /2  d_{x_j} T(XP) )
%  = div ( a/2 grad(T) )
% ( Laplace - Beltrami operator )
% T should be a sampled scalar function on a grid of n x m points
% size(T) should be n*m
%

%% Parameters
% [n,m]=size(T);
dX=model.dX;
BOX=model.kriging.BOX;
n = round((BOX(2,1)-BOX(1,1))/dX(1))+1;
m = round((BOX(2,2)-BOX(1,2))/dX(2))+1;

%% Compute the gradient of T
T=reshape(T,[n m]);
[dTx,dTy]=gradient(T',model.dX(1),model.dX(2));
T(:,:,1)=dTx';
T(:,:,2)=dTy';% n m 2
clear dTx dTy

%% Compute the local variance-covariance matrix a

% Grid XP 
x= linspace(BOX(1,1),BOX(2,1),n);
y = linspace(BOX(1,2),BOX(2,2),m);
[x,y]=ndgrid(x,y);
XP=stk_dataframe ([x(:) y(:)]);

% Local variance-covariance matrix
a = 1/2* bsxfun(@times, fct_unity_approx2(XP,BOX)' , ...
       fct_a_prior(XP,model)); % nx d d
a=reshape(a,[n m 2 2]);% n m d d
clear XP

%% Application of the matrix
T = bsxfun(@times, a, T);
clear a
T = sum(T,3);% n m 1 d

%% Compute the divergence
T = divergence(x',y',T(:,:,:,1)',T(:,:,:,2)')';% n m
% T = divergence(y',x',T(:,:,:,1)',T(:,:,:,2)')';% n m
T=T(:);

