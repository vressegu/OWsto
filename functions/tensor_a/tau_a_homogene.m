function T = tau_a_homogene(T, model)
% Compute \tau(T) = sum_{i,j} d_{x_i} ( a_{ij}(XP) /2  d_{x_j} T(XP) )
%  = div ( a/2 grad(T) )
% ( Laplace - Beltrami operator )
% T should be a sampled scalar function on a grid of n x m points
% size(T) should be n*m
%

%% Parameters
% [n,m]=size(T);
dX=model.grid.dX;
MX=model.grid.MX;
BOX=model.kriging.BOX;
% n = round((BOX(2,1)-BOX(1,1))/dX(1))+1;
% m = round((BOX(2,2)-BOX(1,2))/dX(2))+1;

%% Compute the gradient of T
T=reshape(T,MX);
[dTx,dTy]=gradient(T',model.grid.dX(1),model.grid.dX(2));
T(:,:,1)=dTx';
T(:,:,2)=dTy';% n m 2
clear dTx dTy

%% Compute the local variance-covariance matrix a

% Grid XP 
x= dX(1)*(0:MX(1)-1)+BOX(1,1);
y= dX(2)*(0:MX(2)-1)+BOX(1,2);
% x= dX(1)*(1:MX(1));
% y= dX(2)*(1:MX(2));
% x= linspace(BOX(1,1),BOX(2,1),MX(1));
% y = linspace(BOX(1,2),BOX(2,2),MX(2));
[x,y]=ndgrid(x,y);
XP=stk_dataframe ([x(:) y(:)]);

% Local variance-covariance matrix
a = 1/2 * fct_a_prior(XP,model); % nx d d
a=reshape(a,[MX 2 2]);% n m d d
clear XP

%% Application of the matrix
T = bsxfun(@times, a, T);
clear a
T = sum(T,3);% n m 1 d

%% Compute the divergence
T = divergence(x',y',T(:,:,:,1)',T(:,:,:,2)')';% n m
% T = divergence(y',x',T(:,:,:,1)',T(:,:,:,2)')';% n m
%%
% figure
% pcolor(T)
% keyboard;
% %%
T=T(:);

