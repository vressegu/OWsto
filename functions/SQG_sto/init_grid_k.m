function model = init_grid_k (model)
% Create a grid in the Fourier space
%

if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

%%  Damped Fourier grid
meth_anti_alias=model.advection.meth_anti_alias;
PX=model.grid.MX/2;
if strcmp(meth_anti_alias,'deriv_LowPass')
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
        .* fct_unity_approx5(model.grid.MX(1));
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1] ...
        .* fct_unity_approx5(model.grid.MX(2));
else
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
end
[kx,ky]=ndgrid(kx,ky);
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
k2=kx.^2+ky.^2;
k2(PX(1)+1,:)=0;
k2(:,PX(2)+1)=0;
k=sqrt(k2);
% k=sqrt(kx.^2+ky.^2);
% k(PX(1)+1,:)=0;
% k(:,PX(2)+1)=0;

% Specific operators
on_k = 1./k;
on_k ( k==0 ) = 0;
% F_k = sqrt(1+k.^4/model.k_c^4) ;
% on_F_k = 1./F_k;

%% Save
model.grid.k.kx=kx;
model.grid.k.ky=ky;
model.grid.k.k2=k2;
model.grid.k.k=k;
% model.grid.k.on_F_k=on_F_k;
model.grid.k.on_k=on_k;
clear k kx ky

%% Unstable Fourier grid
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
[kx,ky]=ndgrid(kx,ky);
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
k2=kx.^2+ky.^2;
k2(PX(1)+1,:)=0;
k2(:,PX(2)+1)=0;
% k=sqrt(k2);
% % k=sqrt(kx.^2+ky.^2);
% % k(PX(1)+1,:)=0;
% % k(:,PX(2)+1)=0;

%% Save
% model.grid.k_HV.kx=kx;
% model.grid.k_HV.ky=ky;
% model.grid.k_HV.k=k;
model.grid.k_HV.k2=k2;

