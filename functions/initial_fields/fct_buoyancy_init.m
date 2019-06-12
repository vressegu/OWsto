function [fft_sst,model] = fct_buoyancy_init(model,resolution)
% Create an initial buoyancy field
%

%% Grid
n=resolution;
m=n;
Lx=1e6;
dx=Lx/n;
dy=dx;
x= dx*(0:n-1);
y= dy*(0:m-1);
model.grid.origin=[0 0];
model.grid.x_ref=x;
model.grid.y_ref=y;
[x,y]=ndgrid(x,y);
model.grid.dX=[dx dy];
MX=[n m];
model.grid.MX=MX;

%% Spatial buoyancy field
switch model.type_data
    case 'Vortices'
        b_S = init_Vortices(model,x,y);
    case 'Perturbed_vortices'
        b_S = init_Perturbed_vortices(model,x,y);
    case 'Spectrum'
        b_S = init_Spectrum(model);
    case 'Back_grad'
        model.slope_b_ini = - 5/3;
        b_S = 1/4 * init_Spectrum(model);
%         b_S = init_Back_grad(model,x,y);
%         b_S = b_S + 1/10 * init_Spectrum(model);
    otherwise
        error('this type of initial condition is unknown')
end

%% Fourier transform of the buoyancy field
fft_sst = fft2(b_S);

