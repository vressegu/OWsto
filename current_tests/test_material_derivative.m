%% Test material derivative

init

model.type_data='test';
init_model
N_t = 20;
n = 128;
% n = 256;
model.grid.dX =[1/n 1/n];
model.grid.MX=[n n];
x=model.grid.dX(1)*( 0:(model.grid.MX(1)-1));
y=model.grid.dX(2)*( 0:(model.grid.MX(2)-1));
model.grid.x=x;
model.grid.y=y;
[x,y]=ndgrid(x,y);
model.advection.coef_diff = 0;
model.advection.HV.val = 0;
model.advection.HV.order = 0;
model.advection.meth_anti_alias = 'none';
model = init_grid_k (model);

% w0=w(:,:,:,1);
% w0(:,:,1) = cos( 2*pi * x )*cos( 2*pi * y );
% w0(:,:,2) = 0;
w0(:,:,2) = sin( 2*pi * x );
w0(:,:,1) = 1;
std_w = sqrt(mean(w0(:).^2));

% CFL 
dt = min(model.grid.dX)/std_w;
% dt = 0.1 * min(model.grid.dX)/std_w;
model.advection.dt_adv = dt;
% model.advection.dt_adv =1;
time = dt*(0:(N_t-1));
T_w = 10*dt;
% modulation_w = cos(2*pi/T_w * time );
modulation_w = ones(size(time));
w = bsxfun( @times, w0 , permute( modulation_w,[1 3 4 2] ));

b0 =cos(2*pi*x);
fft_b = nan( [model.grid.MX 1 N_t]);
fft_b(:,:,:,1) = fft2(b0); 
for t=1:(N_t-1)
%     fct_plot_velocity(model,w(:,:,:,t));drawnow;
    fft_b(:,:,:,t+1) = ...
        RK4_fft_advection(model, fft_b(:,:,:,t), w(:,:,:,t));
    if mod(t,1)==0
        fct_plot_f(model,real(ifft2(fft_b(:,:,:,t+1))));drawnow;
    end
    %     pause;
end
b = real(ifft2(fft_b));


D_b = material_deriv_mat(b,w,model.advection.dt_adv,model.grid.dX);

std_ini = sqrt(mean(b0(:).^2));
erreur = abs(D_b) / std_ini;
mean(erreur(:))


% d_t_f = deriv_fft_advection(model, fft_f, w);
