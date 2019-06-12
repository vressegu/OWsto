function inv_tau_f = inv_J_tau_homogene_fft(model,f)
% Inverse the operator (f J - 1/(\rho) \tau ) with the RHS f
%

%% Grid
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=model.grid.MX/2;

kx=1/(model.grid.MX(1)*model.grid.dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
[kx,ky]=ndgrid(kx,ky);
k2 = model.a0/2 * (2*pi)^2*(kx.^2+ky.^2) ;
% k2 = a0/(2*model.coriolis.f0) * (2*pi)^2*(kx.^2+ky.^2) ;
clear kx ky
%% Pretreatement
f=f(:,:,1)+1i*f(:,:,2);
fftf=fft2(f);

%% Inversion in Fourier space
inv_tau_f = fftf ./ ( 1i*model.coriolis.f0 + k2 );

% Prevent aliasing
inv_tau_f(PX(1),:) = 0;
inv_tau_f(:,PX(2)) = 0;

% Force homogene Dirichlet condition
inv_tau_f(1,1) = 0;

%% Post treatement
inv_tau_f=ifft2(inv_tau_f);
inv_tau_f(:,:,2)=imag(inv_tau_f);
inv_tau_f(:,:,1)=real(inv_tau_f(:,:,1));

