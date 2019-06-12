function [v, fft_phi] = fct_SQG2(fft_T,model)
% Compute the Fourier transform of the velocity and of the streamfunction from
% the Fourier trasform of the temperature

%% Grid
M=prod(model.grid.MX);
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=model.grid.MX/2;
kx=2*pi/(model.grid.MX(1)*model.grid.dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
ky=2*pi/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];

[kx,ky]=ndgrid(kx,ky);
k=sqrt(kx.^2+ky.^2);

%% SQG relation
fft_phi = - fft_T ./ k ;
fft_phi(1,1)=0;
warning('jets are removed from the temperature before computing SQG balance');
fft_phi(1,:)=0;
fft_phi(:,1)=0;
fft_phi(PX(1)+1,:)=0;
fft_phi(:,PX(2)+1)=0;
fft_v(:,:,1) = - 1i* ky .* fft_phi;
fft_v(:,:,2) =  1i* kx .* fft_phi;
v = ifft2(fft_v);

warning('need to filter before estimating the SQG coefficient');

%% Coefficient fitted on the kriging nugget
coef = (model.physical_constant.g/model.coriolis.f0).^2 *...
    ( model.nugget / (model.l_pixel_nugget/2)^2 ) ...
    /( 1/(2*prod(model.grid.MX)) * sum(abs(v(:).^2)) );
%     /( 1/(2*(prod(model.grid.MX))^2) * sum(abs(fft_v(:).^2)) );
coef = sqrt(coef);

fft_phi=coef*fft_phi;
v=coef*v;

% phi = ifft2(fft_phi);
% sum(abs(phi(:).^2))/ (1/prod(model.grid.MX)*sum(abs(fft_phi(:).^2)))
% keyboard;

