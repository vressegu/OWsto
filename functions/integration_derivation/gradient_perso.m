function grad_f = gradient_perso(grid, f)
% Compute the orthogonal gradient in pseudo spectral
%

%% Fourier transform
fft_f = fft2(f);

%% Grid
PX=grid.MX/2;
if ~isfield(grid,'k')
    if any( mod(grid.MX,2)~=0)
        error('the number of grid points by axis need to be even');
    end
    kx=1/(grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
        .* fct_unity_approx5_(grid.MX(1));
    ky=1/(grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
        .* fct_unity_approx5_(grid.MX(2));
    [kx,ky]=ndgrid(kx,ky);
    kx=2*pi/grid.dX(1)*kx;
    ky=2*pi/grid.dX(2)*ky;
%     k2=kx.^2+ky.^2;
%     kx(PX(1)+1,:)=0;
%     kx(:,PX(2)+1)=0;
else
    kx=grid.k.kx;
    ky=grid.k.ky;
    k2=grid.k.k2;
end

%% Gradient
fft_grad_f(:,:,1,:) =  bsxfun(@times, + 1i * kx , fft_f );
fft_grad_f(:,:,2,:) =  bsxfun(@times, + 1i * ky , fft_f);

%% Inverse Fourier transform
fft_grad_f(PX(1)+1,:,:,:)=0;
fft_grad_f(:,PX(2)+1,:,:)=0;
grad_f = real(ifft2(fft_grad_f));

end

function t = fct_unity_approx5_(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

slop_size_ratio=6;
% slop_size_ratio=ceil(6*(N_t/1000));
% N_t=1000;

t=ones(1,N_t);
P_t=N_t/2;
sslop=ceil(N_t/slop_size_ratio);
t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
% t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

t(P_t+1)=0;

end

