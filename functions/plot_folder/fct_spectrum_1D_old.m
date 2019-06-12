function [spectrum,kappa,name_plot] = fct_spectrum_1D_old(model,t,ft,color)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end

% % figure;plot(abs((ft(1,:))));
% figure;plot(abs([ft(1,:) ft(1,:)]));
% % figure;plot(real(ifft(ft(1,:))));
% % figure;imagesc(real(ifft2(ft')));axis equal;axis xy

% Square modulus of the Fourier transform
ft=abs(ft).^2;

if size(ft,3)>1
    ft = sum(ft,3);
end

% Get parameters
Mx=length(t);
even = (floor(Mx/2)==Mx/2);
P=Mx/2;
dX=t(2)-t(1);
if any(size(ft)~=[1 Mx])
    error('wrong size');
end
if any( mod(Mx,2)~=0)
    error('the number of grid points by axis need to be even');
end

% ft = ones(size(ft));

sum(ft)/(Mx^2)

% Remove aliasing
if even
    ft(P+1)=0;
end


figure;plot(ft)
nrj = sum(ft)/(Mx^2)


%% Wave vector
if even
    kx = 1/Mx*[ 0:(P-1) 0 (1-P):-1] ;
else
    v1 =  0:(P-1);
    v2 = -1:-1:(1-P);
    v2 = v2(end:-1:1);
    kx = 1/Mx*[v1 v2] ;
end
kx = 2*pi/dX * kx;
[kx,ky]=ndgrid(kx,kx);
k=sqrt(kx.^2+ky.^2);
if even
    k(P+1,:)=inf;
    k(:,P+1)=inf;
end
k=k(:);

%% Discrete wave number
kidx=1/Mx * (0:(P-1)) ;
% kidx=2*pi*kidx;
kidx=2*pi/dX*kidx;

%% Masks associated with the rings of iso wave number
d_kappa = kidx(2) - kidx(1);
idx = sparse( bsxfun(@le,kidx, k ) );
idx = idx & sparse( bsxfun(@lt,k, [ kidx(2:end) kidx(end)+d_kappa ] ) );
idxref=idx;

% figure;
% for j=1:size(idxref,2)
%     imagesc(reshape(idxref(:,j),[Mx Mx])');axis xy;axis equal
%     pause(0.3);
% end

idxref=double(sum(idxref,1));
figure;plot(idxref);

%% Spectrum
% Integration over the rings of iso wave number
% spectrum = 2*pi* kidx .* ft(1:floor(P));
spectrum = idxref .* ft(1:floor(P));
% spectrum = idxref' * ft(:);

% Division by Mx because of the Parseval theorem for
% 1D discrete Fourier transform
spectrum = 1/Mx * spectrum;

% Division by (Mx^2) again in order to the integration 
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/Mx^2 * spectrum;

sum(spectrum)/nrj
sum(idxref)/Mx^2

% Division by the wave number step
% d_kappa = 2*pi/(dX*Mx);
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

sum(spectrum*d_kappa)
keyboard;

%% Plot
% idx_not_inf=~(isinf(log10(spectrum(2:end))) ...
%     | spectrum(2:end)<1e-4*max(spectrum(2:end)) | isinf(kidx(2:end)'));
% line1= -5/3 * log10(kidx(2:end))  ;
% offset = -1 + mean(  log10(spectrum(idx_not_inf)')  - line1(idx_not_inf));
% line1 = line1 + offset;
% ref=10.^line1;
% % if isfield(model,'spectrum_theo')
% %     loglog(kidx(2:end),spectrum_theo(2:end),'r');
% %     hold on;
% % end
% % loglog(kidx(2:end),ref,'--k');
% % hold on;
name_plot = loglog(kidx(2:end) , spectrum(2:end) ,color);
ax=axis;
% ax(4)=max([spectrum(2:end); ref']);
ax(4)=max(spectrum(2:end));
min_ax= 10 ^(-5/3 * log10(kidx(2)*512/2) + offset) ;
ax(3)=min(spectrum(2:end)); 
% ax(3)=6e-2*(kidx(2)/kidx(end))*min([max(spectrum); max(ref)']);
% ax(3) = min( [ax(3) min(ref) min_ax]);
% ax(3) = max( [ax(3) ax(4)*1e-6 ]);
% ax(3) = max( [ax(3) ax(4)*1e-10 ]);
ax(1:2)=kidx(2)*[1 min(model.grid.MX)/2];
if ax(4)>0
    axis(ax)
end
