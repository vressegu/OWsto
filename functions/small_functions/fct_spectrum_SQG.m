function [w_SQG,f_SQG,spectrum_SQG,spectrum] = fct_spectrum_SQG(model,ft,fft_T)
% Compute the spectrum of a function
%

ft_w = ft;

ft=abs(ft).^2;
ft=sum(ft,3);

MX=model.grid.MX;
dX=model.grid.dX;

if any(size(ft)~=MX)
    error('wrong size');
end
% [Mx,My]=size(ft);
% MX=[Mx My];
% M=floor(sqrt(Mx*My));

% M=min(MX/2);
M=min(MX);
PX=MX/2;

persistent idxref MXref dXref kidx kxref kyref kref

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    %% Define wave number
    if any( mod(MX,2)~=0)
        error('the number of grid points by axis need to be even');
    end
    ft(PX(1),:)=0;
    ft(:,PX(2))=0;
    kx=1/(MX(1)*dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(MX(2)*dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    %     kx=1/Mx*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
    %     ky=1/My*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
%     kxref=kx;kyref=ky;
    [kx,ky]=ndgrid(kx,ky);
    kxref=kx;
    kyref=ky;
    k=sqrt(kx.^2+ky.^2);
    %     k=log(kx.^2+ky.^2)/log(10)/2; clear kx ky
    %     % k=sqrt(kxref.^2+kyref.^2); clear kxref kyref
    k=k(:);
    kref=k;
    %     % kref=repmat(kref,[1 M]);
    
    %% Order it
    kidx=linspace(0, ...
        1/max( dX(1), dX(2) )/2, ...
        M);
%     kidx=linspace(1/max(MX(1)*dX(1), MX(2)*dX(2)), ...
%         1/max( dX(1), dX(2) )/2, ...
%         M);
%     %     kidx=linspace(-max(log(MX(1)*dX(1)), ...
%     %                        log(MX(2)*dX(2))), ...
%     %                         -max( log(dX(1)), log(dX(2)) ) -log(2), ...
%     %                          M)/log(10);
%     % %     kidx=linspace(-max(log(MX(1)*dX(1)), ...
%     % %                        log(MX(2)*dX(2))), ...
%     % %         log(1/(dX(1))^2+1/(dX(2))^2)/2 -log(2), ...
%     % %         M)/log(10);
%     % % %     kidx=linspace(-2*max(log(model.grid.MX(1)*model.grid.dX(1)), ...
%     % % %         log(model.grid.MX(2)*model.grid.dX(2))), ...
%     % % %         log(1/(model.grid.dX(1))^2+1/(model.grid.dX(2))^2) -2*log(2), ...
%     % % %         min(MX/8));
%     % %     % kidx=linspace(-max(log(Mx),log(My)),-log(2),M+1);
%     % %     % kidx=(0:(M-1))/(2*M);
%     % %     % kidx=[kidx 1/2+eps];
%     % %     % kidx=repmat(kidx,[size(kref,1) 1]);
    idx = sparse( bsxfun(@le,kidx(:,1:end-1), kref ) );
    idx = idx & sparse( bsxfun(@lt,kref, kidx(:,2:end)  ) );
    idxref=idx;
end

%% Spectrum
spectrum_w = idxref' * ft(:);
if nargin>2
    spectrum_w2 = idxref' * ft2(:);
end

%% Reference

% % iiref=10;
% iiref=11;
iiref=8:15;
% k0 = 10.^(-5);
k0 = 10.^(-5.5);
k_inf = 10.^(-4) / 4;
% k_inf = 10.^(-4) / 2;

reference_spectrum = ((kidx(1:end-1)).^2 ).^(-5/6) ;
% reference_spectrum = (k0^2 + (kidx(1:end-1)).^2 ).^(-5/6) ;

reference_spectrum = reference_spectrum * mean(spectrum_w(iiref)'./reference_spectrum(iiref));
% reference_spectrum = reference_spectrum * spectrum_w(iiref) / reference_spectrum(iiref);

%% Spectre of random small-scale velocity
Gamma_SQG = ones(1,length(kidx-1));

% Pass-band filter
idx1 = (kidx(1:end-1) <= k0);
idx3 = (kidx(1:end-1) > k_inf);
idx2 = ~ (idx1 | idx3);
unit_approx = fct_unity_approx_(sum(idx2));
% figure;plot(unit_approx)
% % keyboard;

Gamma_SQG(idx1 | idx3)=0;
Gamma_SQG(idx2) = Gamma_SQG(idx2) .* unit_approx;

f_SQG = sqrt( Gamma_SQG ./ (2*pi*kidx(1:end-1)).^3 );
d_kappa = kidx(2)-kidx(1);
f_SQG = f_SQG * 1/ sqrt(prod(MX.*dX)*d_kappa); % Due to discretisation

filter = interp1(kidx(1:end-1),f_SQG,kref);

[fft_w_SQG,f_SQG] = fct_SQG(fft_T,model);
f_w_SQG(:,:,1) = filter .* fft_w_SQG(:,:,1);
f_w_SQG(:,:,2) = filter .* fft_w_SQG(:,:,2);
f_SQG = filter .* f_SQG;
n_f_W_SQG = sum(abs(f_w_SQG(:)).^2);
w_filter(:,:,1) = filter .* ft(:,:,1);
w_filter(:,:,2) = filter .* ft(:,:,2);
n_w_filter = sum(abs(w_filter(:)).^2);
model.advection.SQG.coef = sqrt(n_w_filter/n_f_W_SQG);
f_SQG=model.advection.SQG.coef*f_SQG;
clear f_w_SQG w_filter

% Cleaning
f_SQG(kref<=k0)=0;
f_SQG(kref>k_inf)=0;
f_SQG=reshape(f_SQG,MX);
% Antialiasing
marge=3;
f_SQG(PX(1)+(-marge:marge),:)=0;
f_SQG(:,PX(2)+(-marge:marge))=0;

% Influence of the complex brownian variance
f_SQG = 1/sqrt(prod(MX))*f_SQG;

w_SQG(:,:,1)= - 2*1i*pi* ( - kyref ) .* f_SQG;
w_SQG(:,:,2)= - 2*1i*pi* ( + kxref ) .* f_SQG;

ft_SQG=abs(w_SQG).^2;
ft_SQG=sum(ft_SQG,3);
spectrum_SQG = idxref' * ft_SQG(:);

% Calcul of energy
% One has to divid by prod(model.grid.MX) because of the form of Parseval
% theorem for TFD
nrj_w = 1/prod(model.grid.MX) * sum(spectrum_w)

var_w = 1/prod(model.grid.MX) * nrj_w

% Influence of the complex brownian variance
warning('necessary?')
spectrum_SQG = prod(MX)*spectrum_SQG;

%% Plot

figure; pcolor(ft_SQG);shading flat
title('Square modulus of Fourier transorfm of w_SQG');

figure;
subplot(2,1,1)
% % iiref2=11;
% % line3= (-5/3*(log(kidx(1:end-1)) - log(kidx(iiref2))) + log(spectrum_w(iiref2)))/log(10);
% line3= -5/3*log(kidx(1:end-1)) ;
% line3 = line3 + (log(mean(spectrum_w(iiref)))-log(mean(exp(line3(iiref)))));
% line3=line3/log(10);
hold on;
% plot(log(kidx(1:end-1))/log(10),line3,'r')
plot(log(kidx(1:end-1))/log(10),log(reference_spectrum)/log(10),'r')
plot(log(kidx(1:end-1))/log(10),log(spectrum_w)/log(10))
if nargin>2
    plot(log(kidx(1:end-1))/log(10),log(spectrum_w2)/log(10),'c')
end
% plot(log(kidx(1:end-1))/log(10),log(Gamma_SQG)/log(10),'g')
plot(log(kidx(1:end-1))/log(10),log(spectrum_SQG)/log(10),'k')
% plot(log(kidx(1:end-1))/log(10),2*log(f_SQG)/log(10),'g')
hold off
ax=axis;
ax(3)=2;
axis(ax);
title('Spectrum in log scale')

%% Plot
subplot(2,1,2)
% % figure;
% % iiref=10;
% % % % k0 = exp(-6);
% % % k0 = 10.^(-5);
% % % iiref=30;
% line3 = -5/3*(log(kidx(1:end-1)) - log(kidx(iiref2))) + log(spectrum_w(iiref2));
% % reference_spectrum = (k0^2 + (kidx(1:end-1)).^2 ).^(-5/6) ;
% % reference_spectrum = reference_spectrum * spectrum_w(iiref) / reference_spectrum(iiref);
% % % line4= -5/6*(log(k0^2 + (kidx(1:end-1)).^2 ) - 2*log(kidx(iiref))) + log(spectrum_w(iiref));
% % % line4= -5/6*log( (kidx(1:end-1)/kidx(iiref)).^2 +1 ) + log(spectrum_w(iiref));
hold on;
% plot(kidx(1:end-1),exp(line3),'r')
plot(kidx(1:end-1),spectrum_w)
if nargin>2
    plot(kidx(1:end-1),spectrum_w2,'c')
end
% plot(kidx(1:end-1),exp(line4),'b')
plot(kidx(1:end-1),reference_spectrum,'r')
% plot(kidx(1:end-1),Gamma_SQG,'g')
plot(kidx(1:end-1),spectrum_SQG,'k')
% plot(kidx(1:end-1),f_SQG.^2,'g')
hold off
ax=axis;
ax(4)=2e9;
axis(ax);
title('Spectrum in real scale')


end

function t = fct_unity_approx_(N_t)
% XP must be a 2 x n matrix
% the result is a vector of size n
%

slop_size_ratio=(N_t/10);
sslop=ceil(N_t/slop_size_ratio);
sslop=8;

t=ones(1,N_t);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

end
