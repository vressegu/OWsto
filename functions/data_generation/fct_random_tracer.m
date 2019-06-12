function [sigma_on_sq_dt,f_sigma,a0_on_dt,spectrum_sigma,spectrum] = fct_random_tracer(model,ft,ft2)
% Compute the spectrum of a function
%

if nargin >1
    ft=mean(abs(ft).^2,4);
    ft=sum(ft,3);
    if any(size(ft)~=MX)
        error('wrong size');
    end
end

MX=model.grid.MX;
dX=model.grid.dX;

if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;
ft(PX(1),:)=0;
ft(:,PX(2))=0;

persistent idxref MXref dXref kidx kxref kyref kref

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)

    MXref=MX;
    dXref=dX;
    %% Define wave number
    
    kxref=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    kyref=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    kxref=2*pi/model.grid.dX(1)*kxref;
    kyref=2*pi/model.grid.dX(2)*kyref;
    [kxref,kyref]=ndgrid(kxref,kyref);
    kref=sqrt(kxref.^2+kyref.^2);
%     kref(PX(1)+1,:)=max(kref(:));
%     kref(:,PX(2)+1)=max(kref(:));
    kref(PX(1)+1,:)=inf;
    kref(:,PX(2)+1)=inf;
    kref=kref(:);
    
    
    %% Order it
    
    M_kappa=min(MX);
    P_kappa= M_kappa/2;
    d_kappa = max(1./dX);
    kidx=1/(M_kappa)* (0:(P_kappa-1)) ;
    kidx=2*pi*d_kappa*kidx;
    idx = sparse( bsxfun(@lt,kidx(1:end-1), kref ) );
    idx = idx & sparse( bsxfun(@le,kref, kidx(2:end)*(1+10*eps) ) );
    idxref=idx;
end

%% Spectrum
if nargin>1
    spectrum_w = idxref' * ft(:);
end

%% Reference
%     k_inf = kidx(min(PX));
%     k0 = 1/2*k_inf;
alpha = model.slope_spectrum;
k_inf = model.k_min;
k0 = model.k_max;

if k_inf > kidx(min(PX))
    k_inf = kidx(min(PX));
end

reference_spectrum = kidx(2:end) .^ ( - alpha) ;
if nargin > 1
    idx_not_inf=~(isinf(log10(spectrum_w))| spectrum_w<1e-4*max(spectrum_w) | isinf(kidx(2:end)'));
    reference_spectrum = reference_spectrum * 10 .^( ...
        mean(  log10(spectrum_w(idx_not_inf)')  - log10( reference_spectrum(idx_not_inf)) ) );
end
%% Spectre of random small-scale velocity
Gamma_sigma = reference_spectrum;

% Pass-band filter
idx1 = (kidx(2:end) <= k0);
idx3 = (kidx(2:end) > k_inf);
idx2 = ~ (idx1 | idx3);
unit_approx = fct_unity_approx_(sum(idx2));

Gamma_sigma(idx1 | idx3)=0;
Gamma_sigma(idx2) = Gamma_sigma(idx2) .* unit_approx;

f_sigma = sqrt( Gamma_sigma ./ ( kidx(2:end) ) );
% f_sigma = sqrt( Gamma_sigma ./ ( kidx(2:end) .^3 ) );
d_kappa = 1/(2*pi) * ( kidx(2)-kidx(1) );
f_sigma = f_sigma * 1/ sqrt(prod(MX.*dX)*d_kappa); % Due to discretisation
f_sigma = interp1(kidx,[0 f_sigma],kref);

% Cleaning
f_sigma(kref<=k0)=0;
f_sigma(kref>k_inf)=0;
f_sigma=reshape(f_sigma,MX);
% Antialiasing
f_sigma(PX(1)+1,:)=0;
f_sigma(:,PX(2)+1)=0;

% figure; pcolor(abs(f_sigma));shading flat
% keyboard;

% Influence of the complex brownian variance
f_sigma = 1/sqrt(prod(MX))*f_sigma;

% sigma_on_sq_dt(:,:,1)= 1i * ( - kyref ) .* f_sigma;
% sigma_on_sq_dt(:,:,2)= 1i * ( + kxref ) .* f_sigma;
sigma_on_sq_dt = f_sigma;

ft_sigma=abs(sigma_on_sq_dt).^2;
ft_sigma=sum(ft_sigma,3);
spectrum_sigma = idxref' * ft_sigma(:);

% Calcul of energy
% One has to divid by prod(model.grid.MX) because of the form of Parseval
% theorem for Discrete Fourier Transform
a0_on_dt = 1/prod(model.grid.MX) * sum(spectrum_sigma);

% Influence of the complex brownian variance
spectrum_sigma = prod(MX)*spectrum_sigma;

%% Plot

figure; pcolor(ft_sigma);shading flat
title('Square modulus of Fourier transorfm of sigma_on_sq_dt');

%%
figure4=figure(4);
widthtemp = 12;
heighttemp = 6;
set(figure4,'Units','inches', ...
    'Position',[0 0 widthtemp heighttemp], ...
    'PaperPositionMode','auto');

loglog(kidx(2:end),reference_spectrum,'k')
hold on;
if nargin > 1
    loglog(kidx(2:end),spectrum_w)
end
loglog(kidx(2:end),spectrum_sigma,'r')

hold off
ax=axis;
if nargin > 1
    ax(4)=max([spectrum_w; reference_spectrum' ; spectrum_sigma]);
    ax(4)=ax(4)*2;
    ax(3)=(kidx(2)/kidx(end))*min([max(spectrum_w); max(reference_spectrum); max(spectrum_sigma)]);
else
    ax(4)=max([reference_spectrum' ; spectrum_sigma]);
    ax(4)=ax(4)*2;
    ax(3)=(kidx(2)/kidx(end))*min([ max(reference_spectrum); max(spectrum_sigma)]);
end
ax(3) = min( [ax(3) min(reference_spectrum)]);
ax(1:2)=kidx([2 end]);
if ax(4)>0
    axis(ax)
end

set(gca,'XGrid','on','XTickMode','manual');
width = 9;
height = 3;
set(figure4,'Units','inches', ...
    'Position',[0 0 width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')
taille_police = 12;
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('$|\hat{f}(\kappa)|^2$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('$\kappa$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
title(['Spectrum of $w$ and $\sigma dB_t$'],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
        


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
