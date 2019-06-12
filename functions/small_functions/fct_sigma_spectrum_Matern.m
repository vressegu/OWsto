function [sigma_on_sq_dt,f_sigma,a0_on_dt,spectrum_sigma] = ...
    fct_sigma_spectrum_Matern(model)
% - sigma_on_sq_dt is the Fourier transform of the kernel \tilde sigma up to a multiplicative
% constant
% - f_sigma is the Fourier transform of the associted streamfunction
% - a0_on_dt measures the total energy of the field and is generally used to
% set the muliplicative constant
% - spectrum_sigma is the spectrum
%

% Get parameters
MX=model.grid.MX;
dX=model.grid.dX;
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=MX/2;

%% Wave vector
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
kx=2*pi/model.grid.dX(1)*kx;
ky=2*pi/model.grid.dX(2)*ky;
d_kx = kx(2)-kx(1);
d_ky = ky(2)-ky(1);
[kx,ky]=ndgrid(kx,ky);
k=sqrt(kx.^2+ky.^2);
k(PX(1)+1,:)=inf;
k(:,PX(2)+1)=inf;
k=k(:);

%% Wave number
M_kappa=min(MX);
P_kappa= M_kappa/2;
d_kappa = max(1./dX);
kappa=1/(M_kappa)* (0:(P_kappa-1)) ;
kappa=2*pi*d_kappa*kappa;

%% Masks associated with the rings of iso wave number
d_kappa = kappa(2) - kappa(1);
idx = sparse( bsxfun(@le,kappa, k ) );
idx = idx & sparse( bsxfun(@lt,k, [ kappa(2:end) kappa(end)+d_kappa ] ) );

%% 1D Spectrum

% Largest wave number
% k_M = kappa(min(PX));
k_M = model.sigma.kappamax_on_resolution * kappa(min(PX));

% Smallest wave number
k_m = model.sigma.kappamin_on_kappamax * k_M;

% warning('division by 2 pi');
% k_M = k_M *(2*pi);

% Spectrum with slope alpha
% reference_spectrum = kappa .^ model.sigma.slope_sigma ;
reference_spectrum = ( 1/k_m * kappa(2:end) ) .^ (model.sigma.slope_sigma-3) ;
reference_spectrum( kappa(2:end) < k_m ) =1;
reference_spectrum = reference_spectrum .* kappa(2:end) .^ 3 ;

%% Spectre of random small-scale velocity
% f_sigma = reference_spectrum;
% 
% % Parseval ( * prod(model.grid.MX) ) 
% % and integrated over the space ( * prod(model.grid.MX) )
% f_sigma = prod(model.grid.MX)^2 * f_sigma;
% 
% % Band-pass filter
% idx1 = (kappa(2:end) <= k_m);
% idx3 = (kappa(2:end) > k_M);
% idx2 = ~ (idx1 | idx3);
% unit_approx = fct_unity_approx_(sum(idx2));
% f_sigma(idx1 | idx3)=0;
% f_sigma(idx2) = f_sigma(idx2) .* unit_approx;
% 
% % Division by k^2 to get the spectrum of the streamfunction
% f_sigma = f_sigma ./ ( kappa(2:end) .^2 );
% 
% % From omnidirectional spectrum to Fourier tranform square modulus
% % Division by k in dimension 2
% f_sigma = f_sigma ./ ( kappa(2:end) );
% 
% % Influence of discretisation
% f_sigma = f_sigma / ( prod(MX.*dX) /(2*pi) ) ;
% 
% % From square modulus to modulus
% f_sigma = sqrt( f_sigma );
% 
% % From 1D function to 2D function
% f_sigma = interp1(kappa,[0 f_sigma],k);

V0 =model.sigma.V0;
beta = - model.sigma.slope_sigma;
A2 = 2*pi * (beta+1) * V0 / (k_m^2);
A2 = 1/prod(model.grid.dX) * A2;
% The factor 1/prod(dX) is to go from continuous to discrete Fourier
% transform
f_sigma = sqrt(A2) * (1 + (1/k_m * k).^2 ).^(-(beta+3)/4);

% Cleaning
% f_sigma(k<=k_m)=0;
f_sigma(k>k_M)=0;
f_sigma=reshape(f_sigma,MX);

% Antialiasing
f_sigma(PX(1)+1,:)=0;
f_sigma(:,PX(2)+1)=0;

% % Influence of the complex brownian variance
% f_sigma = 1/sqrt(prod(MX))*f_sigma;

% Orthogonal gradient (from streamfunction to velocity)
sigma_on_sq_dt(:,:,1)= 1i * ( - ky ) .* f_sigma;
sigma_on_sq_dt(:,:,2)= 1i * ( + kx ) .* f_sigma;

% Compute the real spectrum of sigma dBt
ft_sigma=abs(sigma_on_sq_dt).^2;
ft_sigma=sum(ft_sigma,3);
spectrum_sigma = idx' * ft_sigma(:);

% % Influence of the complex brownian variance
% spectrum_sigma = prod(MX)*spectrum_sigma;

% One has to divid by prod(model.grid.MX) because of the form of Parseval
% theorem for discrete Fourier transform
spectrum_sigma = 1/prod(model.grid.MX) * spectrum_sigma;

% Calcul of energy
% a0_on_dt = 1/prod(model.grid.MX) * sum(spectrum_sigma);
a0_on_dt = 1/2 * sum(spectrum_sigma);
% % the factor d=2 is for the dimension d of the space R^d
% % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt

% Influence of the complex brownian variance
spectrum_sigma = prod(MX)*spectrum_sigma;

% Division by prod(model.grid.MX) in order to the integration 
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum_sigma = 1/prod(model.grid.MX) * spectrum_sigma;

% Test the mean energy
% mean_nrj = sum(spectrum_sigma)

% Integration over the ring kappa
spectrum_sigma = spectrum_sigma * (d_kx * d_ky)/d_kappa;
warning('modif integration: should modif energy budget');

% Division by the wave number step
spectrum_sigma = spectrum_sigma / d_kappa;

% Compute the offset to superimpose the slop - beta and the spectrum of w
% reference_spectrum = 2*pi * A2 / (prod(model.grid.MX) * d_kappa ) ...
reference_spectrum = 2*pi * A2 / (prod(model.grid.MX) * d_kappa ) ...
    * reference_spectrum;
% warning('debug reference spectrum');
% idx_not_inf=~(isinf(log10(spectrum_sigma(2:end)))| ...
%     spectrum_sigma(2:end)<1e-4*max(spectrum_sigma(2:end)) | isinf(kappa(2:end)'));
% idx_not_inf= [ false; idx_not_inf ];
% reference_spectrum = reference_spectrum * 10 .^( ...
%     mean( log10(spectrum_sigma(idx_not_inf)') ...
%     - log10( reference_spectrum(idx_not_inf(2:end))) ) );

%% Plot spectrum
taille_police = 12;
X0 = [10 20];

figure10=figure(10);
widthtemp = 12;
heighttemp = 6;
set(figure10,'Units','inches', ...
    'Position',[X0 widthtemp heighttemp], ...
    'PaperPositionMode','auto');
loglog(kappa(2:end),reference_spectrum,'k')
hold on;
% loglog(kappa(2:end),spectrum_w(2:end))
% if nargin>2
%     loglog(kappa(2:end),spectrum_w2(2:end),'c')
% end
loglog(kappa(2:end),spectrum_sigma(2:end),'r')
hold off
% ax=axis;
% ax(4)=max([spectrum_w; reference_spectrum' ; spectrum_sigma]);
% ax(4)=ax(4)*2;
% ax(3)=(kappa(2)/kappa(end))*min(max(spectrum_sigma));
% % ax(3)=(kappa(2)/kappa(end))*min([max(spectrum_w); ...
% %     max(reference_spectrum); max(spectrum_sigma)]);
% % ax(3) = min( [ax(3) min(reference_spectrum)]);
% ax(1:2)=kappa([2 end]);
% if ax(4)>0
%     axis(ax)
% end
set(gca,'XGrid','on','XTickMode','manual');
width = 9;
height = 3;
set(figure10,'Units','inches', ...
    'Position',[X0 width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')
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
xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
title('Initial spectrum of $\sigma dB_t$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')

%% Save plot
drawnow

% folder_simu = model.folder.folder_simu;
% eval( ['print -depsc ' folder_simu '/spectrum_sigma_dB_t.eps']);

end

function t = fct_unity_approx_(N_t)
% Approximation of unity
%

sslop=8;
t=ones(1,N_t);
t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;

end
