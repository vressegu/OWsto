function [spectrum,kidx,name_plot] = fct_spectrum_1D(t,ft,color)
% Compute the spectrum of a function and superimposed a slope -5/3
%

% Color by default
if nargin < 3
    color='b';
end

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
% if any( mod(Mx,2)~=0)
%     error('the number of grid points by axis need to be even');
% end

% Remove aliasing
if even
    ft(P+1)=0;
end


%% Wave vector
% if even
%     kx = 1/Mx*[ 0:(P-1) 0 (1-P):-1] ;
% else
%     v1 =  0:(P-1);
%     v2 = -1:-1:(1-P);
%     v2 = v2(end:-1:1);
%     kx = 1/Mx*[v1 v2] ;
% end
% kx = 2*pi/dX * kx;
% [kx,ky]=ndgrid(kx,kx);
% k=sqrt(kx.^2+ky.^2);
% if even
%     k(P+1,:)=inf;
%     k(:,P+1)=inf;
% end
% k=k(:);

%% Discrete wave number
kidx=1/Mx * (0:(P-1)) ;
kidx=2*pi/dX*kidx;

%% Masks associated with the rings of iso wave number
% d_kappa = kidx(2) - kidx(1);
% idx = sparse( bsxfun(@le,kidx, k ) );
% idx = idx & sparse( bsxfun(@lt,k, [ kidx(2:end) kidx(end)+d_kappa ] ) );
% idxref=idx;
% 
% % figure;
% % for j=1:size(idxref,2)
% %     imagesc(reshape(idxref(:,j),[Mx Mx])');axis xy;axis equal
% %     pause(0.3);
% % end
% 
% idxref=double(sum(idxref,1));
% figure;plot(idxref);

%% Spectrum
spectrum = ft(1:floor(P));

% Multiplication by 2 for positive and negative frequencies
spectrum(2:end) = 2 * spectrum(2:end);

% Division by Mx because of the Parseval theorem for
% 1D discrete Fourier transform
spectrum = 1/Mx * spectrum;

% Division by Mx again in order to the integration 
% of the spectrum over the wave number yields the energy of the
% buoyancy averaged (not just integrated) over the space
spectrum = 1/Mx * spectrum;

% Division by the wave number step
d_kappa = kidx(2)-kidx(1);
spectrum = spectrum / d_kappa;

% Check energy
% sum(spectrum*d_kappa)

%% Plot
if ~isempty(color)
    name_plot = loglog(kidx(2:end) , spectrum(2:end) ,color);
    ax=axis;
    % ax(4)=max(spectrum(2:end));
    % ax(3)=min(spectrum(2:end));
    ax(1:2)=kidx([2 end]);
    if ax(4)>0
        axis(ax)
    end
end
