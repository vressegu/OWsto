function [spectrum] = fct_spectrum(model,ft)
% Compute the spectrum of a function
%

ft=abs(ft).^2;

MX=model.grid.MX;
dX=model.grid.dX;

if any(size(ft)~=MX)
    error('wrong size');
end
% [Mx,My]=size(ft);
% MX=[Mx My];
% M=floor(sqrt(Mx*My));
M=min(MX/2);
% M=min(MX/4);
persistent idxref MXref dXref kidx

if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
        || any(dXref ~= dX)
    MXref=MX;
    dXref=dX;
    %% Define wave number
    if any( mod(MX,2)~=0)
        error('the number of grid points by axis need to be even');
    end
    PX=MX/2;
    ft(PX(1),:)=0;
    ft(:,PX(2))=0;
    kx=1/(MX(1)*dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(MX(2)*dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
%     kx=1/Mx*[ 0:(PX(1)-1) 0 (1-PX(1)):-1];
%     ky=1/My*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    [kx,ky]=ndgrid(kx,ky);
%     k=log(kx.^2+ky.^2)/log(10); clear kx ky
    k=log(kx.^2+ky.^2)/log(10)/2; clear kx ky
    % /2 is for the sqrt of the norm 2
%     % k=sqrt(kx.^2+ky.^2); clear kx ky
    k=k(:);
    % k=repmat(k,[1 M]);
    
    %% Order it
    kidx=linspace(-max(log(MX(1)*dX(1)), ...
                       log(MX(2)*dX(2))), ...
                    log(1/max( dX(1), dX(2) )/2)*(1+eps), ...
                     M-1)/log(10);
    kidx = [-inf kidx];
%     kidx=linspace(-max(log(MX(1)*dX(1)), ...
%                        log(MX(2)*dX(2))), ...
%                         -max( log(dX(1)), log(dX(2)) ) -log(2), ...
%                          M)/log(10);
%     kidx=linspace(-max(log(MX(1)*dX(1)), ...
%                        log(MX(2)*dX(2))), ...
%         log(1/(dX(1))^2+1/(dX(2))^2)/2 -log(2), ...
%         M)/log(10);
% %     kidx=linspace(-2*max(log(model.grid.MX(1)*model.grid.dX(1)), ...
% %         log(model.grid.MX(2)*model.grid.dX(2))), ...
% %         log(1/(model.grid.dX(1))^2+1/(model.grid.dX(2))^2) -2*log(2), ...
% %         min(MX/8));
%     % kidx=linspace(-max(log(Mx),log(My)),-log(2),M+1);
%     % kidx=(0:(M-1))/(2*M);
%     % kidx=[kidx 1/2+eps];
%     % kidx=repmat(kidx,[size(k,1) 1]);
    idx = sparse( bsxfun(@le,kidx(:,1:end-1), k ) );
    idx = idx & sparse( bsxfun(@lt,k, kidx(:,2:end)  ) );
    idxref=idx;
end

%% Spectrum
spectrum = idxref' * ft(:);

% dkidx=kidx(2)-kidx(1);

% %% Time step choice
% dx_mesure = 1e4;
% k_mesure = -log(2*dx_mesure)
% % idx2 = sparse( bsxfun(@le,kidx(:,1:end-1), k_mesure ) );
% % idx2 = idx2 & sparse( bsxfun(@lt,k_mesure, kidx(:,2:end)  ) );
% idx2 = sparse( bsxfun(@le,k_mesure, k ) );
% idx2 = idx2 & sparse( bsxfun(@lt,k, k_mesure) );
% idx2ref=reshape(idx2,MX);
% spectrum2 = idx2ref .* ft;
% U=ifft2(spectrum2);


%% Threshold
% threshold=log(max(spectrum(1:end-1)))/log(10)-6;
% % threshold=log(max(spectrum(1:end-1)))-14;
% % 10^(-threshold)
% idx_threshold = find(log(spectrum)/log(10) > threshold);
% idx_threshold=idx_threshold(end);
% lMAX=10.^(-kidx(idx_threshold+[0 1]) )
% % 10^(-kidx(idx_threshold+1) )

%% Plot
idx_not_inf=~(isinf(log(spectrum))| spectrum<1 | isinf(kidx(1:end-1)'));
iiref=20:40;
% iiref=8:15;
% % iiref=30;
line1= -5*kidx(1:end-1) ;
line1 = line1 + mean(  log(spectrum(idx_not_inf)')/log(10) - line1(idx_not_inf));
% line1 = line1 + (log(mean(spectrum(iiref))) ...
%     -log(mean(10.^(line1(iiref)))))/log(10);
line2= -2*kidx(1:end-1) ;
line2 = line2 + mean( log(spectrum(idx_not_inf)')/log(10) - line2(idx_not_inf));
% line2 = line2 + (log(mean(spectrum(iiref))) ...
%     -log(mean(10.^(line2(iiref)))))/log(10);
% % % line1= -5*((kidx(2:end) - kidx(iiref))) + log(spectrum(iiref))/log(10);
% % % line2= -2*((kidx(2:end) - kidx(iiref))) + log(spectrum(iiref))/log(10);
plot((log(2*pi)/log(10)+ kidx(1:end-1)),line1,'k')
plot((log(2*pi)/log(10)+ kidx(1:end-1)),line2,'k')
plot((log(2*pi)/log(10)+ kidx(1:end-1)),log(spectrum)/log(10))
% keyboard;
% plot((kidx(1:end-1)),line1,'k')
% plot((kidx(1:end-1)),line2,'k')
% plot((kidx(1:end-1)),log(spectrum)/log(10))
% % plot((kidx(2:end)),line1,'k')
% % plot((kidx(2:end)),line2,'k')
% % plot((kidx(2:end)),log(spectrum)/log(10))
% % % axis equal
% % % hold on; plot(kidx(2:end),ones(1,M-1)*threshold,'r');hold off
% % % % figure; plot((kidx(2:end-1)),log(spectrum(1:end-1))/log(10))
% % % % hold on; plot(kidx(2:end-1),ones(1,M-2)*threshold,'r');hold off
% % % % keyboard;
% % % % % figure; plot(log(kidx(2:end-1)),log(spectrum(2:end)))
% % % % % figure; plot(log(kidx(1:end-1)),log(spectrum))


