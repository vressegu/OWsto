function filter = fct_design_filter_aniso(model)
% Design a Gaussian filter which enable to recover the initial smooth field
%


% Get parameters
MX=model.grid.MX;
PX=MX/2;
dX=model.grid.dX;
if any( mod(MX,2)~=0)
    error('the number of grid points by axis need to be even');
end

% persistent idxref MXref dXref kidx
% 
% if ~exist('MXref','var') ||  isempty (MXref) || any(MXref ~= MX) ...
%         || any(dXref ~= dX)
%     MXref=MX;
%     dXref=dX;
    
    % Remove aliasing
%     ft(PX(1),:)=0;
%     ft(:,PX(2))=0;
    
    %% Wave vector
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    kx=2*pi/model.grid.dX(1)*kx;
    ky=2*pi/model.grid.dX(2)*ky;
    [kx,ky]=ndgrid(kx,ky);
    k_v(:,:,1)=kx;
    k_v(:,:,2)=ky;
%     k=sqrt(kx.^2+ky.^2);
    k_v(PX(1)+1,:,:)=inf;
    k_v(:,PX(2)+1,:)=inf;
% %     k=k(:);
% end

L2_mat = model.spectrum_theo_aniso0.coef2 - model.spectrum_theo_aniso.coef2;
CST = model.spectrum_theo_aniso0.coef1 / model.spectrum_theo_aniso.coef1;
 
k2_mat = bsxfun (@times , k_v, permute(k_v , [1 2 4 3]));
temp = bsxfun (@times , permute(L2_mat,[3 4 1 2]) , k2_mat );
temp = sum(sum(temp,4),3);

filter = CST * exp( -1/2 * temp );

filter = sqrt(filter);
filter(PX(1)+1,:)=0;
filter(:,PX(2)+1)=0;
