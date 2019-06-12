function sst=fct_sst_toy2(model)
% Create a simple temperature gradient
%

odg_b = 1e-3;
% % sigma = 0.002;
% sigma = 0.03;
sigma = 0.1;
% % sigma = 0.3;

% sst=1/model.grid.MX(2)*(0:(model.grid.MX(2)-1));
my2=model.grid.MX(2)/2;
y = 1/2:(my2-1/2);
y  = y / (2*my2);
% y  = y /my2;
y = [ -y(end:-1:1) y];

Ly = model.grid.dX(2) * model.grid.MX(2);
y = y * Ly;
sigma = sigma *Ly;

sst = exp( - 1/2 * y.^2 / sigma^2 );

% sst = sst - mean(sst(:));

% sst = odg_b * (2*sst-1);
% % sst=[ones(1,my2-n) zeros(1,2*n) ones(1,model.grid.MX(2)-my2-n)];
% % % sst=[ones(1,my2) zeros(1,model.grid.MX(2)-my2)];
sst=repmat(sst,[model.grid.MX(1) 1]);

% % x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% % y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% % imagesc(x,y,sst');axis xy;
% % % hold on;
% % % keyboard;
% 
% 
% mx2=model.grid.MX(1)/2;
% x = 1/2:(mx2-1/2);
% x  = x / (2*mx2);
% x = [ -x(end:-1:1) x];
% [x,y]=ndgrid(x,y);
% 
% Lx = model.grid.dX(1) * model.grid.MX(1);
% n_sst_theo = sqrt(pi)*sigma*Lx - 2*pi*sigma^2*Lx/Ly
% % n_sst_theo = sqrt(pi)*sigma*Lx  % without averaging
% n_sst_exp = sum(sst(:).^2)*prod(model.grid.dX)
% % T0=real(ifft2(fft2(sst)));
% % n_T0_exp = sum(T0(:).^2)*prod(model.grid.dX)
% 
% 
% grad_sst = gradient_perso(model.grid, sst);
% grady = grad_sst(:,:,2);
% grady_theo = - y/sigma^2 .* sst;
% % imagesc(abs((grady-grady_theo)./grady_theo)');axis xy;
% % imagesc(grady_theo');axis xy;
% % % [dxsst,dysst]=gradient(sst',model.grid.dX(1),model.grid.dX(2));
% n_grad_sst_theo = sqrt(pi)/(2*sigma)*Lx
% % n_grad_sst_exp = (sum(dxsst(:).^2)+sum(dysst(:).^2)) ...
% %     *prod(model.grid.dX)
% n_grad_sst_exp = sum(grad_sst(:).^2) ...
%     *prod(model.grid.dX)

