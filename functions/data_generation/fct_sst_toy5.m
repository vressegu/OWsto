function sst=fct_sst_toy5(model)
% Create a simple temperature gradient
%

model.odg_b=1;
model.k_inf_on_k1 = 4;
model.slope_b_ini = -5;
sst = init_Spectrum(model);

figure;imagesc(sst');axis xy;axis equal;
% keyboard;

% odg_b = 1e-3;
% 
% [sigma_on_sq_dt,f_sigma,a0_on_dt,spectrum_sigma,spectrum] = fct_random_tracer(model,ft,ft2)
% 
% mx2=model.grid.MX(1)/2;
% my2=model.grid.MX(2)/2;
% sigma = 0.3;
% y = 1/2:(my2-1/2);
% y  = y /my2;
% y = [ -y(end:-1:1) y];
% x = 1/2:(mx2-1/2);
% x  = x /mx2;
% x = [ -x(end:-1:1) x];
% [x,y]=ndgrid(x,y);
% sst = exp( - 1/2 * (x.^2 + y.^2) / sigma^2 );
% 
% sst = odg_b * (2*sst-1);
% % % sst=[ones(1,my2-n) zeros(1,2*n) ones(1,model.grid.MX(2)-my2-n)];
% % % % sst=[ones(1,my2) zeros(1,model.grid.MX(2)-my2)];
% % sst=repmat(sst,[model.grid.MX(1) 1]);
% 
% % x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% % y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% % imagesc(x,y,sst');axis xy;
% % hold on;
% % keyboard;