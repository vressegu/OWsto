function sst=fct_sst_toy(model)
% Create a simple temperature gradient
%

sst=1/model.grid.MX(2)*(0:(model.grid.MX(2)-1));
my2=ceil(model.grid.MX(2)/2);
n=10;
sst=[ones(1,my2-n) zeros(1,2*n) ones(1,model.grid.MX(2)-my2-n)];
% sst=[ones(1,my2) zeros(1,model.grid.MX(2)-my2)];
sst=repmat(sst(end:-1:1),[model.grid.MX(1) 1]);

% x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
% y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
% imagesc(x,y,sst');axis xy;
% hold on;
% keyboard;