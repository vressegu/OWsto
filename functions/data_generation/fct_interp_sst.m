function sst_interp=fct_interp_sst(model,XP)


time_step_choose=5;
load('data/simu2008_30.mat','sst','dx','dy');
% load('data/simu2008_30.mat','sst','dx','dy');
sstref=sst(3:end,:,time_step_choose);
clear sst;
sstref=sstref-mean(mean(sstref));

sstref=sstref';

nsub=2^7;

[Mx,My]=size(sstref);
xref=dx*(-Mx:nsub:2*Mx-1);
yref=dy*(-My:nsub:2*My-1);
% xref=dx*(0:nsub:Mx-1);
% yref=dy*(0:nsub:My-1);
sstref=sstref(1:nsub:end,1:nsub:end);
% sstref=[ sstref sstref sstref ];
sstref=[ sstref; sstref; sstref ];
msstref = repmat(sstref(:,1),[1 My/nsub]);
Msstref = repmat(sstref(:,end),[1 My/nsub]);
sstref=[ msstref sstref Msstref ];

XP=double(XP);
% x=dx*(0:MXI(1)-1);
% y=dy*(0:MXI(1)-1);

sst_interp = interp2(xref,yref,sstref', XP(:,1), XP(:,2),'spline');

% image(10*reshape(sst_interp,[Mx/2 My/2])+20);
% keyboard;