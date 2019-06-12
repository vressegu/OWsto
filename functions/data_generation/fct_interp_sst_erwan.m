function temp=fct_interp_sst_erwan(model)
% Filter the SQG temperature field
% 

load(['/Users/vressegu/Documents/manu/matlab_files/erwan1sqg' ...
    num2str(model.grid.MX(1)) 'all.mat'],'temp');
temp=temp(:,:,model.id_t);
temp=temp-mean(mean(temp));
% temp=temp';

dx=model.grid.dX(1);
dy=model.grid.dX(2);
[Mx,My]=size(temp);

figure;imagesc(temp);axis xy;

% filter
x=dx*(0:(Mx-1));
y=dy*(0:(My-1));
[x,y]=ndgrid(x,y);
r21=(x-x(1)).^2+(y-y(1)).^2;
r22=(x-(x(end)+dx)).^2+(y-y(1)).^2;
r23=(x-(x(1))).^2+(y-(y(end)+dy)).^2;
r24=(x-(x(end)+dx)).^2+(y-(y(end)+dy)).^2;
% r2=x.^2+y.^2;
sigma_filter= 2*pi*2^(-2);
% sigma_filter= 2*pi*2^(-3);
% % sigma_filter= 2*pi*2^(-4);
filter = exp(-r21/(2*sigma_filter^2)) ... 
        + exp(-r22/(2*sigma_filter^2))... 
        + exp(-r23/(2*sigma_filter^2))... 
        + exp(-r24/(2*sigma_filter^2));
% filter = exp(-r2/(2*sigma_filter^2));
filter=filter/sum(filter(:));
% figure; imagesc((filter));
filter = fft2(filter);

temp=fft2(temp);
temp=temp.*filter;
temp=ifft2(temp);

figure;imagesc(temp);axis xy;

