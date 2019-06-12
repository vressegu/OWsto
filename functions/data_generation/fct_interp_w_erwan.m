function w=fct_interp_w_erwan(model)
% Filter the SQG temperature field
% 

load(['/Users/vressegu/Documents/manu/matlab_files/erwan1sqg' ...
    num2str(model.grid.MX(1)) 'all.mat'],'w');
w=squeeze(w(:,:,model.id_t,:));
% w=permute(w,[2 1 3]);

dx=model.grid.dX(1);
dy=model.grid.dX(2);
[Mx,My,d]=size(w);

figure;imagesc(w(:,:,1).^2+w(:,:,2).^2);axis xy;

% filter
x=dx*(0:(Mx-1));
y=dy*(0:(My-1));
[x,y]=ndgrid(x,y);
r21=(x-x(1)).^2+(y-y(1)).^2;
r22=(x-(x(end)+dx)).^2+(y-y(1)).^2;
r23=(x-(x(1))).^2+(y-(y(end)+dy)).^2;
r24=(x-(x(end)+dx)).^2+(y-(y(end)+dy)).^2;
% r2=x.^2+y.^2;
% % % sigma_filter= 2^(-3);
sigma_filter= dx*2^3
% sigma_filter= 2*pi*2^(-5)
% % sigma_filter= 2*pi*2^(-4)
dx
filter = exp(-r21/(2*sigma_filter^2)) ... 
        + exp(-r22/(2*sigma_filter^2))... 
        + exp(-r23/(2*sigma_filter^2))... 
        + exp(-r24/(2*sigma_filter^2));
% filter = exp(-r2/(2*sigma_filter^2));
filter=filter/sum(filter(:));
% figure; imagesc((filter));
filter = fft2(filter);

w(:,:,1)=fft2(w(:,:,1)).*filter;
w(:,:,2)=fft2(w(:,:,2)).*filter;
w(:,:,1)=ifft2(w(:,:,1));
w(:,:,2)=ifft2(w(:,:,2));

figure;imagesc(w(:,:,1).^2+w(:,:,2).^2);axis xy;
keyboard;

