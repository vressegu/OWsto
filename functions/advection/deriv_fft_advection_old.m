function d_fft_T_adv = deriv_fft_advection_old(model, fft_T, w)
% function d_fft_T_adv = deriv_fft_advection(model, fft_T, fft_w)
% Advection of T with the speed w, in Fourrier space
%

meth_anti_alias=model.advection.meth_anti_alias;

%% Grid
M=prod(model.grid.MX);
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=model.grid.MX/2;
kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
    .* fct_unity_approx5(model.grid.MX(1));
ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
    .* fct_unity_approx5(model.grid.MX(2));
% if strcmp(meth_anti_alias,'deriv_LowPass')
%     kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
%         .* fct_unity_approx5(model.grid.MX(1));
%     ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
%         .* fct_unity_approx5(model.grid.MX(2));
% else
%     kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
%     ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% end
% % % kx=1/(model.grid.MX(1))*[ 0:(PX(1)/2) zeros(1,PX(1)-1) (-PX(1))/2:-1];
% % % ky=1/(model.grid.MX(2))*[ 0:(PX(2)/2) zeros(1,PX(2)-1) (-PX(2))/2:-1];
% % kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
% % ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% % % % kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) -PX(1):-1];
% % % % ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) -PX(2):-1];
% % % % % kx=1/(model.grid.MX(1))*(0:(model.grid.MX(1)-1));
% % % % % ky=1/(model.grid.MX(2))*(0:(model.grid.MX(2)-1));
% % % % % kx=1/(model.grid.dX(1)*model.grid.MX(1))*(0:(model.grid.MX(1)-1));
% % % % % ky=1/(model.grid.dX(2)*model.grid.MX(2))*(0:(model.grid.MX(2)-1));
[kx,ky]=ndgrid(kx,ky);

%%
% w=ifft2(fft_w);
T=ifft2(fft_T);
T=real(T);

wT=bsxfun(@times,w,T);

fft_wT=fft2(wT);clear wT
adv1x=fft_wT(:,:,1);
adv1y=fft_wT(:,:,2);clear fft_wT

% figure(1)
% contour(abs(adv1x))
% figure(2)
% contour(abs(adv1y))
% keyboard;

%% Advection term
% close all

%     % Plot
%     figure(1)
%     pcolor(abs( fft_T ));
%     T_adv = ifft2(fft_T );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
%     figure(2)
%     pcolor(real(T_adv))
% 
%     % Plot
%     figure(4)
%     pcolor(abs( fft_w(:,:,1) ));
%     T_adv = ifft2(fft_w(:,:,1) );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
%     figure(5)
%     pcolor(real(T_adv))
%     
% keyboard;

% adv1x = 1/M * conv2periodic(fft_T,fft_w(:,:,1));

%     % Plot
%     figure(8)
%     pcolor(abs( adv1x ));
% %     T_adv = ifft2(adv1x  );
% %     norm_table(imag( T_adv ))/norm_table(real(T_adv))
% %     figure(9)
% %     pcolor(real(T_adv))
% %     figure(10)
% %     pcolor(imag(T_adv))
% %     figure(11)
% %     pcolor(angle(T_adv))
% % %     figure(12)
% % %     pcolor(angle(adv1x))
% %     norm_table(T_adv - ifft2(fft_T).*ifft2(fft_w(:,:,1)))/norm_table(T_adv)
%     
% keyboard;
    
adv1x = 2*1i*pi/model.grid.dX(1)* kx .* adv1x;
% adv1x = 2*1i*pi* diag(kx(:)) * adv1x(:);

%     % Plot
%     figure(6)
%     contour(abs( adv1x ));
%     T_adv = ifft2(adv1x );
% %     pcolor(abs( reshape(adv1x,model.grid.MX) ));
% %     T_adv = ifft2(reshape(adv1x,model.grid.MX) );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
% %     figure(11)
% %     pcolor(mod(angle(T_adv),pi))
% %     figure(7)
% %     pcolor(real(T_adv))
    
% keyboard;
%     
%     
%     % Plot
%     figure(4)
%     pcolor(abs( fft_w(:,:,2) ));
%     T_adv = ifft2(fft_w(:,:,2) );
%     figure(5)
%     pcolor(real(T_adv))
%     
% keyboard;
    
% adv1y = 1/M * conv2periodic(fft_T,fft_w(:,:,2));

%     % Plot
%     figure(8)
%     pcolor(abs( adv1y ));
%     T_adv = ifft2(adv1y  );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
%     figure(9)
%     pcolor(real(T_adv))
%     
% keyboard;
    
adv1y = 2*1i*pi/model.grid.dX(2)* ky .* adv1y;
% adv1y = 2*1i*pi* diag(ky(:)) * adv1y(:);

%     % Plot
%     figure(4)
%     pcolor(abs( adv1y ));
%     T_adv = ifft2(adv1y );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
%     figure(5)
%     pcolor(real(T_adv))
%     


%% Diffusion
if strcmp(meth_anti_alias,'Lax_Wendroff') || ...
        strcmp(model.advection.type_adv,'geostrophic_adv_filtered') || ...
        strcmp(model.advection.type_adv,'geostrophic_sto_homogene')
    adv2 = model.advection.coef_diff * ...
        - (2*pi)^2 * ( (1/model.grid.dX(1)* kx).^2 + (1/model.grid.dX(2)* ky).^2 ) ...
        .* fft_T;
else
    adv2=0;
end

%%
    
d_fft_T_adv=adv1x+adv1y+adv2; clear adv1x adv1y adv2

%%
if strcmp(meth_anti_alias,'fct_LowPass')
    d_fft_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',d_fft_T_adv);
    d_fft_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,d_fft_T_adv);
end

%%
% % keyboard;
% d_fft_T_adv(PX(1)+1,:)=0;
% d_fft_T_adv(:,PX(2)+1)=0;
% % d_fft_T_adv(PX(1)/2+1:3*PX(1)/2,:)=0;
% % d_fft_T_adv(:,PX(2)/2+1:3*PX(2)/2)=0;
%%

%     % Plot
%     figure(4)
%     pcolor(abs( d_fft_T_adv));
%     T_adv = ifft2( d_fft_T_adv );
%     norm_table(imag( T_adv ))/norm_table(real(T_adv))
%     figure(5)
%     pcolor(real(T_adv))
%     
% keyboard;

% d_fft_T_adv=reshape(d_fft_T_adv,model.grid.MX);

% %% Diffusion term
% if nargin > 3
% end

