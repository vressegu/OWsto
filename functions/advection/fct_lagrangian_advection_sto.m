function [fft_T_adv] = fct_lagrangian_advection_sto(model, fft_T0, w)
% Advection of T with the speed w, with lagrangian transport
%

My = size(fft_T0,2);
if model.mirror
    My=My/2;
end

dt=model.advection.dt_adv;
N_t = ceil(model.advection.advection_duration/dt);

T0=real(ifft2(fft_T0));

% %% Anti aliasing filter
%
% fft_T0= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',fft_T0);
% fft_T0= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,fft_T0);
% w=fft2(w);
% w(:,:,1)= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',w(:,:,1));
% w(:,:,1)= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,w(:,:,1));
% w(:,:,2)= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',w(:,:,2));
% w(:,:,2)= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,w(:,:,2));
% w=real(ifft2(w));

%% Grid
% Moving grid
MX=model.grid.MX;
x=model.grid.dX(1)* (0:(MX(1)-1)) ;
y=model.grid.dX(2)* (0:(MX(2)-1)) ;
X0ref{1}=x;
X0ref{2}=y;
[x,y]=ndgrid(x,y);
X(:,1)=x(:);
X(:,2)=y(:);
Xref=X;
% Reference grid
nbp=3;
X0{1}=model.grid.dX(1)* ((-nbp):(MX(1)+nbp-1)) ;
X0{2}=model.grid.dX(2)* ((-nbp):(MX(2)+nbp-1)) ;
% x=model.grid.MX(1)* ((-nbp):(MX(1)+nbp-1)) ;
% y=model.grid.MX(2)* ((-nbp):(MX(2)+nbp-1)) ;
% % x=model.grid.MX(1)* ((-MX(1)):(2*MX(1)-1)) ;
% % Py=floor(MX(2)/2);
% % y=model.grid.MX(2)* ((-Py):((MX(2)+Py)-1)) ;
% [x,y]=ndgrid(x,y);
% X0(:,1)=x(:);
% X0(:,2)=y(:);

%% Lagrangian advection

% Xref=X;
w=w(:,:,1)+1i*w(:,:,2);
w=[w(end-nbp+1:end,:);w;w(1:nbp,:)];
w=[repmat(w(:,1),[1 nbp]) w repmat(w(:,end),[1 nbp])];
% w=reshape(w,[(MX(2)+2*nbp)*(MX(1)+2*nbp)]);
% w=[w;w;w];
% w=[w(:,end-Py+1:end) w w(:,1:Py)];
% % w=reshape(w,[(MX(2)+2*Py)*3*MX(1) ]);
% w=w(:);
w(:,:,2)=imag(w);
w(:,:,1)=real(w(:,:,1));

T0_ini=T0;
T0=[T0(end-nbp+1:end,:);T0;T0(1:nbp,:)];
T0=[repmat(T0(:,1),[1 nbp]) T0 repmat(T0(:,end),[1 nbp])];
% T0=reshape(T0,[(MX(2)+2*nbp)*(MX(1)+2*nbp)]);
% T0=[T0;T0;T0];
% T0=[T0(:,end-Py+1:end) T0 T0(:,1:Py)];
% % T0=reshape(T0,[(MX(2)+2*Py)*3*MX(1) ]);
% T0=T0(:);

T_adv = T0_ini;

%%
% warning('wfs')
% N_t=0;
%%

for t=1:N_t
    [X,dX] = RK4_advection_lagrangienne(model,X, -w, X0);
end

phi_X0 = X - Xref;
delta_t = N_t*model.advection.dt_adv;
w_moy = 1/delta_t * phi_X0;
a_T=zeros([size(X) 2]);
X=Xref;
pause;
for t=1:N_t
    [X,dX] = RK4_advection_lagrangienne(model,X, -w, X0);
    dX_moy = -w_moy * model.advection.dt_adv;
    dX_sigma = dX - dX_moy;
    a_T = a_T + bsxfun(@times,dX,permute(dX,[1 3 2]));
    delta_t = t*model.advection.dt_adv;
    a=1/delta_t*a_T;
    tr_a = a(:,1,1)+a(:,2,2);
    det_a = a(:,1,1)./a(:,2,2);
    %     det_a = a(:,1,1).*a(:,2,2);
    
    if mod(t,100)==0
        %     if mod(t*dt,3600*24*0.1)==0
        % Plot
        tr_a=reshape(tr_a,[MX ]);
        det_a=reshape(det_a,[MX ]);
        
        tr_a2 = interp2(X0ref{1},X0ref{2},tr_a',X(:,1),X(:,2));
        tr_a2=reshape(tr_a2,MX);
        
        
        figure(3)
        subplot(3,1,1)
        image(5e-3*tr_a(:,1:My)')
        subplot(3,1,2)
        image(5e-3*tr_a2(:,1:My)')
        subplot(3,1,3)
        image(30*det_a(:,1:My)')
        %         image(5e-6*det_a(:,1:My)')
        
        
        Xplot=reshape(X,[MX 2]);
        Xplot=Xplot(1:5:end,1:20:end,:);
        Xplot=reshape(Xplot,[size(Xplot,1)*size(Xplot,2) 2]);
        Xplotref=reshape(Xref,[MX 2]);
        Xplotref=Xplotref(1:5:end,1:20:end,:);
        Xplotref=reshape(Xplotref,[size(Xplotref,1)*size(Xplotref,2) 2]);
        
        T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
        %         tr_a_2 = interp2(X0ref{1},X0ref{2},tr_a',X(:,1),X(:,2));
        %T_adv = interp2(X0(:,1),X0(:,2),T0,X(:,1),X(:,2));
        T_adv=reshape(T_adv,MX);
        %         tr_a2=reshape(tra_2,MX);
        fft_T_adv=fft2(T_adv);
        
        figure(4)
        subplot(2,1,1)
        fct_spectrum(model, fft_T0);
        hold on;
        fct_spectrum(model,fft_T_adv);
        hold off
        ax = axis;
        ax(3)=0;
        axis(ax);
        
        subplot(2,1,2)
        image(10*T_adv(:,1:My)'+20)
        
        figure(5)
%         subplot(3,1,3)
        plot(Xplotref(:,1),Xplotref(:,2),'.')
        hold on;
        plot(Xplot(:,1),Xplot(:,2),'.r')
        hold off
        ax=axis;
        ax(4)=ax(4)/2;
        axis(ax);
        drawnow
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        
    end
end
%%
% T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
% %         T_adv = interp2(X0(:,1),X0(:,2),T0,X(:,1),X(:,2));
% T_adv=reshape(T_adv,MX);

%%

PX=model.grid.MX/2;
kx=1/(model.grid.MX(1)*model.grid.dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
ky=1/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
% kx=1/(model.grid.MX(1)*model.grid.dX(1))*[ 0: (MX(1)-1)] ;
% ky=1/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(MX(2)-1) ];
a_T=reshape(a_T,[MX 2 2]);
if model.advection.a_average
    a_T = 1/prod(MX) * squeeze(sum(sum(a_T,2),1))
    T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
    %         T_adv = interp2(X0(:,1),X0(:,2),T0,X(:,1),X(:,2));
    T_adv=reshape(T_adv,MX);
    fft_T_adv=fft2(T_adv);
    [kx,ky]=ndgrid(kx,ky);
    fft_T_adv = fft_T_adv .* exp( -1/2 * (2*pi)^2 * ( ...
                2*kx.*ky*a_T(1,2) ...
                + kx.^2*a_T(1,1) + ky.^2*a_T(2,2) ));
            
else
%     PX=model.grid.MX/2;
%     kx=1/(model.grid.MX(1)*model.grid.dX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
%     ky=1/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
%     % kx=1/(model.grid.MX(1)*model.grid.dX(1))*[ 0: (MX(1)-1)] ;
%     % ky=1/(model.grid.MX(2)*model.grid.dX(2))*[ 0:(MX(2)-1) ];
    X=reshape(X,[MX 2]);
%     a_T=reshape(a_T,[MX 2 2]);
    fft_T_adv=zeros(size(T0_ini));
    warning('the jacobian should be taken into account');
    for id_kx=1:MX(1)
        for id_ky=1:MX(2)
            %         fft_temp = T0_ini .* exp(-1i*2*pi* ...
            %             (kx(id_kx)*X(:,:,1)+ky(id_ky)*X(:,:,2)));
            fft_temp = T0_ini .* exp( -1i*2*pi* ...
                (kx(id_kx)*X(:,:,1)+ky(id_ky)*X(:,:,2)) ...
                -1/2 * (2*pi)^2 * ( ...
                2*kx(id_kx)*ky(id_ky)*a_T(:,:,1,2) ...
                + kx(id_kx)^2*a_T(:,:,1,1) + ky(id_ky)^2*a_T(:,:,2,2) ));
            fft_T_adv(id_kx,id_ky) = sum(sum(fft_temp));
        end
    end
end
fft_T_adv(PX(1)+1,:)=0;
fft_T_adv(:,PX(2)+1)=0;
%                   n=1
%%

% fft_adv=fft2(T_adv);

% % fft_adv=fft_adv(1:5,1:5)
% % fft_T=fft_T(1:5,1:5)
%
% % norm_table((fft_adv-fft_T)./fft_adv)
% % norm_table((fft_adv-fft_T)./fft_adv)/numel(fft_adv)
%
% err=abs((fft_adv-fft_T)./fft_adv);
% idx = abs(fft_adv(:));
% idx = (idx<max(idx)*1e-3);
% err(idx)=0;
%
% norm_table(err)
% norm_table(err)/numel(err)
%
% keyboard;