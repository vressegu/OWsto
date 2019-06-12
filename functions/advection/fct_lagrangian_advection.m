function [fft_T_adv,model,X,X0] = fct_lagrangian_advection(model, fft_T0, w)
% function [fft_T_adv,X,X0] = fct_lagrangian_advection(model, fft_T0, w)
% Advection of T with the speed w, with lagrangian transport
%

%% Folder to save plots and files
% model.folder.folder_simu = [ 'images/SQG_MU/' model.type_data ];

% Create the folders
fct_create_folder_plots(model)

% Colormap
load('BuYlRd.mat');
model.folder.colormap = BuYlRd; clear BuYlRd

% Create the folders for filtered field
model_filtered = model;
model_filtered.folder.folder_simu = [ model.folder.folder_simu ...
    '/filtered'];
% model_filtered.folder.folder_simu = [ 'images/SQG_MU/' model.type_data '/filtered'];
fct_create_folder_filtered_plots(model_filtered)

% Create the folders for anisotropic filtered field
% model_filtered_aniso = model;
% model_filtered.folder.folder_simu = [ model.folder.folder_simu ...
%     '/filtered_aniso'];
% % model_filtered_aniso.folder.folder_simu = [ 'images/SQG_MU/' model.type_data '/filtered_aniso'];
% fct_create_folder_filtered_plots(model_filtered_aniso)

N_ech = 1;

%%
%         v_ortho = w(:,:,[2 1]);
%         v_ortho(:,:,1) = - v_ortho(:,:,1);
%         v_ortho=v_ortho(3:end-2,3:end-2,:);
% %         v_ortho(2:end-1,2:end-1,:)=v_ortho(2:end-1,2:end-1,:);
%
%         vort = 1/(2*model.grid.dX(1)) ...
%             *(w(3:end,2:end-1,2)-w(2:end-1,1:end-2,2)) ...
%                 - 1/(2*model.grid.dX(2)) ...
%             *(w(2:end-1,3:end,1)-w(2:end-1,1:end-2,1));% Mx-2 My-2
%         nabla_vort(:,:,1) = 1/(2*model.grid.dX(1)) ...
%             *(vort(3:end,2:end-1)-vort(2:end-1,1:end-2)); % Mx-4 My-4 2
%         nabla_vort(:,:,2) = 1/(2*model.grid.dX(2)) ...
%             *(vort(2:end-1,3:end)-vort(2:end-1,1:end-2));% Mx-4 My-4 2
%         sq_a = sum(v_ortho.*nabla_vort,3);
%
%         model.l_pixel_nugget=10e3;
%         Ld=exp(-model.kriging.param(2));
%         l_cri = model.l_pixel_nugget;
% %         l_cri = sqrt(prod(model.gris.dX));
%         time_adv = 2*pi* sqrt(sum(v_ortho.^2,3))/l_cri ./ abs(sq_a);
%         %         time_adv = time_adv/3600/24;
%         time_adv = min(time_adv(:))
%
%         model.advection.advection_duration=time_adv;

model.advection.advection_duration=fct_time_adv(model,w);
% % model.advection.advection_duration=2*model.advection.advection_duration;
% model.advection.advection_duration=model.advection.advection_duration*10;

%         if t*dt >time_adv
%             keyboard;
%         end
% %         figure(1)
% %         imagesc(x,y,time_adv(:,1:My)');
% %         keyboard;

%%

My = size(fft_T0,2);
if model.mirror
    My=My/2;
end

dt=model.advection.dt_adv;
N_t = ceil(model.advection.advection_duration/dt);
fprintf(['Time of advection : ' num2str(N_t*dt/3600/24) ' days \n']);

T0=real(ifft2(fft_T0));

% on_tau_local = fct_on_tau(model,w_eul,w_Lag);

% T0 = zeros(size(T0));
% T0(:,1:ceil(size(T0,2)/3))=1;


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
% X0ref{1}=x;
% X0ref{2}=y;
[x,y]=ndgrid(x,y);
X(:,1)=x(:);
X(:,2)=y(:);
X_forward = X;
% Xref=X;

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

%% Stretching time scale
[on_tau2_global_folding,on_tau_local_folding,rate_switch] =...
    fct_on_tau_2_folding(model,w);
[on_tau2_global_stretching, on_tau_local_stretching] = ...
     fct_on_tau_2_stretching(model,w);
on_tau2_global = fct_plot_tau(model,...
    on_tau_local_folding,on_tau_local_stretching,rate_switch);
% keyboard;

%% Lagrangian advection

% Xref=X;
w_ini = w;
w = fct_mirror_on_border(w,nbp);

grad_T0 = gradient_perso(model.grid, T0);
norm_T0 = sum(T0(:).^2).* prod(model.grid.dX);
norm_grad_T0 = sum(grad_T0(:).^2).*prod(model.grid.dX);
model.spectrum_theo0.coef1 = 2 * sqrt( (2*pi)^3 * norm_T0^3 / norm_grad_T0 );
model.spectrum_theo0.coef2 =  norm_T0 / norm_grad_T0 ;

m_H_T0 = bsxfun( @times, grad_T0 , ...
    permute( grad_T0, [ 1 2 4 3] ) );
m_H_T0 = squeeze(sum(sum(m_H_T0,2),1)) .* prod(model.grid.dX);

model.spectrum_theo_aniso0.coef1 = 2*pi * norm_T0^2 / sqrt(det(m_H_T0));
model.spectrum_theo_aniso0.coef2 =  norm_T0 * (m_H_T0'*m_H_T0)\m_H_T0' ;

T0=[T0(end-nbp+1:end,:);T0;T0(1:nbp,:)];
T0=[T0(:,end-nbp+1:end) T0 T0(:,1:nbp)];
% T0=[repmat(T0(:,1),[1 nbp]) T0 repmat(T0(:,end),[1 nbp])];

if strcmp(model.type_data,'erwan')
    true_temp=fct_real_sst_erwan(model);
end

tt_last=-inf;
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
if model.mirror
    y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
else
    y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
end

%% Choice of the variance tensor a
% Variance tensor
model.sigma.a0 = 2 * model.physical_constant.f0 / model.sigma.k_c^2;
% Diffusion coefficient
model.advection.coef_diff = 1/2 * model.sigma.a0;

%% Fourier transform of the kernel \tilde sigma
if ~isinf(model.sigma.k_c)
    % Fourier transform of the kernel \tilde sigma up to a multiplicative
    % constant
    [sigma_on_sq_dt, ~, missed_var_small_scale_spectrum ] ...
        = fct_sigma_spectrum(model,fft2(w_ini));
    clear w_ini;
    % sigma_on_sq_dt will be used to simulate sigma d B_t
    % missed_var_small_scale_spectrum will be used to set the mulitplicative
    % constant
    
    % Muliplicative constant of the kernel \tilde sigma
    model.sigma.a0_on_dt = model.sigma.a0 / dt;
    sigma_on_sq_dt = sqrt(2*model.sigma.a0_on_dt/missed_var_small_scale_spectrum) ...
        * sigma_on_sq_dt;
    % the factor d=2 is for the dimension d of the space R^d
    % the variance of sigma_dBt/dt is tr(a)/dt = 2 a0 /dt
    clear missed_var_small_scale_spectrum
else
    sigma_on_sq_dt = 0;
end



% fct_plot(model,fft_T0,'0');

dot_red = bsxfun(@plus, [mean(x) mean(y)] , ...
    [ 2.40 0 ; 2.60 0 ; 2.80 0]*1e4 );
%     [ 2.20 0 ; 2.50 0 ; 2.80 0]*1e4 );
%     [ 2.19 0 ; 2.47 0 ; 2.72 0]*1e4 );
dot_red_savex=dot_red(:,1);
dot_red_savey=dot_red(:,2);

v_norm_T = norm_T0;
v_norm_grad_T = norm_grad_T0;
v_norm_grad_T_estim = norm_grad_T0;
v_alpha2_m = 0;
v_alpha2_m_estim = 0;

vect_t_plot = 0;
dt = model.advection.dt_adv;

w0 = w;
w_ini_transp = permute( w_ini , [ 1 2 4 3]);
cov_st = squeeze( sum(sum( bsxfun(@times, w_ini , w_ini_transp ) , 2),1) );
cov_scalaire_st = cov_st(1,1) + cov_st(2,2);

for t=1:N_t
    
    if isinf(model.sigma.k_c) % Deterministic case
        sigma_dBt_dt = 0;
    else % Stochastic case
        % Fourier transform of white noise
        dBt_C_on_sq_dt = fft2( randn( [ model.grid.MX 1 N_ech]));
        % Multiplication by the Fourier transform of the kernel \tilde \sigma
        fft_sigma_dBt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
        clear dBt_C_on_sq_dt
        % Homogeneous velocity field
        sigma_dBt_dt = real(ifft2(fft_sigma_dBt));
        clear fft_sigma_dBt
        sigma_dBt_dt = fct_mirror_on_border(sigma_dBt_dt,nbp);
        w = w + sigma_dBt_dt;
    end
    
    % Backward advection (to get T(x,t) on a regular grid)
    [X,dX] = RK4_advection_lagrangienne(model,X, -w, X0);
    
%     Forward advection (to get T(x_0,t) on a regular grid)
    [X_forward,dX_forward] = ...
        RK4_advection_lagrangienne(model,X_forward, +w, X0);
    
    dot_red = RK4_advection_lagrangienne(model,dot_red, w, X0);
    clear w_adv
    w_adv(:,:,1) = interp2(X0{1},X0{2},w0(:,:,1)',X(:,1),X(:,2));
    w_adv(:,:,2) = interp2(X0{1},X0{2},w0(:,:,2)',X(:,1),X(:,2));
    w_adv = reshape(w_adv,[MX 2]);
    w_adv = permute( w_adv , [1 2 4 3]);
    
%     clear w_adv_forward
%     w_adv_forward(:,:,1) = interp2(X_forward(:,1),X_forward(:,2),...
%         w0(:,:,1)',X0{1},X0{2});
%     w_adv_forward(:,:,2) = interp2(X_forward(:,1),X_forward(:,2),...
%         w0(:,:,2)',X0{1},X0{2});
%     w_adv_forward = reshape(w_adv_forward,[MX 2]);
%     w_adv_forward = permute( w_adv_forward , [1 2 4 3]);
    
    cov_st_ = squeeze( sum(sum( bsxfun(@times, w_ini , w_adv ), 2),1) );
    cov_st_ = cov_st_ / prod(model.grid.dX .* model.grid.MX);
%     clear w_adv
    cov_st(:,:,t+1) = cov_st_ ;
    cov_scalaire_st_ = cov_st_(1,1) + cov_st_(2,2);
    cov_scalaire_st(t+1) = cov_scalaire_st_;
    
    tt = floor(t *dt/ (3600*24)); % Number of days
    if tt > tt_last
        tt_last = tt;
%     if (t - t_last_plot) >= 100
%         %     if (t - t_last_plot) >= 10
        t_last_plot=t;
        day = num2str(floor(t*dt/(3600*24)));
        vect_t_plot = [vect_t_plot t];
        
        %% Plot of the covariance of the lagrangian velocity
%         fct_plot_cov_v(model,cov_scalaire_st,t,dt);
        
        
        %         Xplot=reshape(X,[MX 2]);
        %         Xplot=Xplot(1:10:end,1:20:end,:);
        %         Xplot=reshape(Xplot,[size(Xplot,1)*size(Xplot,2) 2]);
        %         Xplotref=reshape(Xref,[MX 2]);
        %         Xplotref=Xplotref(1:10:end,1:20:end,:);
        %         Xplotref=reshape(Xplotref,[size(Xplotref,1)*size(Xplotref,2) 2]);
%         cov_st = 
        
        %% Tracer 
        
        % Tracer interpolation
        T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
        T_adv=reshape(T_adv,MX);
        fft_T_adv=fft2(T_adv);
        
        % Plots of the tracer
        fct_plot(model,fft_T_adv,day);
        
        % Tracer interpolation forward
%         T_adv_forward = griddata(X_forward(:,1),X_forward(:,2),T0',...
%             X0{1},X0{2});
% %         T_adv_forward = TriScatteredInterp(X_forward(:,1),X_forward(:,2),T0',...
% %             X0{1},X0{2});
% %         T_adv_forward = interp2(X_forward(:,1),X_forward(:,2),T0',...
% %             X0{1},X0{2});
%         T_adv_forward=reshape(T_adv_forward,MX);
%         fft_T_adv_forward=fft2(T_adv_forward);
        
%         % Plots of the tracer forward
%         fct_plot(model,fft_T_adv_forward,day);
        
        % Statistic description of the tracer
        grad_T_adv = gradient_perso(model.grid, T_adv);
        %         warning('centering to compute coefficients');
        %         T_adv_centered = T_adv - mean(T_adv(:));
        T_adv_centered = T_adv;
        
        norm_T = sum(T_adv_centered(:).^2) .* prod(model.grid.dX);
        norm_grad_T = sum(grad_T_adv(:).^2) .* prod(model.grid.dX);
        m_H_T = bsxfun( @times, grad_T_adv , ...
            permute( grad_T_adv, [ 1 2 4 3] ) );
        m_H_T = squeeze(sum(sum(m_H_T,2),1)) .* prod(model.grid.dX);
        %         trace = m_H_T(1,1)+m_H_T(2,2)
        %         norm_grad_T
        %         keyboard;
        
        norm_grad_T_estim = norm_grad_T0 * ( 1+ (t *dt)^2 * on_tau2_global );

        model.spectrum_theo.coef1 = ...
            2 * sqrt( (2*pi)^3 * norm_T0^3 / norm_grad_T_estim );
        model.spectrum_theo.coef2 =  norm_T0 / norm_grad_T_estim ;
%         model.spectrum_theo.coef1 = ...
%             2 * sqrt( (2*pi)^3 * norm_T^3 / norm_grad_T_estim );
%         model.spectrum_theo.coef2 =  norm_T / norm_grad_T_estim ;
        
%         model.spectrum_theo.coef1 = ...
%             2 * sqrt( (2*pi)^3 * norm_T^3 / norm_grad_T );
%         model.spectrum_theo.coef2 =  norm_T / norm_grad_T ;
        
        model.spectrum_theo_aniso.coef1 = 2*pi * norm_T^2 / sqrt(det(m_H_T));
        model.spectrum_theo_aniso.coef2 =  norm_T * (m_H_T'*m_H_T)\m_H_T' ;
        
        v_norm_T = [v_norm_T norm_T];
        v_norm_grad_T = [v_norm_grad_T norm_grad_T];
        v_norm_grad_T_estim = [v_norm_grad_T_estim norm_grad_T_estim];
        
%         %         figure(20);
%         %         plot(vect_t_plot*dt,v_norm_T);
%         % %         hold on;
%         % %         plot(vect_t_plot*dt,norm_T0*ones(length(vect_t_plot)));
%         % %         hold off;
%         figure(21);
%         loglog(vect_t_plot*dt,v_norm_grad_T);
%         %         plot(vect_t_plot*dt,v_norm_grad_T);
%         % %         hold on;
%         % %         plot(vect_t_plot*dt,norm_grad_T0*ones(length(vect_t_plot)));
%         % %         hold off;
        
        %         keyboard;
        
        %% Stretching and Mixing dignostic
        
        % Gradient of the inverse flow
        nabla_phi = fct_nabla_phi(model,X);
        
        % Evolution matrix interpolation
        nabla_phi_mem = nabla_phi;

        if model.plot.alpha_adv_back
            s=size(nabla_phi(:,:,1,1));
            for i=1:2
                for j=1:2
                    nabla_phi_adv_back=nabla_phi(:,:,i,j);
                    
                    nabla_phi_adv_back = [ ...
                        nabla_phi_adv_back(end-nbp+1:end,:); ...
                        nabla_phi_adv_back ; ...
                        nabla_phi_adv_back(1:nbp,:) ];
                    nabla_phi_adv_back = [ ...
                        nabla_phi_adv_back(:,end-nbp+1:end) ...
                        nabla_phi_adv_back ...
                        nabla_phi_adv_back(:,1:nbp) ];
                    
%         T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
%         T_adv=reshape(T_adv,MX);
        
%                     nabla_phi_adv_back=nabla_phi_adv_back(:);
                    nabla_phi_adv_back = interp2(X0{1},X0{2},...
                        nabla_phi_adv_back',X_forward(:,1),X_forward(:,2));
                    nabla_phi_adv_back=reshape(nabla_phi_adv_back,s);
                    nabla_phi(:,:,i,j)=nabla_phi_adv_back;
                end
            end
        end
        
        % Compute alpha, beta, etc
        [alpha2_m, cax_alpha] = fct_mezic2(model,nabla_phi,t);
        alpha2_m_estim = (t *dt)^2 * on_tau2_global;
%         alpha2_m = fct_mezic2(model,X,t);
% %         fct_mezic3(model,X,t,dot_red);
        v_alpha2_m = [v_alpha2_m alpha2_m];
        v_alpha2_m_estim = [v_alpha2_m_estim alpha2_m_estim];
        %         figure(22);
        %         loglog(vect_t_plot*dt,v_alpha2_m);
        %         %         plot(vect_t_plot*dt,v_alpha2_m./(vect_t_plot.^2));
        %         %         plot(vect_t_plot*dt,v_alpha2_m);
        
        if eval(day) < 3
            on_tau2_global = fct_plot_tau(model,...
                on_tau_local_folding,on_tau_local_stretching,...
                rate_switch,cax_alpha);
        end
        
        %% Stretching time scale
%         on_tau2_global = fct_on_tau_2(model,w_ini);
%         on_tau_local = fct_on_tau(model,nabla_phi,w_ini,w_adv);
        
%         on_tau_local = fct_on_tau_forward(model,nabla_phi,w_ini,w_adv);
        
        %% Plot of time evolution of tracer and alpha
        fct_plot_tracer_alpha_comp(model,v_alpha2_m,v_alpha2_m_estim,...
            v_norm_grad_T/norm_grad_T0,...
            v_norm_grad_T_estim/norm_grad_T0,dt,vect_t_plot);   
%         fct_plot_tracer_alpha(model,v_alpha2_m,v_norm_grad_T,...
%                                             norm_grad_T0,dt,vect_t_plot);        
        
        %% Anisotropic filtering of the tracer
        %         % Filter anisotropic
        %         filter_aniso = fct_design_filter_aniso(model);
        %         fft_T_filtered_aniso = fft_T_adv .* filter_aniso;
        % %         figure;imagesc(filter)
        % %         figure;imagesc(real(ifft2(fft_T_filtered)))
        %         fct_plot(model_filtered_aniso,fft_T_filtered_aniso,day);
        
        
        %% Filtered tracer
        if model.plot.filtered_field
            % Filter isotropic
            filter = fct_design_filter(model);
            fft_T_filtered = fft_T_adv .* filter;
            %         figure;imagesc(filter)
            %         figure;imagesc(real(ifft2(fft_T_filtered)))
            fct_plot(model_filtered,fft_T_filtered,day);
        end
        
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
        
        
        %         figure;
        %         plot3(dot_red_savex(1,:),dot_red_savey(1,:),(1:(size(dot_red_savey,2))));
        % %         keyboard;
        
        
    end
end

%%
% % save('dot_red_save_draft.mat','dot_red_savex','dot_red_savey');
% load('dot_red_save_draft.mat');
%
% figure;
% %         dot_red_savex(1,:)=dot_red_savex(1,:)-dot_red_savex(2,:);
% %         dot_red_savex(3,:)=dot_red_savex(3,:)-dot_red_savex(2,:);
% %         dot_red_savey(1,:)=dot_red_savey(1,:)-dot_red_savey(2,:);
% %         dot_red_savey(3,:)=dot_red_savey(3,:)-dot_red_savey(2,:);
% %         dot_red_savex(2,:)=zeros(size(dot_red_savex(2,:)));
% %         dot_red_savey(2,:)=zeros(size(dot_red_savey(2,:)));
% nn=1e3;
% %         plot3(dot_red_savex',dot_red_savey',repmat((1:(size(dot_red_savey,2))),[3 1])');
% plot3(dot_red_savex(:,1:nn)',dot_red_savey(:,1:nn)',repmat((1:nn),[3 1])');
% %         plot3(dot_red_savex',dot_red_savey',repmat((1:(size(dot_red_savey,2))),[3 1])');
% % %         hold on;
% % %         plot3(dot_red_savex(1,:),dot_red_savey(1,:),(1:(size(dot_red_savey,2))));
% keyboard;
%%


% T_adv = interp2(X(:,1),X(:,2),T0',X0{1},X0{2});
T_adv = interp2(X0{1},X0{2},T0',X(:,1),X(:,2));
%         T_adv = interp2(X0(:,1),X0(:,2),T0,X(:,1),X(:,2));
T_adv=reshape(T_adv,MX);

fft_adv=fft2(T_adv);