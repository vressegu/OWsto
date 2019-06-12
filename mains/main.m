%% Test of fct_a

init;

% choice_a = 'spectrum_geo';
% choice_a = 'averaged_lyapunov';
% choice_a = 'kriging';
choice_a = 'var_w';
% type_data = 'klein'
% type_w = 'klein'
type_data = 'toy2'; 
type_w = 'toy8';
k_c = inf;
% k_c = 1/3;
alpha_adv_back=true;

centered = false;

% Spectrum slope of sigma dBt
slope_sigma = - 5/3; 

% Rate between the smallest and the largest wave number of sigma dBt
kappamin_on_kappamax = 1/10;

plot_filtered_field = true;

param_SQG.bool=false;

type_adv = 'geostrophic_adv_lagrangian_alone';
% type_adv = 'geostrophic_adv_lagrangian_alone' ; % 'geostrophic_adv_filtered';
% 'UQ' 'UQ2'

% 'geostrophic_adv_alone' 'geostrophic_adv_lagrangian_alone'
% 'geostrophic_adv_lagrangian_sto' 'geostrophic_adv_filtered'
% 'geostrophic_sto_homogene' 'geostrophic_sto' (stochastic version of
% geostrophic equilibrium)
% subsamp_x=1;
% subsamp_y=1;
subsamp_x=4;
subsamp_y=4;
% % subsamp_x=2;
% % subsamp_y=2;
% % % subsamp_x=2;
% % % subsamp_y=4;
% % % subsamp_x=1;
% % % subsamp_y=8;

% diff_in_forward=true;
diff_in_forward=false;
recompute_a_w=false;
recompute_kriging_ssh=false;
recompute_kriging_sst=true;
mirror=strcmp(type_data,'klein');

% N_ech=3;
N_ech=200;
bool_HV=true;
% HV_order=8;

% dt_adv=15*60;% ok
% dt_adv=1*60;%aliasing
% dt_adv=100*60;%aliasing

% dt_adv=3600;
dt_adv=5*3600;
% dt_adv=15*60;
% dt_adv=10; % 2 min
% advection_duration=3*24*3600; % 3 jours
% % advection_duration=4*24*3600;
% % advection_duration=16*24*3600; % 16 jours comme guillaume

%% Height kriging

model = fct_w_a_geostrophic_stochastic_kriging(recompute_kriging_ssh);

% load('model_globalp2dx');
% model=model_save;

model.plot.alpha_adv_back = alpha_adv_back;
model.mirror=mirror;clear mirror
model.plot.filtered_field=plot_filtered_field;clear plot_filtered_field
model.advection.type_adv=type_adv;
model.centered= centered;
model.sigma.k_c = k_c;
model.physical_constant.f0 = 1e-4;
model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
model.sigma.slope_sigma = slope_sigma;
% model.folder.folder_simu = [ 'images/SQG_MU/' type_data ];
model.folder.folder_simu = [ 'images/SQG_MU/T_' type_data '_w_' type_w ];

%% Generating LR SST
model.plot_sst_kriging=false;
if strcmp(type_data,'klein')
    if recompute_kriging_sst
        % Grid
        load('data/simu2008_30.mat','Mx','My','dx','dy');
        My=My-2;
        
        dx=2e3;
        dy=2e3;
        
        dx=subsamp_x*dx;
        dy=subsamp_y*dy;
        %
        % dx=4e3;
        % dy=4e3;
        
        n=Mx/subsamp_x;
        m=(My)/subsamp_y;
        
        % m=(My-2)/subsamp_y;
        x= dx*(0:n-1);
        y= dy*(0:m-1);
        clear Mx My
        model.grid.origin=[0 0];
        model.grid.x_ref=x;
        model.grid.y_ref=y;
        [x,y]=ndgrid(x,y);
        model.grid.dX=[dx dy];
        MX=[n m];
        model.grid.MX=MX;
        
        XP=stk_dataframe ([x(:) y(:)]);


        sst_interp=fct_interp_sst(model,XP);
        model_save=model;
        save('kriging_sst_file','model_save','sst_interp');
        
    else
        model_ref=model;
        load('kriging_sst_file.mat');
        % model.l_pixel_nugget=dx_under_trace/10;
        model_ref.sst.kriging=model_save.sst.kriging;
        model_ref.sst.obs.x=model_save.sst.obs.x;
        model_ref.sst.obs.T=model_save.sst.obs.T;
        model=model_ref;
        if any(model_save.grid.MX ~= model.grid.MX)
            sst_interp=fct_interp_sst(model,XP);
            model_save.grid.MX= model.grid.MX;
            save('kriging_sst_file','model_save','sst_interp');
        end
    end
else
%     model.grid.MX = ceil(model.grid.MX/5);
%     model.grid.MX(2)=model.grid.MX(2)/2;
%     
%     n_detail=4;
% %     n_detail=2;
%     model.grid.MX=n_detail*model.grid.MX;
%     model.grid.dX=model.grid.dX/n_detail;
    
    model.grid.MX= 2^8*[1 1];
%     model.grid.MX= 2^9*[1 1];
%     model.grid.MX= 2^7*[1 1];
    model.grid.dX=2e3*[1 1];
    
    switch type_data
        case 'toy'
            sst_interp = fct_sst_toy(model);
        case 'toy2'
            sst_interp = fct_sst_toy2(model);
        case 'toy3'
            sst_interp = fct_sst_toy3(model);
        case 'toy4'
            sst_interp = fct_sst_toy4(model);
        case 'toy5'
            sst_interp = fct_sst_toy5(model);
        otherwise
            error('unkown data')
    end
    
    x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
    y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
    model.grid.x_ref=x;
    model.grid.y_ref=y;
    
%     x=model.grid.dX(1)*(0:(model.grid.MX(1)-1));
%     y=model.grid.dX(2)*(0:(model.grid.MX(2)-1));
%     imagesc(x,y,sst_interp');axis xy;
%     hold off;
%     keyboard;
end

if model.centered
    sst_interp = sst_interp - mean(sst_interp(:));    
end

if model.mirror
    sst_interp = reshape( sst_interp ,model.grid.MX);
    sst_interp = [sst_interp sst_interp(:,end:-1:1) ];
else
    sst_interp = reshape( sst_interp ,model.grid.MX);
%     sst_interp = reshape( ...
%         fct_unity_approx3(XP,model.kriging.BOX)' .* ...
%         sst_interp ...
%         ,model.grid.MX);
end

fft_sst_interp = fft2(sst_interp);

%% Compute a and w
if strcmp(type_data,'klein')
    if recompute_a_w
        switch type_adv
            case 'geostrophic_adv_alone'
                w = model.physical_constant.g / model.coriolis.f0 ...
                    .* multitrans(multi_k_x( fct_grad_h(XP,model) ));
                a = 0;% nx x d x d
                w=reshape(w,[MX 2]);
            case 'geostrophic_adv_filtered'
                w = model.physical_constant.g / model.coriolis.f0 ...
                    .* multitrans( multi_k_x(fct_grad_h(XP,model) ));
                
                w=reshape(w,[MX 2]);
            case 'geostrophic_sto_homogene'
                load(['adv_model_geostrophic_adv_filtered_' num2str(n) ...
                    '_x_' num2str(m)] ,'w','w','lambda_RMS');
                
                g_grad_h = model.coriolis.f0 ...
                    * multi_k_x(w')';
                g_grad_h =reshape(g_grad_h,[MX 2]);
                
                % Choice of a
                model.l_pixel_nugget=10e3;
                Ld=exp(-model.kriging.param(2));
                model.a0 = sqrt(2) * lambda_RMS * model.l_pixel_nugget^2;
                w = inv_J_tau_homogene_fft(model,g_grad_h);
                
            case 'geostrophic_sto'
                w = inv_J_tau_FD(model,XP); % nx x d
                w =reshape(w,[MX 2]);
                a = fct_a(XP,model);% nx x d x d
                a =reshape(a,[MX 2 2]);
                div_a(:,:,1) = divergence(x',y',a(:,:,1,1)',a(:,:,2,1)')';% n m
                div_a(:,:,2) = divergence(x',y',a(:,:,1,2)',a(:,:,2,2)')';% n m
                w = w - div_a/2; % MX x d
            otherwise
                error('wrong case');
        end
        
        %% Strain tensor
        
        [dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
        [dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
        s=(dxUx+dyUy)/2;
        d=(dyUx+dxUy)/2;
        clear s ddxUx dxUy dyUx dyUy
        lambda = sqrt(s.^2+d.^2);
        
        max(lambda(:))
        lambda_RMS = sqrt(mean(lambda(:).^2))
        
        % keyboard;
        
        %% Save
        model_save=model;
        save(['adv_model_' type_adv '_' num2str(n) '_x_' num2str(m)], ...
            'model_save','w','w','a');
    else
        
        type_adv=model.advection.type_adv;
        type_advref=type_adv;
        if strcmp(type_adv,'geostrophic_adv_lagrangian_alone') ...
                || strcmp(type_adv,'geostrophic_adv_lagrangian_sto')
            type_adv='geostrophic_adv_alone';
        elseif strcmp(type_adv,'UQ') || strcmp(type_adv,'UQ2')
            type_adv='geostrophic_adv_filtered';
        end
        load(['adv_model_' type_adv '_' num2str(n) '_x_' num2str(m)], ...
            'model_save','w','w','a','lambda_RMS');
        %     load(['adv_model_' num2str(n) '_x_' num2str(m)], ...
        %         'model_save','fft_filter','w','a');
        model.advection=model_save.advection;
        model.advection.type_adv=type_advref;
        %     model.advection.type_adv=type_adv;
    end
    w = reshape(w ,[model.grid.MX 2]);
    
else
    switch type_w
        case 'toy'
            w = fct_w_toy(model);
        case 'toy2'
            w = fct_w_toy2(model); % Pure folding
        case 'toy4'
            w = fct_w_toy4(model); % random field
        case 'toy5'
            w = fct_w_toy5(model); % Slightly squared vortex
        case 'toy6'
            w = fct_w_toy6(model); % Ellipsoid
        case 'toy7'
            w = fct_w_toy7(model); % Pure stretching
        case 'toy8'
            w = fct_w_toy8(model); % better ellispsoid
        case 'toy9'
            w = fct_w_toy9(model); % stretching + ellispsoid
        otherwise
            error('unknown type of data');
    end
%     w=w;
end
% warning('w or w');

% mirror
if model.mirror
%     w = reshape(w ,[model.grid.MX 2]);
    w_star1=[w(:,:,1) w(:,end:-1:1,1)];
    w_star2=[ w(:,:,2) -w(:,end:-1:1,2)];
    clear w
    w(:,:,1)=w_star1; clear w_star1
    w(:,:,2)=w_star2; clear w_star1
    model.grid.MX(2) = 2*model.grid.MX(2) ;
end

%% Choice of the variance tensor a
% average Lyapunov
if model.mirror
    %     My=My/2;
    model.grid.y_ref = [model.grid.y_ref ...
        model.grid.y_ref(end)+model.grid.dX(2)+model.grid.y_ref];
end
[dxUx,dyUx] = gradient(w(:,:,1)',model.grid.x_ref,model.grid.y_ref);
[dxUy,dyUy] = gradient(w(:,:,2)',model.grid.x_ref,model.grid.y_ref);
s=(dxUx+dyUy)/2;
d=(dyUx+dxUy)/2;
lambda = sqrt(s.^2+d.^2);
model.advection.lambda_RMS = sqrt(mean(lambda(:).^2));
clear s ddxUx dxUy dyUx dyUy lambda

model.advection.choice_a=choice_a;
switch model.advection.choice_a
    case 'averaged_lyapunov'
        % Lyapunov criterion for diffusion
        model.a0 = sqrt(2) * model.advection.lambda_RMS * model.l_pixel_nugget^2;
        
    case 'spectrum_geo'
        % With missed variance in the velocity spectrum
        [~, ~, missed_var_small_scale_spectrum ] ...
            = fct_sigma_spectrum(model,fft2(reshape(w,[model.grid.MX 2])));
        model.a0_on_dt = missed_var_small_scale_spectrum;
        
    case 'kriging'
        % Nugget of the Kriging
        model.l_pixel_nugget=10e3;
        model.nugget = exp(model.kriging.lognoisevariance);
        % model.nugget = (1e-2)^2;
        model.a0_on_dt =   (model.physical_constant.g/model.coriolis.f0).^2 *...
            ( model.nugget / (model.l_pixel_nugget/2)^2 );
    case 'var_w'
        model.a0_on_dt =  max(w(:).^2);
end

%% Hyperviscosity
model.advection.HV.bool=bool_HV;
if model.advection.HV.bool
    model.advection.HV.order=8;
    %     model.advection.HV.order=4;
    model.advection.HV.val= ...
        40 * model.advection.lambda_RMS * ...
        (mean(model.grid.dX)/pi)^model.advection.HV.order;
end

%% Choice of time step
% w=w;

dX2=(model.grid.dX /pi).^2;
% dX2= model.grid.dX .^2;
if strcmp(model.advection.type_adv, 'geostrophic_adv_alone') || ...
        strcmp(model.advection.type_adv, 'geostrophic_adv_lagrangian_alone')
    bound1=inf;
else
    if strcmp(model.advection.choice_a, 'averaged_lyapunov')
        bound1=2/model.a0*prod(dX2)/sum(dX2);
    else
        bound1=sqrt(2/model.a0_on_dt*prod(dX2)/sum(dX2));
    end
end
dX=permute(model.grid.dX,[1 3 2]);
bound2=sum(bsxfun(@times,abs(w),pi./dX),3);
% bound2=sum(bsxfun(@times,abs(w),1./dX),3);
bound2=max(bound2(:));
bound2=1/bound2/2;
% bound4=sum(sqrt(model.a0_on_dt)*pi./dX,3);
% % bound4=sum(sqrt(model.a0_on_dt)./dX,3);
% bound4=1/bound4/2;
if  strcmp(model.advection.type_adv, 'geostrophic_adv_lagrangian_alone')
    bound3=inf;
else
    if ~model.advection.HV.bool
        bound3=inf;
    else
        bound3=1/model.advection.HV.val*(prod(dX2)/sum(dX2)) ^ ...
            (model.advection.HV.order/2);
    end
end
bound_dt = min([bound1 bound2 bound3]);
clear dX dX2
model.advection.dt_adv = bound_dt;
if strcmp(model.advection.choice_a, 'averaged_lyapunov')
    model.a0_on_dt = model.a0 / model.advection.dt_adv;
else
    model.a0=model.a0_on_dt*model.advection.dt_adv;
end

model.advection.coef_diff =  1/2 * model.a0;

%% Advection
model.advection.SQG=param_SQG;
model.type_data=type_data;
fprintf('advection \n')
[fft_T_adv1,model] = fct_lagrangian_advection(model, fft_sst_interp, w );
save('res_temp1','fft_T_adv1')
keyboard;