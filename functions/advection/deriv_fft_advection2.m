function d_fft_T_adv = deriv_fft_advection2(model, fft_T, w)
% function d_fft_T_adv = deriv_fft_advection(model, fft_T, fft_w)
% Advection of T with the speed w, in Fourrier space
%

meth_anti_alias=model.advection.meth_anti_alias;
if model.advection.evol_variance
    fft_Var_T=fft_T(:,:,2);
    fft_T=fft_T(:,:,1);
end

%% Grid
M=prod(model.grid.MX);
if any( mod(model.grid.MX,2)~=0)
    error('the number of grid points by axis need to be even');
end
PX=model.grid.MX/2;
if strcmp(meth_anti_alias,'deriv_LowPass')
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ...
        .* fct_unity_approx5(model.grid.MX(1));
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1]...
        .* fct_unity_approx5(model.grid.MX(2));
else
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
end
[kx,ky]=ndgrid(kx,ky);

% %% A supprimer
% fft_T=1./sqrt(kx.^2+ky.^2);
% fft_T(1,1)=0;

%% Advection term

adv1x = -2*1i*pi/model.grid.dX(1)* kx .* fft_T;
adv1y = -2*1i*pi/model.grid.dX(2)* ky .* fft_T;
% adv1x = 2*1i*pi/model.grid.dX(1)* kx .* fft_T;
% adv1y = 2*1i*pi/model.grid.dX(2)* ky .* fft_T;
gradT(:,:,1)=real(ifft2(adv1x));
gradT(:,:,2)=real(ifft2(adv1y));
wgradT=sum(bsxfun(@times,w,gradT),3);
adv1=fft2(wgradT);clear wgradT
adv1(PX(1)+1,:)=0;
adv1(:,PX(2)+1)=0;

if model.advection.evol_variance
    adv1_Varx = -2*1i*pi/model.grid.dX(1)* kx .* fft_Var_T;
    adv1_Vary = -2*1i*pi/model.grid.dX(2)* ky .* fft_Var_T;
%     adv1_Varx = 2*1i*pi/model.grid.dX(1)* kx .* fft_Var_T;
%     adv1_Vary = 2*1i*pi/model.grid.dX(2)* ky .* fft_Var_T;
    gradT_Var(:,:,1)=real(ifft2(adv1_Varx));
    gradT_Var(:,:,2)=real(ifft2(adv1_Vary));
    wgrad_VarT=sum(bsxfun(@times,w,gradT_Var),3);
    adv1_Var=fft2(wgrad_VarT);clear wgrad_VarT
    
    %% Source term
    % if model.advection.evol_variance
    adv3_Var = model.a0 * sum(gradT.^2,3);
%     adv3_Var = 2*model.advection.coef_diff * sum(gradT.^2,3);
    adv3_Var = fft2(adv3_Var);
    adv3_Var(PX(1)+1,:)=0;
    adv3_Var(:,PX(2)+1)=0;
%     % else
%     %     adv3=0;
%     if ~strcmp(meth_anti_alias,'fct_LowPass')
%         adv3_Var= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',adv3_Var);
%         adv3_Var= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,adv3_Var);
%     end
end
clear gradT


%% Diffusion
if strcmp(meth_anti_alias,'Lax_Wendroff') || ...
        strcmp(model.advection.type_adv,'geostrophic_adv_filtered') || ...
        strcmp(model.advection.type_adv,'geostrophic_sto_homogene')|| ...
        ( strcmp(model.advection.type_adv,'UQ') && ...
        strcmp(model.advection.step,'finite_variation') ) || ...
        ( strcmp(model.advection.type_adv,'UQ2') && ...
        strcmp(model.advection.step,'finite_variation') )
    
    kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
    ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
    [kx,ky]=ndgrid(kx,ky);
    
    adv2 = model.advection.coef_diff * ...
        - (2*pi)^2 * ( (1/model.grid.dX(1)* kx).^2 + (1/model.grid.dX(2)* ky).^2 ) ...
        .* fft_T;
    
    if model.advection.evol_variance
        adv2_Var = model.advection.coef_diff * ...
            - (2*pi)^2 * ( (1/model.grid.dX(1)* kx).^2 + (1/model.grid.dX(2)* ky).^2 ) ...
            .* fft_Var_T;
    end
else
    adv2=0;
    adv2_Var = 0;
end

%% Hyperviscosity
if model.advection.HV.bool 
% if model.advection.HV.bool && ...
%          ( ( strcmp(model.advection.type_adv,'UQ') && ...
%              strcmp(model.advection.step,'finite_variation') ) || ...
%            ( strcmp(model.advection.type_adv,'UQ2') && ...
%              strcmp(model.advection.step,'finite_variation') ) )
    
%         kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
%         ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
%         [kx,ky]=ndgrid(kx,ky);
    
    adv4 = model.advection.HV.val * ...
        - ( (2*pi)^2 * ( (1/model.grid.dX(1)* kx).^2 + (1/model.grid.dX(2)* ky).^2 ) ) ...
        .^ (model.advection.HV.order/2) ...
        .* fft_T;
    
    if model.advection.evol_variance
        adv4_Var = model.advection.HV.val * ...
            - ( (2*pi)^2 * ( (1/model.grid.dX(1)* kx).^2 + (1/model.grid.dX(2)* ky).^2 ) ) ...
            .^ (model.advection.HV.order/2) ...
            .* fft_Var_T;
    end
else
    adv4=0;
    adv4_Var = 0;
end


%% Summing


% figure;fct_spectrum(model,fft_T);
% keyboard;
% figure;fct_spectrum(model,adv1*model.advection.dt_adv);
% hold on;
% fct_spectrum(model,adv2*model.advection.dt_adv);
% fct_spectrum(model,adv4*model.advection.dt_adv);
% hold off;
% keyboard;

d_fft_T_adv=adv1+adv2+adv4; clear adv1 adv2 adv4
% d_fft_T_adv=adv1+adv2; clear adv1 adv2 adv3
% % d_fft_T_adv=adv1x+adv1y+adv2; clear adv1x adv1y adv2
if model.advection.evol_variance
    d_fft_Var_T_adv=adv1_Var+adv2_Var+adv3_Var+adv4_Var;
    clear adv1_Var adv2_Var adv3_Var adv4_Var
    %     d_fft_Var_T_adv=adv1_Var+adv2_Var+adv3_Var; clear adv1_Var adv2_Var adv3_Var
end

%% Anti-aliasing filter
if strcmp(meth_anti_alias,'fct_LowPass')
    d_fft_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',d_fft_T_adv);
    d_fft_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,d_fft_T_adv);
    
    if model.advection.evol_variance
        d_fft_Var_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',d_fft_Var_T_adv);
        d_fft_Var_T_adv= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,d_fft_Var_T_adv);
    end
end

%%

if model.advection.evol_variance
    d_fft_T_adv(:,:,2)=d_fft_Var_T_adv;
end

