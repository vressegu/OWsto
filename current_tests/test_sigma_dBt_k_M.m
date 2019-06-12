%% Test of fct_a

init;

% beta_v = 0:0.5:5;
% beta_v = 2.5;
beta_v = 5/3;
kappamax_on_resolution_v = 1./(8);
% kappamax_on_resolution_v = 1./(5);
% kappamax_on_resolution_v = 1./(1:10);

k_m = 1e-5;
% % k_m = 3e-5;

for beta = beta_v(end:-1:1)
    for kappamax_on_resolution = kappamax_on_resolution_v
        beta
        kappamax_on_resolution
        
        test_mean_nrj = false;
        test_val_a = false;
        test_grad_sigma = true;
        
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
        % alpha_adv_back=true;
        
        % centered = false;
        
        % Spectrum slope of sigma dBt
        slope_sigma = - beta;
        % slope_sigma = - 5/3;
        
        % Rate between the smallest and the largest wave number of sigma dBt
        % kappamin_on_kappamax = 1/50;
        
        % The largest wave number of sigma dBt
        % kappamax_on_resolution = 1/3;
        
        V0 = 1;
        
        dt =1;
        % a0=10;
        
        %% Grid
        
        N_ech=500;
        
        model.grid.MX= 2^8*[1 1];
        model.grid.dX=2e3*[1 1];
        
        % Wave numbers
        M_kappa=min(model.grid.MX);
        P_kappa= M_kappa/2;
        kappa=1/(M_kappa)* (0:(P_kappa-1)) ;
        kappa= 2*pi*max(1./model.grid.dX) * kappa;
        
        % Wave vector
        PX=model.grid.MX/2;
        kx=1/(model.grid.MX(1))*[ 0:(PX(1)-1) 0 (1-PX(1)):-1] ;
        ky=1/(model.grid.MX(2))*[ 0:(PX(2)-1) 0 (1-PX(2)):-1];
        kx=2*pi/model.grid.dX(1)*kx;
        ky=2*pi/model.grid.dX(2)*ky;
        [kx,ky]=ndgrid(kx,ky);
        k=sqrt(kx.^2+ky.^2);
        k(PX(1)+1,:)=inf;
        k(:,PX(2)+1)=inf;
        
        % Spatial grid
        x = model.grid.dX(1) * ( 0:(model.grid.MX(1)-1) );
        y = model.grid.dX(2) * ( 0:(model.grid.MX(2)-1) );
        [x,y]=ndgrid(x,y);
        norm_x = sqrt(x.^2+y.^2);
        
        %% Parameters
        
        % model = fct_w_a_geostrophic_stochastic_kriging(recompute_kriging_ssh);
        
        model.sigma.V0 = V0;
        % model.plot.alpha_adv_back = alpha_adv_back;
        % model.mirror=mirror;clear mirror
        % model.plot.filtered_field=plot_filtered_field;clear plot_filtered_field
        % model.advection.type_adv=type_adv;
        % model.centered= centered;
        % model.sigma.k_c = k_c;
        model.physical_constant.f0 = 1e-4;
        % model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
        model.sigma.kappamax_on_resolution = kappamax_on_resolution;
        model.sigma.slope_sigma = slope_sigma;
        % model.folder.folder_simu = [ 'images/SQG_MU/' type_data ];
        model.folder.folder_simu = [ 'images/SQG_MU/T_' type_data '_w_' type_w ];
        
        % model.sigma.a0 = a0;
        model.advection.dt_adv = dt;
        
        % Largest wave number
        k_M = kappamax_on_resolution * kappa(P_kappa);
        % k_M = kappa(P_kappa);
        
        % Smallest wave number
        kappamin_on_kappamax = k_m / k_M;
        % k_M / k_m
        model.sigma.kappamin_on_kappamax = kappamin_on_kappamax;
        % k_m = model.sigma.kappamin_on_kappamax * k_M;
        
        beta = - slope_sigma;
        
        %% Covariance
        cov_psi = V0 * 2^((1-beta)/2)/gamma((beta+1)/2) * ...
            ( k_m * norm_x ).^((beta+1)/2) .* ...
            besselk( (beta+1)/2 , k_m * norm_x );
        cov_psi(1,1) = V0;
        
        cov_psi_period = [ cov_psi(end-1:-1:2,:); cov_psi ];
        cov_psi_period = [ cov_psi_period(:,end-1:-1:2) cov_psi_period ];
        grid_modif = model.grid;
        grid_modif.MX = 2 * grid_modif.MX-2;
        
        % figure;imagesc(cov_psi_period');axis xy;axis equal;colorbar;
        
        cov_a_period = gradient_ortho_perso(grid_modif, cov_psi_period);
        % figure;imagesc(cov_a(:,:,1,1)');axis xy;axis equal;colorbar;
        cov_a_period = permute(cov_a_period,[1 2 4 3]);
        % figure;imagesc(cov_a(:,:,1,1)');axis xy;axis equal;colorbar;
        cov_a_period = - gradient_ortho_perso(grid_modif, cov_a_period);
        
        % figure;imagesc(cov_a(:,:,1,1)');axis xy;axis equal;colorbar;
        
        cov_a = cov_a_period( (model.grid.MX(1))-1:2*model.grid.MX(1)-2, ...
            (model.grid.MX(2))-1:2*model.grid.MX(2)-2,:,:);
        
        % figure;imagesc(cov_a(:,:,1,1)');axis xy;axis equal;colorbar;
        
        if test_val_a
            a_0_from_cov = cov_a(1,1,1,1)
        end
        
        nab_cov_a_period = reshape(cov_a_period,[ grid_modif.MX 1 4]);
        nab_cov_a_period = gradient_perso(grid_modif, nab_cov_a_period);
        nab_cov_a_period = reshape(nab_cov_a_period,[ grid_modif.MX 1 8]);
        nab_cov_a_period = - gradient_perso(grid_modif, nab_cov_a_period);
        nab_cov_a_period = reshape(nab_cov_a_period,[ grid_modif.MX 2 2 2 2]);
        
        nab_cov_a = nab_cov_a_period( (model.grid.MX(1))-1:2*model.grid.MX(1)-2, ...
            (model.grid.MX(2))-1:2*model.grid.MX(2)-2,:,:,:,:);
        
        figure;imagesc(nab_cov_a(:,:,1,1,1,1)');axis xy;axis equal;colorbar;
        
        nab_a_0_from_cov = squeeze(nab_cov_a(1,1,:,:,:,:));
        clear nab_cov_a
        
        %% Spectrum
        
        V0 =model.sigma.V0;
        beta = - model.sigma.slope_sigma;
        A2 = 2*pi * (beta+1) * V0 / (k_m^2);
        A2 = 1/prod(model.grid.dX) * A2;
        % The factor 1/prod(dX) is to go from continuous to discrete Fourier
        % transform
        psi_sigma = sqrt(A2) * (1 + (1/k_m * k).^2 ).^(-(beta+3)/4);
        % Antialiasing
        psi_sigma(PX(1)+1,:)=0;
        psi_sigma(:,PX(2)+1)=0;
        
        S_psi = psi_sigma .^2;
        
        cov_psi_theo = real(ifft2(S_psi));
        figure;
        subplot(2,2,1);imagesc(cov_psi');
        axis equal; axis xy;colorbar;cax=caxis;
        subplot(2,2,2);imagesc(cov_psi_theo');
        axis equal; axis xy;colorbar;
        % caxis(cax);
        % keyboard;
        
        % Fourier transform of the kernel \tilde sigma up to a multiplicative
        % constant
        [sigma_on_sq_dt,f_sigma,a0_from_spectrum,spectrum_sigma] = ...
            fct_sigma_spectrum_Matern(model);
        clear w_ini;
        % sigma_on_sq_dt will be used to simulate sigma d B_t
        % missed_var_small_scale_spectrum will be used to set the mulitplicative
        % constant
        
        if test_val_a
            a0_from_spectrum
        end
        
        
        %% Sampling
        % Fourier transform of white noise
        dBt_C_on_sq_dt = sqrt(model.advection.dt_adv) * ...
            fft2( randn( [ model.grid.MX 1 N_ech]));
        % Multiplication by the Fourier transform of the kernel \tilde \sigma
        fft_sigma_dBt = bsxfun(@times,sigma_on_sq_dt,dBt_C_on_sq_dt);
        clear dBt_C_on_sq_dt
        % Homogeneous velocity field
        % sigma_dBt_dt = sqrt(prod(model.grid.MX)) * real(ifft2(fft_sigma_dBt));
        sigma_dBt_dt = real(ifft2(fft_sigma_dBt));
        clear fft_sigma_dBt
        
        if test_val_a
            sigma_dBt_estim = sigma_dBt_dt.^2/model.advection.dt_adv;
            % sigma_dBt_estim = sigma_dBt_dt(1,1,1,:).^2/model.advection.dt_adv;
            a0_estim = mean(sigma_dBt_estim(:))
        end
        
        sigma_dBt_estim = (sigma_dBt_dt(:,:,1,1).^2 + ...
            sigma_dBt_dt(:,:,2,1).^2)   ...
            /model.advection.dt_adv;
        if test_mean_nrj
            mean_nrj_estim = mean(sigma_dBt_estim(:))
        end
        
        nab_sigma_dBt = gradient_perso(model.grid, ...
            reshape(sigma_dBt_dt, [model.grid.MX 1 2*N_ech]));
        nab_sigma_dBt = reshape(nab_sigma_dBt ,[model.grid.MX 2 1 2 1 N_ech]);
        di_sigma_dp_sigma_from_sampl = bsxfun(@times, nab_sigma_dBt, ...
            permute( nab_sigma_dBt , [1 2 4 3 6 5 7]) );
        di_sigma_dp_sigma_from_sampl = mean(di_sigma_dp_sigma_from_sampl,7);
        di_sigma_dp_sigma_from_sampl = mean(di_sigma_dp_sigma_from_sampl,2);
        di_sigma_dp_sigma_from_sampl = mean(di_sigma_dp_sigma_from_sampl,1);
        di_sigma_dp_sigma_from_sampl = squeeze(di_sigma_dp_sigma_from_sampl);
        
        
        %% Theoritical relations
        
        a0_from_formula = V0*k_m^2/(beta-1);
        if test_val_a
            a0_from_formula
        end
        
        if beta ~= 3
            c_inf =  a0_from_formula * k_m^2 / (beta-3);
            % c_inf =  1/4 * a0_from_formula * (k_m^2 / 4) * (beta-1)^2 / (beta-3);
            c_k_M = 1 - (beta^2-1)/8 * (k_M/k_m).^(3-beta);
            c_k_M_from_formula = c_inf * c_k_M ;
        else
            c_k_M_from_formula =  a0_from_formula * k_m^2 * ...
                ( log(k_M/k_m) - 3/4 );
        end
        
        c_k_M_from_cov = squeeze(nab_a_0_from_cov(1,1,1,1));
        c_k_M_from_sampl = squeeze(di_sigma_dp_sigma_from_sampl(1,1,1,1));
        
        err = abs(c_k_M_from_formula-c_k_M_from_sampl)/abs(c_k_M_from_sampl)
%         err = abs(c_k_M_from_formula-c_k_M_from_sampl)/abs(c_k_M_from_sampl-c_inf)
        %     c_k_M_from_formula-c_inf
        %     c_k_M_from_sampl-c_inf
        
        % if test_grad_sigma
        %
        %     for i=1:2
        %         % J e_i
        %         e_i = zeros(2,1);
        %         e_i(i)=1;
        %         J_e_i(1) = - e_i(2);
        %         J_e_i(2) = + e_i(1);
        %         for p=1:2
        %             % J e_p
        %             e_p = zeros(2,1);
        %             e_p(p)=1;
        %             J_e_p(1) = - e_p(2);
        %             J_e_p(2) = + e_p(1);
        %
        %             JeiJep= J_e_i' * J_e_p;
        %             twoPi_JeiJep = JeiJep + JeiJep';
        %
        %             i
        %             p
        %
        %             di_sigma_dp_sigma_from_formula = (i==p) * eye(2) + twoPi_JeiJep;
        %             di_sigma_dp_sigma_from_formula = ...
        %                 c_k_M * di_sigma_dp_sigma_from_formula
        %
        % %             di_sigma_dp_sigma_from_cov = squeeze(nab_a_0_from_cov(i,p,:,:))
        %             di_sigma_dp_sigma_from_sampl = squeeze(nab_a_0_from_cov(i,p,:,:))
        %         end
        %     end
        
        %     lap_a_from_formula = 4 *c_k_M
        %
        %     lap_a_from_cov = squeeze(nab_a_0_from_cov(1,1,:,:))+...
        %         squeeze(nab_a_0_from_cov(2,2,:,:))
        %
        % end
        
    end
end

keyboard;