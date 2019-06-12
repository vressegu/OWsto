function [salt,temp,model] = read_data_gula_tracer_multitime()
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

model.type_data = 'Gula';
model.filtering.smoothing= true;
init_model
vectname_data = {'u','v','salt','temp'};

Mx_LR = 64; % filter of 4 km -> resolution ~ 10km
% Mx_LR = 8; % filter of 30 km -> resolution 50-100km

x= 900e3/1200 * (0:(1200-1));
y= 1050e3/1400 * (0:(1400-1));
% n_y = 1100;
n_y = 900;

% Remove cost
xx= x(501:end);
yy= y(1:n_y);

% xx= x;
% yy= y;

model.grid.x = xx;
model.grid.y = yy;
model.grid.origin = [xx(1) yy(2)];
model.grid.dX = [xx(2)-xx(1) yy(2)-yy(1)];
model.grid.MX = [length(xx) length(yy)];
model.grid.BOX = [ xx(1) xx(end)  ; yy(1) yy(end) ];

[xx,yy]=ndgrid(xx,yy);

Mx_HR = min(model.grid.MX);
smooth_rate = Mx_HR/Mx_LR

t_ini = 380;
% t_final = 380;
t_final = 61*2+380;
% t_final = 738;
N_t = t_final - t_ini +1;
time = t_ini:t_final;
one_day = 24 * 3600;
one_week = 7 * one_day;
model.dt = 0.5 * one_day;

% t2 = 380 + 2*(day-1);

% temp = nan( [model.grid.MX 1 N_t+1]);
% salt = nan( [model.grid.MX 1 N_t+1]);
w = nan( [model.grid.MX 2 N_t+1]);
w_HR = nan( [model.grid.MX 2 N_t+1]);

for k=1:length(vectname_data)
    data_tot = nan( [model.grid.MX 1 N_t+1]);
    data_tot_smooth = nan( [model.grid.MX 1 N_t+1]);
    
    %% Load data
    for k_t = 1:2:N_t
        t_local = time(k_t);
        
        name_data=vectname_data{k};
        %     fprintf([name_data ' \n'])
        data_temp=ncread([ ...
            '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
            '/CHABU_surface/chabu_surf.0' num2str(t_local) '.nc'],name_data );
        data_temp=data_temp(2:1201,2:1401,:,:);
        
        % Remove cost
        data_temp=data_temp(501:end,1:n_y,:,:);
        
        data_temp=double(data_temp);
        
        
        data_tot(:,:,1,[k_t k_t+1])=data_temp;
    end
    eval([ name_data '_HR=data_tot;']);
    
    %% Smoothing
    if model.filtering.smoothing
        for k_t = 1:1:(N_t+1)
            data_temp=data_tot(:,:,1,k_t);
            
            if  Mx_LR < Mx_HR
                % Filter
                sigma_filter = min(model.grid.dX)/2 * Mx_HR/Mx_LR ;
                model.advection.sigma_filter=sigma_filter;
                filter_perso = design_filter_3(model);
                
                % Convolution
                data_temp_smooth = conv2(data_temp,filter_perso,'same');
            else
                warning('the resolution is already too coarse');
            end
            data_tot_smooth(:,:,1,k_t)=data_temp_smooth;
        end
    end
    eval([ name_data '=data_tot_smooth;']);
    
    % 2D Velocity
    if k==2
        w(:,:,1,:)=u;
        w(:,:,2,:)=v;
        clear u v;
        w_HR(:,:,1,:)=u_HR;
        w_HR(:,:,2,:)=v_HR;
        clear u_HR v_HR;
    end
    
    %% Compute growth rate HR
    if k > 2
        std_tracer = std(data_tot(:));
        DTracerDt = material_deriv_mat(...
            data_tot, w_HR, model.dt, model.grid.dX);
        mean_source_HR = mean( abs(DTracerDt(:)) )/(std_tracer/one_week)
        
%         eval([ name_data '=data_tot;']);
        
        grad_tracer = gradient_mat_2(data_tot,model.grid.dX); % Mx My 2 N_t
        log_n_grad = 1/2 * log( sum(grad_tracer.^2,3) );
        growth_rate = material_deriv_mat(...
            log_n_grad, w_HR, model.dt, model.grid.dX);
        
        eval([ 'growth_rate_' name_data '_HR=growth_rate;']);
        eval([ 'source_' name_data '_HR=DTracerDt;']);
    
        %     eval([ name_data '(:,:,k,[k_t k_t+1])=data_temp;']);
    end
    
    %% Compute growth rate LR
    if k > 2
        std_tracer = std(data_tot_smooth(:));
        DTracerDt = material_deriv_mat(...
            data_tot_smooth, w, model.dt, model.grid.dX);
        
        marge = 20;
        mean_source_LR = DTracerDt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,:);
        mean_source_LR = mean( abs(mean_source_LR(:)) )/(std_tracer/one_week)
        
%         eval([ name_data '=data_tot;']);
        
        grad_tracer = gradient_mat_2(data_tot_smooth...
            ,model.grid.dX); % Mx My 2 N_t
        log_n_grad = 1/2 * log( sum(grad_tracer.^2,3) );
        growth_rate = material_deriv_mat(...
            log_n_grad, w, model.dt, model.grid.dX);
        
        eval([ 'growth_rate_' name_data '=growth_rate;']);
        eval([ 'source_' name_data '=DTracerDt;']);
    
        %     eval([ name_data '(:,:,k,[k_t k_t+1])=data_temp;']);
    end
end


%% Save
save([model.folder.folder_simu '/tracer_multi_time.mat'],...
    'temp_HR','salt_HR','w','w_HR',...
    'source_temp','source_salt',...
    'temp','salt',...
    'growth_rate_temp','growth_rate_salt',...
    'growth_rate_temp_HR','growth_rate_salt_HR',...
    'model','-v7.3');


%% Plots tracer
N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t))
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t))
    drawnow;
end


%% Plots of Lagrangian source of tracer : are they really tracers?
N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

std_salt = salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,:);
std_salt = std(std_salt(:));

cax = [-1 1]*20;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,source_salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t)/(std_salt/one_week),...
        cax)
    %     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end


std_temp = temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,:);
std_temp = std(std_temp(:));

for k_t = 1:(N_t+1)
    fct_plot_f(model_,source_temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t)/(std_temp/one_week),...
        cax)
    drawnow;
end

%% Plots tracer gradient growth rate
N_t = size(salt,4)-1;
marge = 20;
model_ = model;
model_.grid.MX = model_.grid.MX - 2*marge;

cax = [-1 1]*100;
for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_salt(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
%     fct_plot_f(model,salt(:,:,:,k_t))
    drawnow;
end

for k_t = 1:(N_t+1)
    fct_plot_f(model_,one_week * growth_rate_temp(...
        (1+marge):end-marge,(1+marge):end-marge,...
        :,k_t),...
        cax)
    drawnow;
end


keyboard;

end
