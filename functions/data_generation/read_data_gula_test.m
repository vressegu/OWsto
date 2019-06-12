function read_data_gula_test()
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

init


model.type_data =  'Gula';
remove_cost = true;
init_model

% day = 30;
        day = 60;
time_day = 'day';
        
        
vectname_data = {'u','v','zeta'};

x= 900e3/1200 * (0:(1200-1));
y= 1050e3/1400 * (0:(1400-1));

% Remove cost
if remove_cost
    xx= x(501:end);
    yy= y(1:900);
else
    xx= x;
    yy= y;
end

model.grid.x = xx;
model.grid.y = yy;
model.grid.origin = [xx(1) yy(2)];
model.grid.dX = [xx(2)-xx(1) yy(2)-yy(1)];
model.grid.MX = [length(xx) length(yy)];
model.grid.BOX = [ xx(1) xx(end)  ; yy(1) yy(end) ];

mask_boundaries = fct_unity_approx6(model.grid.MX(1))' * ...
    fct_unity_approx6(model.grid.MX(2));

[xx,yy]=ndgrid(xx,yy);

t_ini = 380;
t_final = 380;
% t_final = 738;

t2 = 380 + 2*(day-1);

w = nan( [model.grid.MX length(vectname_data)]);

for k=1:length(vectname_data)
    
    name_data=vectname_data{k};
    fprintf([name_data ' \n'])
    ncdisp([ ...
        '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
        '/CHABU_surface/chabu_surf.0' num2str(t2) '.nc'] );
    w_temp=ncread([ ...
        '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
        '/CHABU_surface/chabu_surf.0' num2str(t2) '.nc'],name_data );
    w_temp=w_temp(2:1201,2:1401,:,:);
    
    
    %     % Remove cost
    if remove_cost
        w_temp=w_temp(501:end,1:900,:,:);
    end
    
    
    if strcmp(time_day,'night')
        w(:,:,k)=w_temp(:,:,1);
    elseif strcmp(time_day,'day')
        w(:,:,k)=w_temp(:,:,2);
    else
        warning('unknown day time');
    end
    
    w_temp=w(:,:,k);
    
%     % Remove boundaries
%     w_temp=w_temp .* mask_boundaries;
    
    w_temp= w_temp-mean(w_temp(:));
%     close all
    fct_plot_f(model,double(w_temp));
    
    if k==3
        % Physical parameters
        model = fct_physical_param(model);
        
        % Coriolis
        nabla_ssh = gradient_mat_2( w_temp ,model.grid.dX);
        w(:,:,1,:)= - nabla_ssh(:,:,2,:);
        w(:,:,2,:)= nabla_ssh(:,:,1,:);
        clear nabla_ssh;
        w = ...
            model.physical_constant.g / model.physical_constant.f0 ...
            * w;
        
        fct_plot_f(model,double(w(:,:,1,:)));
        fct_plot_f(model,double(w(:,:,2,:)));
    end
    
end

end
