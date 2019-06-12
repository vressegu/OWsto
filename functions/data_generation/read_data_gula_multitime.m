function [w,model] = read_data_gula_multitime()
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

model.type_data = 'Gula';
init_model
vectname_data = {'u','v'};

x= 900e3/1200 * (0:(1200-1));
y= 1050e3/1400 * (0:(1400-1));

% Remove cost
xx= x(501:end);
yy= y(1:1100);

% xx= x;
% yy= y;

model.grid.x = xx;
model.grid.y = yy;
model.grid.origin = [xx(1) yy(2)];
model.grid.dX = [xx(2)-xx(1) yy(2)-yy(1)];
model.grid.MX = [length(xx) length(yy)];
model.grid.BOX = [ xx(1) xx(end)  ; yy(1) yy(end) ];

[xx,yy]=ndgrid(xx,yy);

t_ini = 380;
% t_final = 380;
t_final = 61*2+380;
% t_final = 738;
N_t = t_final - t_ini +1;
time = t_ini:t_final;

% t2 = 380 + 2*(day-1);

w = nan( [model.grid.MX 2 N_t+1]);
% w = nan( [model.grid.MX 2]);

for k_t = 1:2:N_t
    t_local = time(k_t);
    for k=1:length(vectname_data)
        
        name_data=vectname_data{k};
        %     fprintf([name_data ' \n'])
        w_temp=ncread([ ...
            '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
            '/CHABU_surface/chabu_surf.0' num2str(t_local) '.nc'],name_data );
        w_temp=w_temp(2:1201,2:1401,:,:);
        
        % Remove cost
        w_temp=w_temp(501:end,1:1100,:,:);
        
        %     if strcmp(time_day,'night')
        %         w(:,:,k)=w_temp(:,:,:,1);
        %     elseif strcmp(time_day,'day')
        %         w(:,:,k)=w_temp(:,:,:,2);
        %     else
        w(:,:,k,[k_t k_t+1])=w_temp;
    end
end


save([model.folder.folder_simu '/w_HR_multi_time.mat'],'w','-v7.3');


for k_t = 1:N_t
    fct_plot_velocity_2(model,w(:,:,:,k_t))
    drawnow;
end
keyboard;

end
