function [w,model] = read_data_gula(day,time_day,model)
% Read the velocity of the day (day or night) of the simulation of Gula,
% Molemaker, McWilliams
%

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
t_final = 738;

t2 = 380 + 2*(day-1);

w = nan( [model.grid.MX 2]);

for k=1:length(vectname_data)
    
    name_data=vectname_data{k};
%     fprintf([name_data ' \n'])
    w_temp=ncread([ ...
        '/Volumes/WD_Ressegui/These/from_tour/data_brest/data_gula' ...
        '/CHABU_surface/chabu_surf.0' num2str(t2) '.nc'],name_data );
    w_temp=w_temp(2:1201,2:1401,:,:);
    
    % Remove cost
    w_temp=w_temp(501:end,1:1100,:,:);
    
    if strcmp(time_day,'night')
        w(:,:,k)=w_temp(:,:,:,1);
    elseif strcmp(time_day,'day')
        w(:,:,k)=w_temp(:,:,:,2);
    else
        warning('unknown day time');
    end
end

end
