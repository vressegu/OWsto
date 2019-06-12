function w = fct_w_toy6(model)
% Create an incompressible velocity with circle streamlines and an angular
% velocity which is Gaussian function of the distance at the origin
%

%%
R=0.5;
c=0.1;
% R=1;
% c=0.2;

center1=c*[1 1];
R1=R;
vort1 = fct_vortex(model,center1,R1);

center1=-c*[1 1];
R1=R;
vort2 = fct_vortex(model,center1,R1);

vort = vort1+vort2;

w = vort2velocity(model.grid, vort);

% model.odg_b=1e4;
% model.k_inf_on_k1 = 4;
% model.slope_b_ini = -8;
% % model.slope_b_ini = -5;
% phi = init_Spectrum(model);
% w = gradient_ortho_perso(model.grid, phi);


%% Plots
% sqrt(mean(w(:).^2))
vort = vorticity_perso(model.grid, w);

width=3.5;
height=3;
X0 = [0 0];

% figure(2)
figure('Name','Angular velocity','NumberTitle','off','Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% x=x-mean(x);
y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% y=y-mean(y);
imagesc(x,y,vort');
axis xy; axis equal;
% title('Velocity norm');
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('y(m)',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('x(m)',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
colorbar
title('$Vorticity$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'interpreter','latex',...
    'FontName','Times')
axis([x(1) x(end) y(1) y(end)])

% Create the folders
fct_create_folder_plots(model)

eval( ['print -depsc ' model.folder.folder_simu '/vort_w.eps']);

end


function vort = fct_vortex(model_,center,Radius)

% LX=model_.grid.MX .* model_.grid.dX;

vort = fct_vortex_(model_,center,Radius);
vort = vort + fct_vortex_(model_,center+[2 0],Radius);
vort = vort + fct_vortex_(model_,center-[2 0],Radius);
vort = vort + fct_vortex_(model_,center+[2 2],Radius);
vort = vort + fct_vortex_(model_,center-[2 2],Radius);
vort = vort + fct_vortex_(model_,center+[2 -2],Radius);
vort = vort + fct_vortex_(model_,center+[-2 +2],Radius);
vort = vort + fct_vortex_(model_,center+[0 2],Radius);
vort = vort + fct_vortex_(model_,center-[0 2],Radius);

end


function phi = fct_vortex_(model_,center,Radius)
% function phi = fct_vortex(model_,center,Radius)

LX = model_.grid.dX .* model_.grid.MX;
x = model_.grid.dX(1)*(0:model_.grid.MX(1)-1);
x=x-mean(x);
x=x-LX(1)/2*center(1);
y = model_.grid.dX(2)*(0:model_.grid.MX(2)-1);
y=y-mean(y);
y=y-LX(2)/2*center(2);
r = sqrt(bsxfun(@plus,(x').^2,(y).^2));

sigma = 0.1*sqrt(prod(model_.grid.dX.*model_.grid.MX)); % ~ Rossby radius
radius=0.03*sqrt(prod(model_.grid.dX.*model_.grid.MX));
radius= Radius * radius;

ampli = model_.coriolis.f0/2;
ampli = ampli/8;
theta_dot = ampli * exp(-1/2*(r-radius).^2/sigma^2);

siz = size(theta_dot);
% r=r(:);
theta_dot = theta_dot(:);
theta_dot(r<radius) = ampli ;
phi = reshape(theta_dot,siz);
end

