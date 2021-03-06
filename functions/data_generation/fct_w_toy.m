function w = fct_w_toy(model)
% Create an incompressible velocity with circle streamlines and an angular
% velocity which is Gaussian function of the distance at the origin
%

% w=nan(model.grid.MX);
% x = 1/model.grid.MX(1)*(0:model.grid.MX(1)-1);
% y = 1/model.grid.MX(2)*(0:model.grid.MX(2)-1);
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
x=x-mean(x);
y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
y=y-mean(y);
% [x,y]=ndgrid(x,y);
% r = sqrt(bsxfun(@plus,(x').^2,(2*y).^2));
% y=2*y;
r = sqrt(bsxfun(@plus,(x').^2,(y).^2));
% clear x y
% sigma = 2*0.05*sqrt(prod(model.grid.dX.*model.grid.MX)); % ~ Rossby radius
% radius=0;
sigma = 0.05*sqrt(prod(model.grid.dX.*model.grid.MX)); % ~ Rossby radius
radius=0.2*sqrt(prod(model.grid.dX.*model.grid.MX));
ampli = model.coriolis.f0/2;
theta_dot = ampli * exp(-1/2*(r-radius).^2/sigma^2);

% figure;imagesc(r); axis xy;axis equal; keyboard;
% figure;imagesc(theta_dot); axis equal;  keyboard;

[x,y]=ndgrid(x,y);
w(:,:,1) = -y ; w(:,:,2) = x ;
w = 1/2 * bsxfun(@times,w,theta_dot);


%% Plots
% width=3;
% height=2;
width=2.5;
height=2;
X0=[0 0];


% figure('Name','Angular velocity','NumberTitle','off','Units','inches', ...
%     'Position',[X0(1) X0(2) width height], ...
%     'PaperPositionMode','auto');


phi = -sigma^2/2 * theta_dot;

% x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% figure;imagesc(x,y,phi'); axis equal;
% keyboard;

% figure;
% x=x+x(1);
% y=y+y(1);
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% color=colormap;
% figure;
hold on;
cc=get(gca,'Clim');
% keyboard;
% contourf(x,y,-phi',5);
contour(x,y,phi',2,'k');
% % imagesc(x,y,theta_dot');
% % axis equal;
% % title('Streamfunction');
% % colormap(color);
set(gca,'Clim',cc);
hold off;
% colormap('winter')

% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% ylabel('y',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',11,...
%     'FontName','Times')
% xlabel('x',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% axis([x(1) x(end) y(1) y(end)])
% % % colorbar
% % % title('Deformation of a line',...
% % %     'FontUnits','points',...
% % %     'FontWeight','normal',...
% % %     'FontSize',12,...
% % %     'FontName','Times')
% % eval( ['print -depsc ' pwd '/plots/deformed_line.eps']);

%%
% hold on
% n1=8;n2=10;n3=8;
% R1=radius*1.3;R2_k=radius*1.5;R2_r=radius*1.1;
% % R1=radius*1.2;R2_k=radius*1.4;R2_r=radius*1;
% r_k = [ R1*ones(1,n1) R2_k*ones(1,n2) R1*ones(1,n3)];
% r_r = [ R1*ones(1,n1) R2_r*ones(1,n2) R1*ones(1,n3)];
% theta_k = ampli * exp(-1/2*(r_k-radius).^2/sigma^2);
% theta_r = ampli * exp(-1/2*(r_r-radius).^2/sigma^2);
% delta_t = 3600*24*0.1;
% theta_k = cumsum(theta_k)*delta_t;
% theta_r = cumsum(theta_r)*delta_t;
% x_k=r_k.*cos(theta_k)+mean(x);
% y_k=r_k.*sin(theta_k)+mean(y);
% x_r=r_r.*cos(theta_r)+mean(x);
% y_r=r_r.*sin(theta_r)+mean(y);
% plot(x_r,y_r,'r.-');
% plot(x_k,y_k,'y.-');
% 
% nsubsample=20;
% % quiver(x(1:nsubsample:end),y(1:nsubsample:end), ...
% %     w(1:nsubsample:end,1:nsubsample:end,1)',w(1:nsubsample:end,1:nsubsample:end,2)');
% 
% 
% hold off;
% 
% axis equal;
% axis([mean(x) x(end)-mean(x)/3 mean(y) y(end)-mean(y)/3])
% 
% axis xy;
% % title('Velocity norm');
% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% ylabel('y',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',11,...
%     'FontName','Times')
% xlabel('x',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% % colorbar
% % title('$\dot \theta$',...
% %     'FontUnits','points',...
% %     'FontWeight','normal',...
% %     'FontSize',12,...
% %     'interpreter','latex',...
% %     'FontName','Times')
% % axis([x(1) x(end) y(1) y(end)])
% eval( ['print -depsc ' pwd '/plots/diffused_vorticity.eps']);
% 
% keyboard;

%%

% % figure(2)
% figure('Name','Velocity','NumberTitle','off','Units','inches', ...
%     'Position',[X0(1) X0(2) width height], ...
%     'PaperPositionMode','auto');
% x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% % x=x-mean(x);
% y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% % y=y-mean(y);
% imagesc(x,y,sqrt(w(:,:,1).^2+w(:,:,2).^2)');
% axis xy; axis equal;
% % title('Velocity norm');
% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% ylabel('y',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',11,...
%     'FontName','Times')
% xlabel('x',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% colorbar
% title('Norm of the velocity',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'FontName','Times')
% eval( ['print -depsc ' pwd '/plots/norm_w.eps']);
% keyboard;

%%

%%

% % figure(2)
% figure('Name','Angular velocity','NumberTitle','off','Units','inches', ...
%     'Position',[X0(1) X0(2) width height], ...
%     'PaperPositionMode','auto');
% x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% % x=x-mean(x);
% y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% % y=y-mean(y);
% imagesc(x,y,theta_dot');
% axis xy; axis equal;
% % title('Velocity norm');
% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% ylabel('y',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',11,...
%     'FontName','Times')
% xlabel('x',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',11,...
%     'FontName','Times')
% % colorbar
% title('$\dot \theta$',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',12,...
%     'interpreter','latex',...
%     'FontName','Times')
% axis([x(1) x(end) y(1) y(end)])
% eval( ['print -depsc ' pwd '/plots/angular_w.eps']);

% keyboard;

%%

% % figure;imagesc(x,y,theta_dot); axis equal;
% % title('Angular velocity');
% % 
% % omega= - 2 * r.^2/sigma^2 .* theta_dot;
% % figure;imagesc(x,y,omega); axis equal;
% % title('Vorticity');
% % 
% % % startx = mean(x)+std(x)*(0:1:1);
% % % starty = mean(y)*ones(size(startx));
% % % figure;streamline(x,y,w(:,:,1)',w(:,:,2)',startx,starty); axis equal;
% % % title('Streamfunction');
% 
% 
% % keyboard;
% 
% % figure; quiver(x',y',w(:,:,1)',w(:,:,2)'); axis equal; keyboard;
