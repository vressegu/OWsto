function [on_tau2_global, on_tau_local] = fct_on_tau_2_stretching(model,w_eul)
% Compute the inverse of the local time scale of differntial advection
%

% MX=model.grid.MX;
% x=model.grid.dX(1)* (0:(MX(1)-1)) ;
% x=x-mean(x);
% y=model.grid.dX(2)* (0:(MX(2)-1)) ;
% y=y-mean(y);
% [x,y]=ndgrid(x,y);
% r=sqrt(x.^2+y.^2);
% 
% MX=model.grid.MX;
% w_Lag=permute(w_Lag,[1 2 4 3]);
% nabla_phi=permute(nabla_phi,[1 2 4 3]);

% figure;
% norm_nabla_phi = sum(sum(nabla_phi.^2,4),3);
% imagesc(norm_nabla_phi);
% % imagesc(1./norm_nabla_phi);

% Normalized orthogonal Lagrangian velocity
v_tilde_ortho = fct_ortho(fct_normalize(w_eul));
% v_tilde_ortho = fct_ortho(fct_normalize(w_Lag));

% % figure;
% % norm = sum(v_tilde_ortho.^2,3);
% % % norm = sum(v_tilde_ortho.^2,3);
% % imagesc(norm);
% figure;
% quiver(v_tilde_ortho(:,:,1)',v_tilde_ortho(:,:,2)'); 
% axis xy; axis equal

% Norm of the Eulerian velocity
n_v = sqrt(sum(w_eul.^2,3));
% n_v = n_v ./ r;

% figure;
% imagesc(n_v)

% % Curvature
% w_eul_n = fct_normalize(w_eul);
% w_eul_n_ortho = fct_ortho(w_eul_n);
% 
% % Gradient of the norm of the Eulerian velocity
% grad_w_eul_n = permute( w_eul_n , [1 2 4 3]);
% grad_w_eul_n = gradient_perso(model.grid, grad_w_eul_n );
% grad_w_eul_n = bsxfun(@times, w_eul_n, grad_w_eul_n);
% grad_w_eul_n = sum(grad_w_eul_n,3);
% grad_w_eul_n = permute( grad_w_eul_n , [1 2 4 3]);
% curv = bsxfun(@times, w_eul_n_ortho, grad_w_eul_n);
% curv = sum(curv,3);

% figure;
% % imagesc(curv);
% % imagesc(curv.*r);
% mask = ( n_v > 1e-3 * max(n_v(:)) );
% % imagesc(mask);
% imagesc(curv.*r.*mask);

% % Using local angular velocity instead of velocity
% n_v = n_v .* curv;

% Gradient of the norm of the Eulerian velocity
grad_n_v = gradient_perso(model.grid, n_v );

% figure;
% norm = sum(grad_n_v.^2,3);
% imagesc(norm);

% % Remove boundary pixels
% v_tilde_ortho=v_tilde_ortho(2:MX(1)-1,2:MX(2)-1,:);
% grad_n_v=grad_n_v(2:MX(1)-1,2:MX(2)-1,:);
% MX=MX-2;

% % Inversion of nabla_phi
% on_tau_local = nan([2 MX]);
% % on_tau_local = nan([MX 2]);
% grad_n_v=permute(grad_n_v,[3 1 2]);
% nabla_phi=permute(nabla_phi,[3 4 1 2]);
% for i=1:MX(1)
%     for j=1:MX(2)
%         on_tau_local(:,i,j) = nabla_phi(:,:,i,j) \ grad_n_v(:,i,j);
% %         on_tau_local(i,j,:) = nabla_phi(i,j,:,:) \ grad_n_v(i,j,:);
%     end
% end
% on_tau_local=permute(on_tau_local,[2 3 1]);

% Derivation at the initial time othogonal to the stream direction
on_tau_local = 1/sqrt(2) * sum( v_tilde_ortho .* grad_n_v , 3) ;
% on_tau_local = sum( v_tilde_ortho .* on_tau_local , 3) ;

% figure(25);
% subplot(2,1,1)
% imagesc(v_tilde_ortho(:,:,1)');axis xy;axis equal;colorbar;
% subplot(2,1,2)
% imagesc(v_tilde_ortho(:,:,2)');axis xy;axis equal;colorbar;

% ngr = 1/2*grad_n_v(:,:,1).^2;
% % ngr = 1/2*sum(sum(grad_n_v.^2,4),3);
% figure(25)
% imagesc(ngr');axis xy;axis equal;colorbar;

% grad_v = gradient_perso(model.grid, permute(w_eul,[1 2 4 3]) );
% ngr = 1/2*sum(sum(grad_v.^2,4),3);
% % imagesc(ngr');axis xy;axis equal;colorbar;

% figure(25)
% subplot(2,2,1)
% imagesc(1/2*(grad_v(:,:,1,1).^2)');axis xy;axis equal;colorbar;
% subplot(2,2,2)
% imagesc(1/2*(grad_v(:,:,2,1).^2)');axis xy;axis equal;colorbar;
% subplot(2,2,3)
% imagesc(1/2*(grad_v(:,:,1,2).^2)');axis xy;axis equal;colorbar;
% subplot(2,2,4)
% imagesc(1/2*(grad_v(:,:,2,2).^2)');axis xy;axis equal;colorbar;


on_tau_local = abs(on_tau_local);

% warning('bidouille tau stretch');
% on_tau_local =4*on_tau_local;

on_tau2_global = 1/prod(model.grid.MX) * sum(on_tau_local(:).^2);
tau_global = 1/sqrt(on_tau2_global) /(3600*24)
% keyboard;


%% Plot

% % Grid
% x = model.grid.x_ref;
% y = model.grid.y_ref;
% 
% % Other parameters
% taille_police = 12;
% id_part=1;
% type_data = model.type_data;
% folder_simu = model.folder.folder_simu;
% map = 'default';
% % map = model.folder.colormap;
% 
% width = 3.3;
% %     height = 3.2;
% ax = [x(end)-x(1) y(end)-y(1)] ;
% aspect_ratio = ax(2)/ax(1);
% height = aspect_ratio * width;
% X0=[0 10];
% 
% figure23=figure(23);
% set(figure23,'Units','inches', ...
%     'Position',[X0(1) X0(2) width height], ...
%     'PaperPositionMode','auto');
% % imagesc(x,y,log((on_tau_local.^2)'));
% imagesc(x,y,((on_tau_local.^2)'));
% set(gca,...
%     'Units','normalized',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',taille_police,...
%     'FontName','Times')
% ylabel('y(m)',...
%     'FontUnits','points',...
%     'interpreter','latex',...
%     'FontSize',taille_police,...
%     'FontName','Times')
% xlabel('x(m)',...
%     'interpreter','latex',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'FontSize',taille_police,...
%     'FontName','Times')
% title('\hspace{0.5cm} $1/\tau_s^2$ ',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times')
% axis xy; axis equal
% colormap(map)
% colorbar
% drawnow
% eval( ['print -depsc ' folder_simu '/on_tau_2_stretching.eps']);
% keyboard;


end

function f_ortho = fct_ortho(f)
% Compute the orthogonal vector in each point of the space
%
f_ortho(:,:,1)= - f(:,:,2);
f_ortho(:,:,2)= + f(:,:,1);
end

function g = fct_normalize(g)
% Compute the orthogonal vector in each point of the space
%
ng = sqrt(sum(g.^2,3));
g = bsxfun( @times, 1./ng , g);
end