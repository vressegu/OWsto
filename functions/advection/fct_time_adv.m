function time_adv = fct_time_adv(model,w)
% Estimate the advection time when a filament,
% initialy separated from l_cri = sqrt(prod(model.grid.dX)), is folded
% around a vortice (with one tour)

v_ortho = w(:,:,[2 1]);
v_ortho(:,:,1) = - v_ortho(:,:,1);
v_ortho=v_ortho(3:end-2,3:end-2,:);
%         v_ortho(2:end-1,2:end-1,:)=v_ortho(2:end-1,2:end-1,:);

vort = 1/(2*model.grid.dX(1)) ...
    *(w(3:end,2:end-1,2)-w(1:end-2,2:end-1,2)) ...
    - 1/(2*model.grid.dX(2)) ...
    *(w(2:end-1,3:end,1)-w(2:end-1,1:end-2,1));% Mx-2 My-2
% vort = 1/(2*model.grid.dX(1)) ...
%     *(w(3:end,2:end-1,2)-w(2:end-1,1:end-2,2)) ...
%     - 1/(2*model.grid.dX(2)) ...
%     *(w(2:end-1,3:end,1)-w(2:end-1,1:end-2,1));% Mx-2 My-2
nabla_vort(:,:,1) = 1/(2*model.grid.dX(1)) ...
    *(vort(3:end,2:end-1)-vort(1:end-2,2:end-1)); % Mx-4 My-4 2
%     *(vort(3:end,2:end-1)-vort(2:end-1,1:end-2)); % Mx-4 My-4 2
nabla_vort(:,:,2) = 1/(2*model.grid.dX(2)) ...
    *(vort(2:end-1,3:end)-vort(2:end-1,1:end-2));% Mx-4 My-4 2
time_adv = sum(v_ortho.*nabla_vort,3);
% time_adv = time_adv * 2;

% model.l_pixel_nugget=10e3;
% Ld=exp(-model.kriging.param(2));
% l_cri = model.l_pixel_nugget;
l_cri = sqrt(prod(model.grid.dX));
time_adv = 2*pi* sqrt(sum(v_ortho.^2,3))/l_cri ./ abs(time_adv);
%         time_adv = time_adv/3600/24;


%%
% My = size(w,2)/2;
% x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
% % y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
% y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
% 
% % figure(5)
% width=3;
% height=4;
% % width=2.5;
% % height=2;
% figure('Name','Tracer','NumberTitle','off','Units','inches', ...
%     'Position',[0 0 width height], ...
%     'PaperPositionMode','auto');
% imagesc(x,y,1./(time_adv(:,1:My)'))
% % imagesc(x,y,1./(time_adv'))
% % imagesc(x,y,-log(time_adv'))
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
% %         title('$log(1/T)$',...
% title('1/T',...
%     'FontUnits','points',...
%     'FontWeight','normal',...
%     'interpreter','latex',...
%     'FontSize',12,...
%     'FontName','Times');
%  axis equal;
%  axis xy ;
%  axis([x(1) x(end) y(1) y(end)])
% % fct_w_toy(model);
% % eval( ['print -depsc ' pwd '/plots/t_adv_without_diff_toy.eps']);
% % imagesc(x,y,-log(time_adv(:,1:My)'))

        
        %%
%         
% % [xgrid,ygrid]=ndgrid(x,y);
% % XP=stk_dataframe ([xgrid(:) ygrid(:)]);
% % zp = stk_predict (model.kriging, model.obs.x,model.obs.h, XP);
% % save(['data_temp_' num2str(model.grid.MX(1)) '_x_' num2str(model.grid.MX(2)/2)],'model','zp','XP');
% 
% hold on;
% load(['data_temp_' num2str(model.grid.MX(1)) '_x_' num2str(model.grid.MX(2)/2)])
% ZZ=reshape(double(zp.mean),[model.grid.MX(1) model.grid.MX(2)/2] );
% hold on;
% clim=get(gca,'Clim');
% contour(x,y,ZZ','k')
% set(gca,'Clim',clim);
% hold off
% 
% %         set(gca,...
% %             'Units','normalized',...
% %             'FontUnits','points',...
% %             'FontWeight','normal',...
% %             'FontSize',11,...
% %             'FontName','Times')
% %         ylabel('y',...
% %             'FontUnits','points',...
% %             'interpreter','latex',...
% %             'FontSize',11,...
% %             'FontName','Times')
% %         xlabel('x',...
% %             'FontUnits','points',...
% %             'FontWeight','normal',...
% %             'FontSize',11,...
% %             'FontName','Times')
% % %         title('$log(1/T)$',...
% %         title('Mezic',...
% %             'FontUnits','points',...
% %             'FontWeight','normal',...
% %             'interpreter','latex',...
% %             'FontSize',12,...
% %             'FontName','Times');
% %         axis equal;axis xy;
% %         axis([x(1) x(end) y(1) y(end)])
%         eval( ['print -depsc ' pwd '/plots/time_adv_klein.eps']);
% %         eval( ['print -depsc ' pwd '/plots/mezic_crtiterion.eps']);
% % 
% keyboard;
%%


time_adv = min(time_adv(:));