function fct_plot_tracer_alpha_comp(model,v_alpha2_m,v_alpha2_m_estim,...
            v_norm_grad_T,...
            v_norm_grad_T_estim,dt,vect_t_plot)
% Plot the sqaure of the L^2 norm of the tracer gradient and of alpha as a
% function of time

t = vect_t_plot(end);
day = num2str(floor(t*dt/(3600*24)));

width=9;
height=3;
%         width=12;
%         height=4;
figure10=figure(10);
set(figure10,'Units','inches', ...
    'Position',[8 1 width height], ...
    'PaperPositionMode','auto');

subplot(1,2,1)
loglog(vect_t_plot*dt,v_alpha2_m,'b');
hold on
loglog(vect_t_plot*dt,v_alpha2_m_estim,'r');
hold off
grid on;
title('Mean $\alpha^2$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
ylabel('$\alpha^2$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('$Time (s)$',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')

% iii= (v_norm_grad_T < 1);
% v_norm_grad_T(iii)=1;
% iii= (v_norm_grad_T_estim < 1);
% v_norm_grad_T_estim(iii)=1;

subplot(1,2,2)
loglog(vect_t_plot*dt, ...
    v_norm_grad_T,'b');
hold on;
loglog(vect_t_plot*dt, ...
    v_norm_grad_T_estim,'r');
hold off;
%     v_norm_grad_T ./ norm_grad_T0);
% %             max(zeros(size(v_norm_grad_T)),v_norm_grad_T-norm_grad_T0));
 ax = axis; ax(3) = 0.7; axis(ax);
% ax = axis;  ax(1) =1e4; 
% if ax(2)>ax(1)
%     axis(ax);
% end
grid on;
title('$\| \nabla T \|^2_{L^2}$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times');
set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
% % %         ylabel('$\| \nabla T \|^2_{L^2} - \| \nabla T_0 \|^2_{L^2} $',...
% %         ylabel(['$\| \nabla T \|^2_{(\mathcal{L}^2 (\mathbb{R}^2))^2} ' ...
% %                 '/ \| \nabla T_0 \|^2_{(\mathcal{L}^2 (\mathbb{R}^2))^2} $'],...
% ylabel('$\| \nabla T \|^2_{L^2} / \| \nabla T_0 \|^2_{L^2} $',...
ylabel('$\| \nabla T \|^2_{L^2} / \| \nabla T_0 \|^2_{L^2} $',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',11,...
    'FontName','Times')
xlabel('$Time(s)$',...
    'interpreter','latex',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',11,...
    'FontName','Times')
eval( ['print -depsc ' model.folder.folder_simu ...
    '/evol_time_alpha/' num2str(day) '.eps']);
%         keyboard;