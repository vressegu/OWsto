function [fft_T_adv] = fct_fft_advection(model, fft_T0, w)
% function [fft_T_adv,fft_Var_T_adv] = fct_fft_advection(model, fft_T0, w, fft_Var_T0)
% Advection of T with the speed w, in Fourrier space
%

% figure(7)
% image(3e-3*abs( fft_T0));
% s_fft_T0=fft_T0;
%
% n=model.grid.MX(1);
% m=model.grid.MX(2);
% fft_T0=fft_T0(1:(n/2+1),:);
% fft_w=fft_w(1:(n/2+1),:,:);
% % fft_T0=fft_T0(1:(n/2+1),1:(m/2+1));
% % fft_w=fft_w(1:(n/2+1),1:(m/2+1),:);

My = size(fft_T0,2);
if model.mirror
    My=My/2;
end

dt=model.advection.dt_adv;

model.advection.advection_duration=fct_time_adv_diff(model,w);
% fct_time_adv_diff(model,w)/3600/24
if strcmp(model.type_data,'erwan')
    model.advection.advection_duration = 1e-5*24*3600;
end

N_t = ceil(model.advection.advection_duration/dt);
fprintf(['Time of advection : ' num2str(N_t*dt/3600/24) ' days \n']);
% % warning('w is aliased');
% warning('on pourrait integrer seulement sur un quadrant');

%% Hyperviscosity
if model.advection.HV.bool
    model.advection.HV.order=8;
    %     model.advection.HV.order=4;
    model.advection.HV.val= ...
        40 * model.advection.lambda_RMS * ...
        (mean(model.grid.dX)/pi)^model.advection.HV.order;
end

%% Anti aliasing filter

fft_T0= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',fft_T0);
fft_T0= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,fft_T0);
w=fft2(w);
w(:,:,1)= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',w(:,:,1));
w(:,:,1)= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,w(:,:,1));
w(:,:,2)= bsxfun(@times, fct_unity_approx5(model.grid.MX(1))',w(:,:,2));
w(:,:,2)= bsxfun(@times, fct_unity_approx5(model.grid.MX(2)) ,w(:,:,2));
w=real(ifft2(w));

%%
% model.advection.evol_variance=false;
% if nargin >3
%     fft_Var_T_adv = fft_Var_T0;
%     model_var=model;
%     model_var.evol_variance=true;
% end
fft_T_adv = fft_T0;
t_last_plot=-inf;
x = model.grid.dX(1)*(0:model.grid.MX(1)-1);
if model.mirror
    y = model.grid.dX(2)*(0:model.grid.MX(2)/2-1);
else
    y = model.grid.dX(2)*(0:model.grid.MX(2)-1);
end
% clear fft_T0
for t=1:N_t
    %     pause;
    %     if nargin > 3
    % %         fft_Var_T_adv = RK4_fft_advection(model_var,fft_Var_T_adv, w);
    %         fft_T_adv(:,:,2)=fft_Var_T_adv;
    %     end
    fft_T_adv = RK4_fft_advection(model,fft_T_adv, w);
    
    %         % Plot
    if (t - t_last_plot) >= 10
%     if (t - t_last_plot)*dt >= 3600*24*0.1
%         %     if mod(t,10)==0
%         %     if mod(t*dt,3600*24*0.1)==0
        t_last_plot=t;
        
        if ~ model.advection.evol_variance
            figure(4)
            subplot(2,1,1)
            fct_spectrum(model, fft_T0);
            hold on;
            fct_spectrum(model,fft_T_adv);
            hold off
            ax = axis;
            ax(3)=0;
            axis(ax);
            
            %         contour(abs( fft_T_adv(2:end-1,2:end-1) ));
            
            T_adv = real(ifft2( fft_T_adv ));
            %         figure(5)
            subplot(2,1,2)
            imagesc(x,y,T_adv(:,1:My)')
            %             image(10*T_adv(:,1:My)'+20)
            
            [n_grad_T0, n_T0] = norm_tracer_tot(model, fft_T0);
            [n_grad_T, n_T] = norm_tracer_tot(model, fft_T_adv);
            filter_width = n_grad_T0*n_grad_T /( n_T0 * (n_grad_T - n_grad_T0));
            fft_filter = design_filter2(model,1/sqrt(filter_width));
            fft_T_adv2 = fft_filter.*fft_T_adv;
            
            figure(5)
            subplot(2,1,1)
            fct_spectrum(model, fft_T0);
            hold on;
            fct_spectrum(model,fft_T_adv2);
            hold off
            ax = axis;
            ax(3)=0;
            axis(ax);
            
            T_adv = real(ifft2( fft_T_adv2 ));
            subplot(2,1,2)
            imagesc(x,y,T_adv(:,1:My)')
            axis xy
            
        else
            %             figure(4)
            %             subplot(2,2,1)
            %             fct_spectrum(model, fft_T0(:,:,1));
            %             hold on;
            %             fct_spectrum(model,fft_T_adv(:,:,1));
            %             hold off
            % %             keyboard;
            %             ax = axis;
            %             if ~strcmp(model.type_data,'toy')
            %                 ax(3)=0;
            %             end
            %             axis(ax);
            
            T_adv = real(ifft2( fft_T_adv(:,:,1) ));
            
            %             subplot(2,2,2)
            %             imagesc(x,y,T_adv(:,1:My)')
            %             axis xy
            
            if t_last_plot==1
%                 if strcmp(model.type_data,'toy')
%                     close all
                    width=2.5;
                    height=2;
                    figure('Name','Tracer','NumberTitle','off','Units','inches', ...
                        'Position',[0 0 width height], ...
                        'PaperPositionMode','auto');
%                 else
%                     figure(4)
%                 end
            end
            imagesc(x,y,T_adv(:,1:My)')
            axis xy
            axis equal
            
            if strcmp(model.type_data,'toy')
            set(gca,...
                'Units','normalized',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontSize',11,...
                'FontName','Times')
            ylabel('y',...
                'FontUnits','points',...
                'interpreter','latex',...
                'FontSize',11,...
                'FontName','Times')
            xlabel('x',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontSize',11,...
                'FontName','Times')
            axis([x(1) x(end) y(1) y(end)])
            
% %             colormap(gray);
%             model.grid.MX(2)=model.grid.MX(2)/2;
%             fct_w_toy(model);
%             model.grid.MX(2)=model.grid.MX(2)*2;
            
                %             eval( ['print -depsc ' pwd '/plots/deformed_line_t='...
                %                 num2str(t) '.eps']);
                eval( ['print -depsc ' pwd '/plots/adv_diff_deformed_tracer_t='...
                    num2str(t) '.eps']);
            end
            %%
%             
%             
%             T_Var_adv = real(ifft2( fft_T_adv(:,:,2) ));
%             %             sum((T_Var_adv(:) < 0))/numel(T_Var_adv < 0)
%             T_Var_adv(T_Var_adv< 0)=0;
%             % %             T_Var_adv=sqrt(T_Var_adv(:,1:My));
%             
%             subplot(2,2,3)
%             fct_spectrum(model,fft2(sqrt(T_Var_adv)));
%             
%             %             fct_spectrum(model, fft_T0(:,:,2));
%             %             hold on;
%             %             fct_spectrum(model,fft_T_adv(:,:,2));
%             %             hold off
%             
%             %             ax = axis;
%             %             ax(3)=0;
%             %             axis(ax);
%             subplot(2,2,4)
%             %             image(1e2*T_Var_adv(:,1:My)')
%             imagesc(x,y,sqrt(T_Var_adv(:,1:My)'))
%             axis xy
%             
%             %             image(3e3*T_Var_adv+30)
%             %             keyboard;
%             %             contour(T_Var_adv(:,1:My))
%             %             pcolor(T_Var_adv(:,1:My))
%%
        end
        drawnow
        fprintf([ num2str(t*dt/(24*3600)) ' days of advection \n'])
%         close
    end
end