function fct_plot_f_under_trace(model,grid_trace_index)
% For each trace, 10 measures on each sides are modified to go smoothly to
% zero on boundaries of the trace. The "usefuel points" are the number of
% measures that are not modified. Then, the spectrum along this trace is
% computed and weighted by the number of "usefuel points".
% The "weighted average spectrum" is the sum of all these spectrum divided
% by the total number of "useful points" on all spectrum's trace.
%

plot_and_weight = false;

%% Get paramters

% Grid
x = model.grid.x;
y = model.grid.y;

% Other parameters
taille_police = 12;
id_part=1;
type_data = model.type_data;
% folder_simu = model.folder.folder_simu;
plot_moments = false;
map = model.folder.colormap;


%% Spectrum under each trace

warning('boundaries are removed to compute the spectrum');

if plot_and_weight
    figure4=figure(1);
end


nb_trace = length(grid_trace_index);
length_trace=nan(1,nb_trace);
length_trace_useful = nan(nb_trace,1);
for k=1:nb_trace
    length_trace(k) = size(grid_trace_index{k},2);
    if plot_and_weight
        [mask_boundaries,N_k_useful] = ...
            fct_unity_approx_(length_trace(k));
        length_trace_useful(k)=N_k_useful;
    end
end
[Mx_max,idx_longest_trace]=max(length_trace);
% even = (floor(Mx_max/2)==Mx_max/2);
P_max=Mx_max/2;
% P_max=Mx_max;
trace_k_max = grid_trace_index{idx_longest_trace};
dX=trace_k_max(3,2)-trace_k_max(3,1);
% Discrete wave number
kappa_max=1/Mx_max * (0:(P_max-1)) ;
kappa_max=2*pi/dX*kappa_max;

spectrum_k_interp = nan(nb_trace,floor(P_max));
% spectrum_k_interp = nan(nb_trace,Mx_max);
for k=1:nb_trace
    trace_k = grid_trace_index{k};
    
    %             figure(1);plot(trace_k(3,:),'o-')
    %             figure(2);plot(trace_k(3,:),trace_k(4,:),'o-')
    
    % Remove boundaries
    [mask_boundaries,N_k_useful] = ...
        fct_unity_approx_(length_trace(k));
    length_trace_useful(k)=N_k_useful;
    
    if plot_and_weight
        figure8=figure(8);plot(mask_boundaries,'o-');
        title('trace without boundaries');
    end
    
    trace_k(4,:) = bsxfun(@times, mask_boundaries, trace_k(4,:));
    
    if plot_and_weight
        X0=[3.3 1];
        widthtemp = 12;
        heighttemp = 6;
        
        figure(1);
        set(figure4,'Units','inches', ...
            'Position',[X0(1) X0(2) widthtemp heighttemp], ...
            'PaperPositionMode','auto');
        
        % Check energy
        %     mean(trace_k(4,:).^2)
        
        if k==idx_longest_trace
            color='g';
        else
            color='b';
        end
    else
        color=[];
    end
    
    [spectrum_k,kappa_k] = ...
        fct_spectrum_1D(trace_k(3,:),fft(trace_k(4,:)),color);
    hold on;
    %     if k==idx_longest_trace
    %         kappa_k
    %     end
    %     hold off;
    %     pause;
    
    %     meas_k_reflected = [trace_k(4,:) trace_k(4,end:-1:1)];
    %     rho_k_reflected = trace_k(3,1) + ...
    %         (trace_k(3,2)-trace_k(3,1))* ...
    %         (0:(length(meas_k_reflected)-1));
    %     [spectrum_k,kappa_k] = ...
    %         fct_spectrum_1D(rho_k_reflected,fft(meas_k_reflected),[]);
    
    %     figure(3);plot(rho_k_reflected,meas_k_reflected,'o-')
    %     meas_k_reflected = meas_k_reflected - mean(meas_k_reflected(:));
    %     mean(meas_k_reflected.^2)
    
    %     X0=[3.3 1];
    %     widthtemp = 12;
    %     heighttemp = 6;
    %     set(figure4,'Units','inches', ...
    %         'Position',[X0(1) X0(2) widthtemp heighttemp], ...
    %         'PaperPositionMode','auto');
    %     [spectrum_k,kappa_k] = ...
    %         fct_spectrum_1D(rho_k_reflected,fft(meas_k_reflected),'b');
    %     trace_k(5,:) = kappa_k;
    %     trace_k(6,:) = spectrum_k;
    
    
    if plot_and_weight
        set(gca,'XGrid','on','XTickMode','manual');
        width = 4.5;
        % width = 4;
        height = 3;
        set(figure4,'Units','inches', ...
            'Position',[X0(1) X0(2) width height], ...
            'PaperPositionMode','auto');
        set(gca,'YGrid','on')
        
        set(gca,...
            'Units','normalized',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'FontName','Times')
        ylabel('$\overline{\Gamma}_T$',...
            'FontUnits','points',...
            'interpreter','latex',...
            'FontSize',taille_police,...
            'FontName','Times')
        xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'FontSize',taille_police,...
            'interpreter','latex',...
            'FontName','Times')
        title('Spectrum of each traces',...
            'FontUnits','points',...
            'FontWeight','normal',...
            'interpreter','latex',...
            'FontSize',12,...
            'FontName','Times')
        drawnow
        
        weight = length_trace_useful(k)/sum(length_trace_useful)
    
    end
    
    spectrum_k_interp(k,:) = interp1(...
        kappa_k,spectrum_k, kappa_max,'spline');
    spectrum_k_interp(k,:) = length_trace_useful(k) ...
        * spectrum_k_interp(k,:);
    %     spectrum_k_interp(k,:) = length_trace(k) ...
    %         * spectrum_k_interp(k,:);
    
    %     pause;
    
    grid_trace_index{k} =  trace_k;
end

spectrum_k_interp = 1/sum(length_trace_useful) * ...
    sum(spectrum_k_interp,1);
% spectrum_k_interp = 1/sum(length_trace) * sum(spectrum_k_interp,1);
% kappa_max

% X0=[3.3 1];
% figure5=figure;
% widthtemp = 12;
% heighttemp = 6;
% set(figure5,'Units','inches', ...
%     'Position',[X0(1) X0(2) widthtemp heighttemp], ...
%     'PaperPositionMode','auto');

if plot_and_weight
    hold on;
    loglog(kappa_max(2:end) , spectrum_k_interp(2:end) ,'r');
    % loglog(kappa_max(2:end) , spectrum_k_interp(2:end) ,'b');
    ax=axis;
    ax(4)=max(spectrum_k_interp(2:end));
    ax(3)=min(spectrum_k_interp(2:end));
    ax(1:2)=kappa_max([2 end]);
    if ax(4)>0
        axis(ax)
    end
end

%% Average spectrum

X0=[3.3 1];
figure5=figure(5);
widthtemp = 12;
heighttemp = 6;
set(figure5,'Units','inches', ...
    'Position',[X0(1) X0(2) widthtemp heighttemp], ...
    'PaperPositionMode','auto');

% hold on;
% loglog(kappa_max(2:end) , spectrum_k_interp(2:end) ,'r');
loglog(kappa_max(2:end) , spectrum_k_interp(2:end) ,'b');
ax=axis;
ax(4)=max(spectrum_k_interp(2:end));
ax(3)=min(spectrum_k_interp(2:end));
ax(1:2)=kappa_max([2 end]);
if ax(4)>0
    axis(ax)
end

% fct_spectrum( model,fft_w(:,:,:,id_part),'b');
set(gca,'XGrid','on','XTickMode','manual');
width = 4.5;
% width = 4;
height = 3;
set(figure5,'Units','inches', ...
    'Position',[X0(1) X0(2) width height], ...
    'PaperPositionMode','auto');
set(gca,'YGrid','on')

set(gca,...
    'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'FontName','Times')
ylabel('$\overline{\Gamma}_T$',...
    'FontUnits','points',...
    'interpreter','latex',...
    'FontSize',taille_police,...
    'FontName','Times')
xlabel('$\kappa \bigl ( rad.m^{-1} \bigr )$',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',taille_police,...
    'interpreter','latex',...
    'FontName','Times')
title('Weighted average spectrum under trace',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'interpreter','latex',...
    'FontSize',12,...
    'FontName','Times')
drawnow
% eval( ['print -depsc ' folder_simu '/Spectrum/' day '.eps']);
% keyboard;

if plot_and_weight
    close(figure4);close(figure8);
end

    function [t,N_t_usefull] = fct_unity_approx_(N_t)
        % XP must be a 2 x n matrix
        % the result is a vector of size n
        %
        
        %         % % slop_size_ratio=10;
        % %         slop_size_ratio=6;
        %         slop_size_ratio=3;
        %         % % slop_size_ratio=ceil(6*(N_t/1000));
        %         % % N_t=1000;
        
        t=ones(1,N_t);
        %P_t=N_t/2;
        
        %         sslop=ceil(N_t/slop_size_ratio);
        sslop=10;
        
        
        % t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
        % t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
        % t(P_t+1)=0;
        t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
        t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
        N_t_usefull = max([N_t-2*sslop 0]);
    end

    function [t,N_t_usefull] = fct_unity_approx_old(N_t)
        % XP must be a 2 x n matrix
        % the result is a vector of size n
        %
        
        % % slop_size_ratio=10;
        %         slop_size_ratio=6;
        slop_size_ratio=3;
        % % slop_size_ratio=ceil(6*(N_t/1000));
        % % N_t=1000;
        
        t=ones(1,N_t);
        %P_t=N_t/2;
        sslop=ceil(N_t/slop_size_ratio);
        % t((P_t-sslop+1):P_t)= (-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
        % t((P_t+2):(P_t+1+sslop))= (tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
        % t(P_t+1)=0;
        t(1:sslop)=(tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) )+1)/2;
        t(end-sslop+1:end)=(-tanh(-3 + 6/(sslop-1)*(0:(sslop-1)) ) +1)/2;
        N_t_usefull = max([N_t-2*sslop 0]);
    end



end

