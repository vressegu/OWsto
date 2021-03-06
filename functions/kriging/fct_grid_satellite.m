function Xn = fct_grid_satellite(model,L,theta,dx_under_trace,dx_maille)
% function Xn = fct_grid_satellite(L,theta,dx_under_trace,dx_maille)
% L = [ Lx Ly ] is the size of the domain
% We assume that the origin (dx,0) is at the middle (in high) left


% if nargin < 2
%     load('data/simu2008_30.mat','ssh','dx','dy');
%     [My,Mx,~]=size(ssh);clear ssh
%     My=My-2;
%     L = [ dx*[0 Mx-1]  ; dy*[0 My-1] ];
%     L_period = [dx*Mx dy*My]; clear Mx My
% %      clear Mx My dx dy
% %     L = [ dx*[1 Mx]  ; dy*[1 My] ]; clear Mx My dx dy
% %     % L = 1e6* [ 0.0021 1.0890 ; 0.0019 1.9440];
% %     % %    L = 1e6* [ 0.0021 2.1823 ; -0.4841 0.4841];
% end
% if nargin < 1
%     BOX = L;
% else
%     BOX=BOX';
% end
if nargin < 2
    BOX = model.grid.BOX;
    BOX(1,:) = BOX(1,:) - model.grid.origin(1);
    BOX(2,:) = BOX(2,:) - model.grid.origin(2);
end
if nargin < 3
    %    theta =0;
    % %    theta = atan((BOX(2,2)-BOX(2,1))/(BOX(1,2)-BOX(1,1)))- pi/4;
    %        theta = atan(L(2,2)/L(1,2))- pi/4;
    L_period = model.grid.MX .* model.grid.dX;
    theta = atan(L_period(2)/L_period(1))- pi/4;
end
if nargin < 4
    %     dx_under_trace = 10e3; % 10 km
    dx_maille =  100e3; % 100 km
    % %     n_maille=ceil((BOX(1,2)-BOX(1,1))/dx_maille);
    % %     dx_maille = (BOX(1,2)-BOX(1,1))/n_maille;
    %     n_maille=ceil(L(1,2)/dx_maille);
    %     dx_maille = L(1,2)/n_maille;
    n_maille=ceil(L_period(1)/dx_maille);
    dx_maille = L_period(1)/n_maille;
    
    dx_under_trace = 10e3; % 10 km
    L0=sqrt(L_period(1)^2+L_period(2)^2);
    %     L0=sqrt(L(1,2)^2+L(2,2)^2);
    n_maille_ut=ceil(L0/dx_under_trace);
    dx_under_trace = L0/n_maille_ut;
end

% nb_trace = 3*ceil((L(1,2)-L(1,1))/dx_maille);
% origin_trace = dx_maille *(-nb_trace/3:2*nb_trace/3-1);
% % nb_trace = ceil((L(1,2)-L(1,1))/dx_maille);
% % origin_trace = dx_maille *(0:nb_trace-1);

origin_trace = dx_maille *(-2*n_maille:3*n_maille-1);
% origin_trace = dx_maille *(-n_maille:2*n_maille-1);
nb_trace = length(origin_trace);

Xn=[];
for k = 1:nb_trace;
    Xn=[Xn double_trace(origin_trace(k)) ];
end

% Add the origin
xn=Xn(1,:);
yn=Xn(2,:);
% xn=Xn(1,:)+BOX(1,1);
% yn=Xn(2,:)+BOX(2,1);
% % yn=Xn(2,:)+(BOX(2,1)+BOX(2,2))/2;
% % % xn=Xn(1,:);
% % % yn=Xn(2,:);

%%

% Remove external points

% Length_BOX=(BOX(:,2)-BOX(:,1))/2;
% BOX(:,1)=BOX(:,1)-Length_BOX;
% BOX(:,2)=BOX(:,2)+Length_BOX;
%%

% % BOX(:,1)=max(BOX(:,1)-dx_maille/2,L(:,1));
% % BOX(:,2)=min(BOX(:,2)+dx_maille/2,L(:,2));
% BOX(:,1)=BOX(:,1)-dx_maille/2;
% BOX(:,2)=BOX(:,2)+dx_maille/2;

if model.grid_trace_larger
    BOX(1,1)=BOX(1,1)-3*dx_maille;
    BOX(1,2)=BOX(1,2)+3*dx_maille;
    BOX(2,1)=BOX(2,1)-3*dx_maille;
    BOX(2,2)=BOX(2,2)+3*dx_maille;
%     BOX(1,1)=BOX(1,1)-3*dx_maille;
%     BOX(1,2)=BOX(1,2)+3*dx_maille;
%     BOX(2,1)=BOX(2,1)-dx_maille/2;
%     BOX(2,2)=BOX(2,2)+dx_maille/2;
else
    marge = 0;
    BOX(1,1)=BOX(1,1)+marge*dx_maille;
    BOX(1,2)=BOX(1,2)-marge*dx_maille;
    BOX(2,1)=BOX(2,1)+marge*dx_maille;
    BOX(2,2)=BOX(2,2)-marge*dx_maille;
end

idx = (BOX(1,1)<=xn & xn<=BOX(1,2));
idx = (BOX(2,1)<=yn & yn<=BOX(2,2)) & idx;
% % idx = (BOX(1,1)<=xn & xn<BOX(1,2));
% % idx = (BOX(2,1)<=yn & yn<BOX(2,2)) & idx;

xn=xn(idx);
yn=yn(idx);

xn = xn +model.grid.origin(1);
yn = yn +model.grid.origin(2);

Xn=[xn;yn];
% mXn=max(Xn,[],2);
% keyboard;
% idxXn=any((Xn==mXn));
% Xn(:,idxXn)=[];

    function Xi = double_trace(origin_tr)
        Lmax=2*max(L_period);
        %         Lmax=2*max(L(:,2));
        % %         Lmax=sqrt(2)*max(L(:,2))/2;
        % %         Lmax=sqrt(2)*max(L(:,2)-L(:,1))/2;
        
        % y increases with x
        rho = dx_under_trace * ...
            [0:-1:-Lmax/dx_under_trace 1:Lmax/dx_under_trace];
        %         rho = -Lmax:dx_under_trace:Lmax;
        alpha = pi/4 + theta;
        xi = rho*cos(alpha);
        yi = rho*sin(alpha);
        %         origin_trace_y=1e6*1.9440/2;
        origin_trace_y=0;
        %         Xi = [origin_tr+xi origin_tr+xi+dx_under_trace/2 ; origin_trace_y+yi origin_trace_y-yi];
        
        Xi = [origin_tr+xi ; origin_trace_y+yi ];
        
%         figure;plot(Xi(1,:),Xi(2,:),'.r');
        
        % y decreases with x
        rho = dx_under_trace * ...
            [-1/2:-1:-Lmax/dx_under_trace 1/2:Lmax/dx_under_trace];
        xi = rho*cos(alpha);
        yi = rho*sin(alpha);
        Xi = [ Xi ...
            [origin_tr+xi ; origin_trace_y-yi] ];
        
%         hold on;
%         plot(Xi(1,:),Xi(2,:),'.b');
%         hold off;
        
    end

end