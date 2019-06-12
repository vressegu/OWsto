function high = ssh_data(X,unity_approx,BOX)

persistent xref yref sshref Mx My dx dy L BOXref unity_approxref

if isempty (xref) || isempty (yref) || isempty (sshref) || ...
        isempty (Mx) || isempty (My) || isempty (dx) || ...
        isempty (dy) || isempty(L)|| isempty(BOXref)|| ...
        isempty(unity_approxref) || any(BOX(:)~=BOXref(:)) ...
        || unity_approxref~=unity_approx
    time_step_choose=5;
    load('data/simu2008_30.mat','ssh','dx','dy');
    sshref=ssh(3:end,:,time_step_choose);
%     sshref=ssh(:,:,time_step_choose);
    clear ssh;
    sshref=sshref-mean(mean(sshref));
    
    sshref=sshref';
    
    %     image(50*sshref(:,end:-1:1)'+30)
    %
    %     ssh=ssh-mean(mean(ssh(:,:,1)));
    %     for t=1:30
    %     image(50*ssh(:,:,t)+30)
    % %     image(50*ssh(:,:,t)+30)
    % drawnow
    % pause
    %     end
    
    [Mx,My]=size(sshref);
    %     dX=[dx dy]
    xref=dx*(0:Mx-1);
    yref=dy*(0:My-1);
%     xref=dx*(1:Mx);
%     yref=dy*(1:My);
    %     yref=dy*((-My/2+1/2):(+My/2-1/2));
    L=[min(xref) max(xref) ; min(yref) max(yref) ];
    
    %     BOX=L';
    BOXref=BOX;
    unity_approxref=unity_approx;
    %     if unity_approx
    %         sshref=fct_unity_approx(sshref,xref,yref,BOX);
    %     end
end

%% Boundaries conditions
xi=xref;
yi=yref;
sshi=sshref;
Li=L;
if any(X' < L(:,1) | X' > L(:,2) )
    
%     if any(X' < L(:,1))
        if X(:,1) < L(1,1)
            xi=dx*(-Mx:(Mx-1));
%             xi=dx*(-(Mx-1):Mx);
            sshi=[sshi;sshi];
            Li=[(min(xref)-dx*(Mx)) max(xref) ; min(yref) max(yref) ];
%             Li=[(min(xref)-dx*(Mx+1)) max(xref) ; min(yref) max(yref) ];
        end
        if X(:,1) > L(1,2)
            xi=dx*(0:2*Mx-1);
%             xi=dx*(1:2*Mx);
            sshi=[sshi;sshi];
            Li=[min(xref) (max(xref)+dx*(Mx)) ; min(yref) max(yref) ];
%             Li=[min(xref) (max(xref)+dx*(Mx+1)) ; min(yref) max(yref) ];
        end
%     end
%     if any(X' > L(:,2))
        if X(:,2) < L(2,1)
            high=interp1(xi,sshi(:,1)', X(:,1));
            return
        end
        if X(:,2) > L(2,2)
            high=interp1(xi,sshi(:,end)', X(:,1));
            return
        end
        
%     end
end

if any(X' < Li(:,1) | X' > Li(:,2) )
    error('pas normal');
end
%%

% % if unity_approx
% %     high =interp2(xref,yref,fct_unity_approx(sshref,xref,yref,BOX)', X(:,1), X(:,2));
% % else
%     high =interp2(xref,yref,sshref', X(:,1), X(:,2));
% % end


high =interp2(xi,yi,sshi', X(:,1), X(:,2));
