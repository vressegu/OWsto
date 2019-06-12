function temp = sst_data(X,unity_approx,BOX)

persistent xref yref sstref Mx My dx dy L BOXref unity_approxref

if isempty (xref) || isempty (yref) || isempty (sstref) || ...
        isempty (Mx) || isempty (My) || isempty (dx) || ...
        isempty (dy) || isempty(L)|| isempty(BOXref)|| ...
        isempty(unity_approxref) || any(BOX(:)~=BOXref(:)) ...
        || unity_approxref~=unity_approx
    time_step_choose=5;
    load('data/simu2008_30.mat','sst','dx','dy');
    sstref=sst(3:end,:,time_step_choose);
%     sstref=sst(:,:,time_step_choose);
    clear sst;
    sstref=sstref-mean(mean(sstref));
    
    sstref=sstref';
    
%     image(10*sstref(:,end:-1:1)'+30)
%     
%     sst=sst-mean(mean(sst(:,:,1)));
%     for t=1:30
%     image(50*sst(:,:,t)+30)
% %     image(50*sst(:,:,t)+30)
% drawnow
% pause
%     end
    
    [Mx,My]=size(sstref);
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
%         sstref=fct_unity_approx(sstref,xref,yref,BOX);
%     end
end


%% Boundaries conditions
xi=xref;
yi=yref;
ssti=sstref;
Li=L;
if any(X' < L(:,1) | X' > L(:,2) )
    
%     if any(X' < L(:,1))
        if X(:,1) < L(1,1)
            xi=dx*(-Mx:(Mx-1));
%             xi=dx*(-(Mx-1):Mx);
            ssti=[ssti;ssti];
            Li=[(min(xref)-dx*(Mx)) max(xref) ; min(yref) max(yref) ];
%             Li=[(min(xref)-dx*(Mx+1)) max(xref) ; min(yref) max(yref) ];
        end
        if X(:,1) > L(1,2)
            xi=dx*(0:2*Mx-1);
%             xi=dx*(1:2*Mx);
            ssti=[ssti;ssti];
            Li=[min(xref) (max(xref)+dx*(Mx)) ; min(yref) max(yref) ];
%             Li=[min(xref) (max(xref)+dx*(Mx+1)) ; min(yref) max(yref) ];
        end
        if X(:,2) < L(2,1)
            temp=interp1(xi,ssti(:,1)', X(:,1));
            return
        end
%     end
%     if any(X' > L(:,2))
        if X(:,2) > L(2,2)
            temp=interp1(xi,ssti(:,end)', X(:,1));
            return
        end
        
%     end
end

if any(X' < Li(:,1) | X' > Li(:,2) )
    error('pas normal');
end
%%

% if unity_approx
%     temp =interp2(xref,yref,fct_unity_approx(sstref,xref,yref,BOX)', X(:,1), X(:,2));
% else
if size(X,1)>1
    error('changer les qlqs lignes suivantes');
end
% idx= bsxfun(@and,(abs(xref-X(:,1))/dx<10e-2)' , abs(yref-X(:,2))/dy<10e-2);
idxx= (abs(xi-X(:,1))/dx<10e-2)';
idxy=  abs(yi-X(:,2))/dy<10e-2;
% idx= (xi==X(:,1))' & (yi==X(:,2));
if any(idxx) || any(idxy)
    if any(idxx) && any(idxy)
    temp=ssti(idxx,idxy);
    elseif any(idxx)
        ssty=ssti(idxx,:);
        temp =interp1(yi,ssty, X(:,2));
    else
        sstx=ssti(:,idxy);
        temp =interp1(xi,sstx, X(:,1));
    end
% if any(any(idx))
%     [I,J]=find(idx);
%     temp=ssti(I,J);
else
    temp =interp2(xi,yi,ssti', X(:,1), X(:,2));
end
if isnan(temp)
    error('there is NaN');
end
% end
