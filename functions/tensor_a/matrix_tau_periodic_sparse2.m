function tau = matrix_tau_periodic_sparse2(model)
% Compute the matrix which define the linear operator tau on the grid
% with periodic boudaries conditions on null conditions on y=0
%
warning(['The operator tau is defined for periodic functions '...
    'with homogene Diriclet conditions on y=0. ' ...
    'Thus, its matrix representation should acts on the part of the function ' ...
    'which such as y > 0. ' ...
    'Moreover, to optain tau(f) on as to add 0 for y=0 values.']);
boundary_cond='homo_dirichlet_conditions_on_y';
% 'homo_dirichlet_conditions' 'homo_dirichlet_conditions_on_y'

MX=model.grid.MX;
dX=model.grid.dX;
origin=model.kriging.BOX(1,:);
% origin=model.grid.origin;
dx=dX(1);
dy=dX(2);
Mx=MX(1);
My=MX(2);
MX=MX-[0 1];
M=prod(MX);

%% Compute a on the grid
x=(-2:(Mx+1))*dx+origin(1);
% y=(-2:(My+1))*dy+origin(2);
y=(-1:(My+1))*dy+origin(2);
My=My-1;
% x=(-1:(Mx))*dx+origin(1);
% y=(-1:(My))*dy+origin(2);
[x,y]=ndgrid(x,y);
XP=stk_dataframe ([x(:) y(:)]);
clear x y
a= fct_a(XP,model);% M d d
a = reshape(a,[MX+4 2 2]);% MX+2 d d
% a = reshape(a,[MX+2 2 2]);% MX+2 d d
a=bsxfun(@times,permute(1./(dX'*dX),[3 4 1 2]), a);
clear XP

%% Compute the big matrix
% Initialization
% tau = zeros([M MX]);
tau = sparse([],[],[],M,M);

for I=1:M
    [ii,jj] = f_idxij(I);
    %     tau(I,:,:)= tau_I_1(ii,jj,1) ;
    %     tau(I,:,:)= 1/4* tau_I_1(ii,jj,2);
    %     tau(I,:,:)= (16/9) * tau_I_1(ii,jj,1) - 2/9 * tau_I_2(ii,jj) ...
    %                 + 1/36 * tau_I_1(ii,jj,2);
    tau(I,:)= (16/9) * tau_I_1(ii,jj,1) - 2/9 * tau_I_2(ii,jj) ...
        + 1/36 * tau_I_1(ii,jj,2);
end

tau=1/8*tau;
% tau=1/8*reshape(tau,[M M]);
% 1/2 because of a/2 and 1/2*1/2 because of the 2nd order derivation scheme

%% Boundary conditions
% ze=false(MX);
% switch boundary_cond
%     case homo_dirichlet_conditions
%         % Force =0 on left boundaries
%         ze(1,:)=true;
%         ze(:,1)=true;
%     case homo_dirichlet_conditions_on_y
%         % Force =0 on left boundaries on y
%         ze(:,1)=true;
%     otherwise
%         error('unknown bounday contion');
% end
% ze=ze(:);
% tau(ze,:)=0;

%% Sub functions
    function [idx_i,idx_j] = f_idxij(idx_I)
        idx_j = floor((idx_I-1)/Mx)+1;
        idx_i = idx_I -(idx_j-1)*Mx;
    end

    function tau_I = tau_I_1(i,j,step)
%         tau_I = zeros(MX);
        tau_I = sparse([],[],[],MX(1),MX(2));
        ia=i+2;
        ja=j+2;
        
        tau_I(mod(i-2*step-1,Mx)+1,j)= a(ia-step,ja,1,1);
        tau_I(mod(i-step-1,Mx)+1,mod(j-step-1,My)+1)= (a(ia-step,ja,1,2)+a(ia,ja-step,1,2));
        tau_I(mod(i-step-1,Mx)+1,mod(j+step-1,My)+1)= - (a(ia-step,ja,1,2)+a(ia,ja+step,1,2));
        tau_I(i,mod(j-2*step-1,My)+1)=a(ia,ja-step,2,2);
        tau_I(i,j)= - (a(ia+step,ja,1,1)+a(ia-step,ja,1,1)) ...
            - (a(ia,ja+step,2,2)+a(ia,ja-step,2,2));
        tau_I(i,mod(j+2*step-1,My)+1)= a(ia,ja+step,2,2);
        tau_I(mod(i+step-1,Mx)+1,mod(j-step-1,My)+1)= - (a(ia+step,ja,1,2)+a(ia,ja-step,1,2));
        tau_I(mod(i+step-1,Mx)+1,mod(j+step-1,My)+1)= (a(ia+step,ja,1,2)+a(ia,ja+step,1,2));
        tau_I(mod(i+2*step-1,Mx)+1,j)= a(ia+1,ja,1,1);
        
%         tau_I = permute(tau_I,[3 1 2]); % 1 MX
        tau_I =tau_I(:)';
    end

    function tau_I = tau_I_2(i,j)
%         tau_I = zeros(MX);
        tau_I = sparse([],[],[],MX(1),MX(2));
        ia=i+2;
        ja=j+2;
        
        tau_I(mod(i-3-1,Mx)+1,j)= (a(ia-2,ja,1,1)+a(ia-1,ja,1,1));
        tau_I(mod(i-2-1,Mx)+1,mod(j-1-1,My)+1)= (a(ia-2,ja,1,2)+a(ia,ja-1,1,2));
        tau_I(mod(i-2-1,Mx)+1,mod(j+1-1,My)+1)= -(a(ia-2,ja,1,2)+a(ia,ja+1,1,2));
        tau_I(mod(i-1-1,Mx)+1,mod(j-2-1,My)+1)= (a(ia,ja-2,1,2)+a(ia-1,ja,1,2));
        tau_I(mod(i-1-1,Mx)+1,j)= -(a(ia-2,ja,1,1)+a(ia+1,ja,1,1));
        tau_I(mod(i-1-1,Mx)+1,mod(j+2-1,My)+1)= - (a(ia,ja+2,1,2)+a(ia-1,ja,1,2));
        
        tau_I(i,mod(j-3-1,My)+1)= (a(ia,ja-2,2,2)+a(ia,ja-1,2,2));
        tau_I(i,mod(j-1-1,My)+1)= -(a(ia,ja-2,2,2)+a(ia,ja+1,2,2));
        tau_I(i,mod(j+1-1,My)+1)= -(a(ia,ja+2,2,2)+a(ia,ja-1,2,2));
        tau_I(i,mod(j+3-1,My)+1)= (a(ia,ja+2,2,2)+a(ia,ja+1,2,2));
        
        tau_I(mod(i+1-1,Mx)+1,mod(j-2-1,My)+1)= -(a(ia,ja-2,1,2)+a(ia+1,ja,1,2));
        tau_I(mod(i+1-1,Mx)+1,j)= -(a(ia+2,ja,1,1)+a(ia-1,ja,1,1));
        tau_I(mod(i+1-1,Mx)+1,mod(j+2-1,My)+1)= (a(ia,ja+2,1,2)+a(ia+1,ja,1,2));
        tau_I(mod(i+2-1,Mx)+1,mod(j-1-1,My)+1)= -( a(ia+2,ja,1,2)+a(ia,ja-1,1,2));
        tau_I(mod(i+2-1,Mx)+1,mod(j+1-1,My)+1)= (a(ia+2,ja,1,2)+a(ia,ja+1,1,2));
        tau_I(mod(i+3-1,Mx)+1,j)= (a(ia+2,ja,1,1)+a(ia+1,ja,1,1));
        
%         tau_I = permute(tau_I,[3 1 2]); % 1 MX
        tau_I =tau_I(:)';
    end
end

