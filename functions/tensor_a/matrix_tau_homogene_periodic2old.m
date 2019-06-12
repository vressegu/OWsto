function tau = matrix_tau_homogene_periodic2(model)
% Compute the matrix which define the linear operator tau on the grid
% with periodic boudaries conditions
%

MX=model.grid.MX;
dX=model.grid.dX;
origin=model.kriging.BOX(1,:);
% origin=model.grid.origin;
dx=dX(1);
dy=dX(2);
M=prod(MX);
Mx=MX(1);
My=MX(2);

%% Compute a on the grid
x=(-1:(Mx))*dx+origin(1);
y=(-1:(My))*dy+origin(2);
[x,y]=ndgrid(x,y);
XP=stk_dataframe ([x(:) y(:)]);
clear x y
a= fct_a_prior(XP,model);% M d d
a = reshape(a,[MX+2 2 2]);% MX+2 d d
clear XP

%% Compute the big matrix
% Initialization
tau = zeros([M MX]);

for I=1:M
    [i,j] = f_idxij(I);
    ia=i+1;
    ja=j+1;  
    tau(I,mod(i-3,Mx)+1,j)= a(ia-1,ja,1,1)/(dx^2);
    tau(I,mod(i-2,Mx)+1,mod(j-2,My)+1)= (a(ia-1,ja,1,2)+a(ia,ja-1,1,2))/(dx*dy);
    tau(I,mod(i-2,Mx)+1,mod(j,My)+1)= - (a(ia-1,ja,1,2)+a(ia,ja+1,1,2))/(dx*dy);
    tau(I,i,mod(j-3,My)+1)=a(ia,ja-1,2,2)/(dy^2);
    tau(I,i,j)= - (a(ia+1,ja,1,1)+a(ia-1,ja,1,1))/(dx^2) ...
        - (a(ia,ja+1,2,2)+a(ia,ja-1,2,2))/(dy^2);
    tau(I,i,mod(mod(j,My)+1,My)+1)= a(ia,ja+1,2,2)/(dy^2);
    tau(I,mod(i,Mx)+1,mod(j-2,My)+1)= - (a(ia+1,ja,1,2)+a(ia,ja-1,1,2))/(dx*dy);
    tau(I,mod(i,Mx)+1,mod(j,My)+1)= (a(ia+1,ja,1,2)+a(ia,ja+1,1,2))/(dx*dy);
    tau(I,mod(i+1,Mx)+1,j)= a(ia+1,ja,1,1)/(dx^2);
end

tau=1/8*reshape(tau,[M M]);
% 1/2 because of a/2 and 1/2*1/2 because of the 2nd order derivation scheme

%% Sub functions
%     function idx_I = f_idxI(idx_i,idx_j)
%         idx_I = (idx_i-1)*My + idx_j;
%     end
    function [idx_i,idx_j] = f_idxij(idx_I)
        idx_j = floor((idx_I-1)/Mx)+1;
        idx_i = idx_I -(idx_j-1)*Mx;
%         idx_i = floor((idx_I-1)/My)+1;
%         idx_j = idx_I -(idx_i-1)*My;
        
    end
end

