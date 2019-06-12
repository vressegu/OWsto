function nabla2_f = laplacian_acc_mat(f,dX,type)
% Compute the laplacian of a matrix field f
% We must have size(f) = [m d Mx My (Mz)]
% where Mx My (Mz) are the space dimensions of R^d.
% The result will be of size size(f)
%

%%

siz = size(f);
d=length(siz)-2;
% d=siz(2);
% if ~(length(siz)==2+d)
%     error('wrong size');
% end

nabla2_f=zeros(siz);
for k_dim=1:d
    d2f=diff_mat(f,k_dim,dX,type);
    d2f=diff_mat(d2f,k_dim,dX,type);
    nabla2_f = nabla2_f + d2f;
end
 

end