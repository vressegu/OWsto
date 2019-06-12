function nabla_f = gradient_acc_mat(f,dX,type)
% Compute the gradient of a matrix field f
% We must have size(f) = [m d Mx My (Mz)]
% where Mx My (Mz) are the space dimensions of R^d.
% The result will be of size [ size(f) d]
%

%%

siz = size(f);
d=length(siz)-2;
% d=siz(2);
% if ~(length(siz)==2+d)
%     error('wrong size');
% end
nabla_f=zeros([siz d]);

idx='';
for k_dim=1:d
    idx = [idx ':,'];
end

for k_dim=1:d
    eval(['nabla_f(:,:,' idx 'k_dim)=diff_mat(f,k_dim,dX,type);']);
end
 

end