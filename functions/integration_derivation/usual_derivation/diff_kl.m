function d2a = diff_kl(a,dX)
% a must have the size :  m x d x d x d x Mx x My (x Mz)
% Compute d2a(:,:,k,l,x_1,x_2(,x_3))= d_{x_k} d_{x_l} a(:,:,k,l,x_1,x_2(,x_3))
% ( derivation along the axis x_k and x_l
% Then, the results is sumed along the 3rd and 4th dimension
% The result have the size m x d x Mx x My (x Mz)

siz = size(a); % = [m x d x d x d x Mx x My (x Mz)]
d=length(siz)-4;
% d=siz(3);
% if ~(length(siz)==4+d && siz(3) == siz(4))
%     error('wrong size');
% end

idx='';
for k_dim=1:d
    idx = [idx ',:'];
end

d2a=zeros(siz);

for l=1:d
    eval(['d2a(:,:,:,l' idx ')= diff_l(a(:,:,:,l' idx '),l,dX);']);
end
d2a = permute(d2a,[1 2 4 3 5:ndims(d2a)]);
for k=1:d
    eval(['d2a(:,:,:,k' idx ')= diff_l(d2a(:,:,:,k' idx '),k,dX);']);
end
d2a = permute(d2a,[1 2 4 3 5:ndims(d2a)]);

% d2a = squeeze(sum(sum(d2a,4),3));
d2a = sum(sum(d2a,4),3);
d2a = permute(d2a, [1 2 5:ndims(d2a) 3 4]);

% function da = diff_l(a,l,dX)
% % a must have the size :  m x d x d x 1 x Mx x My (x Mz)
% % Compute da(:,:,:,1,x_1,x_2(,x_3))= d_{x_l} a(:,:,:,1,x_1,x_2(,x_3))
% % 
% 
% siz = size(a);
% d=siz(2);
% if ~(length(siz)==4+d && siz(2)==siz(3) && siz(4) == 1 )
%     error('wrong size');
% end
% dx=dX(l);
% 
% da=zeros(siz);
% if d == 2
%     a = permute(a,[1 2 3 4 5+(l-1) 5+mod(l,2)]); % Put the x_l cordinate in first place
%     da = permute(da,[1 2 3 4 5+(l-1) 5+mod(l,2)]); % Put the x_l cordinate in first place
%     da(:,:,:,:,2:end-1,:) = 1/(2*dx) * (a(:,:,:,:,3:end,:)-a(:,:,:,:,1:end-2,:));
%     da(:,:,:,:,1,:) = 1/(dx) * (a(:,:,:,:,2,:)-a(:,:,:,:,1,:));
%     da(:,:,:,:,end,:) = 1/(dx) * (a(:,:,:,:,end,:)-a(:,:,:,:,end-1,:));
%     da = permute(da,[1 2 3 4 5+mod(1-l,2) 5+mod(2-l,2)]); % Put the first cordinate in first place again
% else
%     a = permute(a,[1 2 3 4 5+(l-1) 5+mod(l,3) 5+mod(l+1,3)]); % Put the x_l cordinate in first place
%     da = permute(da,[1 2 3 4 5+(l-1) 5+mod(l,3) 5+mod(l+1,3)]); % Put the x_l cordinate in first place
%     da(:,:,:,:,2:end-1,:,:) = 1/(2*dx) * (a(:,:,:,:,3:end,:,:)-a(:,:,:,:,1:end-2,:,:));
%     da(:,:,:,:,1,:,:) = 1/(dx) * (a(:,:,:,:,2,:,:)-a(:,:,:,:,1,:,:));
%     da(:,:,:,:,end,:,:) = 1/(dx) * (a(:,:,:,:,end,:,:)-a(:,:,:,:,end-1,:,:));
%     da = permute(da,[1 2 3 4 5+mod(1-l,3) 5+mod(2-l,3) 5+mod(3-l,3)]); % Put the first cordinate in first place again
% 
% end