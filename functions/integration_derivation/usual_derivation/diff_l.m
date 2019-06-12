function da = diff_l(a,l,dX)
% a must have the size :  m x d1 x d x 1 x Mx x My (x Mz)
% Compute da(:,:,:,1,x_1,x_2(,x_3))= d_{x_l} a(:,:,:,1,x_1,x_2(,x_3))
%

siz = size(a);
d=length(siz)-4;
% d=siz(3);
% if ~(length(siz)==4+d && siz(4) == 1 )
%     error('wrong size');
% end
dx=dX(l);

idx='';
for k_dim=2:d
    idx = [idx ',:'];
end

da=zeros(siz);
a = permute(a,[(1:ndims(a)-d) ndims(a)-d+1+mod((l-1):(l-2+d),d)]);
% Put the k coordinate in first place
da = permute(da,[(1:ndims(da)-d) ndims(da)-d+1+mod((l-1):(l-2+d),d)]);
% Put the k coordinate in first place
eval(['da(:,:,:,:,2:end-1' idx ') = 1/(2*dx) * (a(:,:,:,:,3:end' idx ')-a(:,:,:,:,1:end-2' idx '));']);
eval(['da(:,:,:,:,1' idx ') = 1/(dx) * (a(:,:,:,:,2' idx ')-a(:,:,:,:,1' idx '));']);
eval(['da(:,:,:,:,end' idx ') = 1/(dx) * (a(:,:,:,:,end' idx ')-a(:,:,:,:,end-1' idx '));']);
da = permute(da,[(1:ndims(a)-d) ndims(a)-d+1+mod((1-l):(d-l),d)]);
% Put the first coordinate in first place again


% da=zeros(siz);
% if d == 2
%     a = permute(a,[(1:ndims(a)-d) ndims(a)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
%     da = permute(da,[(1:ndims(da)-d) ndims(da)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
% %     a = permute(a,[1 2 3 4 5+(l-1) 5+mod(l,2)]); % Put the x_l coordinate in first place
% %     da = permute(da,[1 2 3 4 5+(l-1) 5+mod(l,2)]); % Put the x_l coordinate in first place
%     da(:,:,:,:,2:end-1,:) = 1/(2*dx) * (a(:,:,:,:,3:end,:)-a(:,:,:,:,1:end-2,:));
%     da(:,:,:,:,1,:) = 1/(dx) * (a(:,:,:,:,2,:)-a(:,:,:,:,1,:));
%     da(:,:,:,:,end,:) = 1/(dx) * (a(:,:,:,:,end,:)-a(:,:,:,:,end-1,:));
%     da = permute(da,[(1:ndims(a)-d) ndims(a)-d+1+mod((1-l):(d-l),d)]); % Put the first coordinate in first place again
% %     da = permute(da,[1 2 3 4 5+mod(1-l,2) 5+mod(2-l,2)]); % Put the first coordinate in first place again
% else
%     a = permute(a,[(1:ndims(a)-d) ndims(a)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
%     da = permute(da,[(1:ndims(da)-d) ndims(da)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
% %     a = permute(a,[1 2 3 4 5+(l-1) 5+mod(l,3) 5+mod(l+1,3)]); % Put the x_l coordinate in first place
% %     da = permute(da,[1 2 3 4 5+(l-1) 5+mod(l,3) 5+mod(l+1,3)]); % Put the x_l coordinate in first place
%     da(:,:,:,:,2:end-1,:,:) = 1/(2*dx) * (a(:,:,:,:,3:end,:,:)-a(:,:,:,:,1:end-2,:,:));
%     da(:,:,:,:,1,:,:) = 1/(dx) * (a(:,:,:,:,2,:,:)-a(:,:,:,:,1,:,:));
%     da(:,:,:,:,end,:,:) = 1/(dx) * (a(:,:,:,:,end,:,:)-a(:,:,:,:,end-1,:,:));
%     da = permute(da,[(1:ndims(a)-d) ndims(a)-d+1+mod((1-l):(d-l),d)]); % Put the first coordinate in first place again
% %     da = permute(da,[1 2 3 4 5+mod(1-l,3) 5+mod(2-l,3) 5+mod(3-l,3)]); % Put the first coordinate in first place again
%
% end