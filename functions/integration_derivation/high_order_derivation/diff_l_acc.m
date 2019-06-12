function da = diff_l_acc(a,l,dX,type)
% a must have the size :  m x d1 x d x 1 x Mx x My (x Mz)
% Compute da(:,:,:,1,x_1,x_2(,x_3))= d_{x_l} a(:,:,:,1,x_1,x_2(,x_3))
%

d1=size(a,2);
d=ndims(a)-4;
% d=size(a,3);
idx='';
for k_dim=2:d
    idx = [idx ',:'];
end

a = permute(a,[(1:ndims(a)-d) ndims(a)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
% m x d1 x d x 1 x Mx x My
a=permute(a,[1:3 5:ndims(a) 4]); % m x d1 x d x Mx x My
da = zeros(size(a)); % m x d1 x d x Mx x My

% [ m, ~,~, Mx, My] = size(a);
size_a = size(a);
m=size_a(1);
Mx = size_a(4);
MY = size_a(5:end); % My (Mz)

dx=dX(l);

if (Mx < 12 && floor((Mx+1)/2)==(Mx+1)/2) || Mx < 2
    error(['There are not enough points along the axis number ' num2str(l) ]);
end


if strcmp(type,'standard')
    % Filter for interrior
    FD13p = [ 1/5544 -3/1155 1/56 -5/63 15/56 -6/7 0 ...
        6/7 -15/56 5/63 -1/56 3/1155 -1/5544 ]';
    FD11p = [ -1/1260 5/504 -5/84 5/21 -5/6 0 5/6 -5/21 5/84 -5/504 1/1260 ]';
    FD9p = [ 1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280 ]';
    FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
    FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
    FD3p = [ -1/2 0 1/2 ]';
    FD2p = [ -1 1]';
elseif strcmp(type,'optimized')
    % Filter for interrior
    
    FD13p = [ +0.001456501759 ...
        -0.011169294114 ...
        +0.045246480208 ...
        -0.133442885327 ...
        +0.337048393268 ...
        -0.907646591371 ...
        0              ...
        +0.907646591371 ...
        -0.337048393268 ...
        +0.133442885327 ...
        -0.045246480208 ...
        +0.011169294114 ...
        -0.001456501759 ]';
    
    FD11p = [ -0.002484594688 ...
        +0.020779405824 ...
        -0.090320001280 ...
        +0.286511173973 ...
        -0.872756993962 ...
        0              ...
        +0.872756993962 ...
        -0.286511173973 ...
        +0.090320001280 ...
        -0.020779405824 ...
        +0.002484594688 ]';
    
    FD9p = [  +0.007650904064 ...
        -0.059463584768 ...
        +0.244678631765 ...
        -0.841570125482 ...
        0               ...
        +0.841570125482 ...
        -0.244678631765 ...
        +0.059463584768 ...
        -0.007650904064 ]';
    
    FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
    FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
    FD3p = [ -1/2 0 1/2 ]';
    FD2p = [ -1 1]';
end
filters={FD2p FD3p FD5p FD7p FD9p FD11p FD13p};

% Boundaries

% Boundaries pixels
neighb1 = zeros([ m d1 d 2 MY 2]);
eval(['neighb1(:,:,:,1' idx ',:) = permute( a(:,:,:,[1 2]' idx ') , ' ...
    ' [ 1:3 ndims(a)+1 5:ndims(a) 4]);']);
eval(['neighb1(:,:,:,2' idx ',:) = permute( a(:,:,:,[end-1 end]' idx '), ' ...
    ' [ 1:3 ndims(a)+1 5:ndims(a) 4]);']);
filter1=permute(filters{1},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x (1 x) len_filter
da1 = bsxfun(@times,neighb1, filter1); % m x d1 x d x 2 x MY x len_filter
clear neighb1;
da1= sum(da1,ndims(da1)); % m x d1 x d x 2 x MY
eval(['da(:,:,:,[1 end]' idx ')= 1/dx * da1;']);

% From 2 to 6 pixels far from the boundaries
for k=2:min(6,Mx/2)
    % Get neighborhoods
    len_filter = 2*(k-1)+1;
    neighb1 = zeros([ m d1 d 2 MY len_filter]);
    eval(['neighb1(:,:,:,1' idx ',:) = permute( a(:,:,:, 1:len_filter ' idx ') , ' ...
        ' [ 1:3 ndims(a)+1 5:ndims(a) 4]);']);
    eval(['neighb1(:,:,:,2' idx ',:) = permute( a(:,:,:, (end-len_filter+1):end ' idx ') , ' ...
        ' [ 1:3 ndims(a)+1 5:ndims(a) 4]);']);
    filter1=permute(filters{k},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x (1 x) len_filter
    da1 = bsxfun(@times,neighb1, filter1); % m x d1 x d x 2 x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1)); % m x d1 x d x 2 x MY
    eval(['da(:,:,:,[k (Mx-(k-1)) ]' idx ')= 1/dx * da1;']);
end

if Mx >= 13
    % Interior
    % Get neighborhoods
    len_filter = 13;
    neighb1 = get_neighborhood_mat_1D( a, len_filter); % m x d1 x d x Mx x MY x len_filter
    eval(['neighb1 = neighb1(:,:,:, 7:(end-6)' idx ',:);']); % m x d1 x d x (Mx-12) x MY x len_filter
    filter1=permute(filters{7},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x (1 x) len_filter
    da1 = bsxfun(@times,neighb1, filter1); % m x d1 x d x (Mx-12) x MY x len_filter
    clear neighb1;
    da1= sum(da1,ndims(da1));  % m x d1 x d x (Mx-12) x MY
    eval(['da(:,:,:,7:(end-6)' idx ')= 1/dx * da1 ;']);
end

da = permute(da, [1:3 ndims(da)+1 4:ndims(da)] ); % m x d1 x d x 1 x Mx x MY
da = permute(da,[(1:ndims(da)-d) ndims(da)-d+1+mod((1-l):(d-l),d)]);
% Put the first coordinate in first place again

end


function neighb = get_neighborhood_mat_1D(A,len)
% Get the len^2-1 neighborhood (in 2D) of each points
% A must have the size :  ? x ? x ? x Mx x My (2D!!!)
%

s=size(A);
neighb = zeros([s len]);% ... x Mx x MY x len
d=s(3);
idx='';
for k_dim=2:d
    idx = [idx ',:'];
end

le=ceil((len-1)/2);

for i=-le:0
    eval(['neighb(:,:,:,1:end+i' idx ',le+1+i)=A(:,:,:,1-i:end' idx ');']);
end
for i =1:le
    eval(['neighb(:,:,:,(1+i):end' idx ',le+1+i)=A(:,:,:,1:end-i' idx ');']);
end

eval(['neighb=neighb(:,:,:,:' idx ',end:-1:1);']);

end


% function da = diff_l_acc(a,l,dX,type)
% % a must have the size :  m x d x d x 1 x Mx x My (x Mz)
% % Compute da(:,:,:,1,x_1,x_2(,x_3))= d_{x_l} a(:,:,:,1,x_1,x_2(,x_3))
% %
%
% d=size(a,2);
% if d==3
%     error('pas encore code');
% end
%
% a = permute(a,[(1:ndims(a)-d) ndims(a)-d+1+mod((l-1):(l-2+d),d)]); % Put the k coordinate in first place
%  % m x d x d x 1 x Mx x My
% a=permute(a,[1:3 5:ndims(a) 4]); % m x d x d x Mx x My
% da = zeros(size(a)); % m x d x d x Mx x My
%
% [ m, ~,~, Mx, My] = size(a);
% dx=dX(l);
%
% if (Mx < 12 && floor((Mx+1)/2)==(Mx+1)/2) || Mx < 2
%     error(['There are not enough points along the axis number ' num2str(l) ]);
% end
%
%
% if strcmp(type,'standard')
%     % Filter for interrior
%     FD13p = [ 1/5544 -3/1155 1/56 -5/63 15/56 -6/7 0 ...
%         6/7 -15/56 5/63 -1/56 3/1155 -1/5544 ]';
%     FD11p = [ -1/1260 5/504 -5/84 5/21 -5/6 0 5/6 -5/21 5/84 -5/504 1/1260 ]';
%     FD9p = [ 1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280 ]';
%     FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
%     FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
%     FD3p = [ -1/2 0 1/2 ]';
%     FD2p = [ -1 1]';
% elseif strcmp(type,'optimized')
%     % Filter for interrior
%
%     FD13p = [ +0.001456501759 ...
%               -0.011169294114 ...
%               +0.045246480208 ...
%               -0.133442885327 ...
%               +0.337048393268 ...
%               -0.907646591371 ...
%                0              ...
%               +0.907646591371 ...
%               -0.337048393268 ...
%               +0.133442885327 ...
%               -0.045246480208 ...
%               +0.011169294114 ...
%               -0.001456501759 ]';
%
%     FD11p = [ -0.002484594688 ...
%               +0.020779405824 ...
%               -0.090320001280 ...
%               +0.286511173973 ...
%               -0.872756993962 ...
%                0              ...
%               +0.872756993962 ...
%               -0.286511173973 ...
%               +0.090320001280 ...
%               -0.020779405824 ...
%               +0.002484594688 ]';
%
%     FD9p = [  +0.007650904064 ...
%               -0.059463584768 ...
%               +0.244678631765 ...
%               -0.841570125482 ...
%               0               ...
%               +0.841570125482 ...
%               -0.244678631765 ...
%               +0.059463584768 ...
%               -0.007650904064 ]';
%
%     FD7p = [ -1/60  3/20 -3/4 0 3/4 -3/20 1/60 ]';
%     FD5p = [ 1/12  -2/3 0 2/3 -1/12 ]';
%     FD3p = [ -1/2 0 1/2 ]';
%     FD2p = [ -1 1]';
% end
% filters={FD2p FD3p FD5p FD7p FD9p FD11p FD13p};
%
% % Boundaries
%
% % Boundaries pixels
% neighb1 = zeros([ m d d 2 My 2]);
% neighb1(:,:,:,1,:,:) = permute( a(:,:,:,[1 2],:) , ...
%                                                 [ 1:3 ndims(a)+1 5:ndims(a) 4]);
% neighb1(:,:,:,2,:,:) = permute( a(:,:,:,[end-1 end],:) , ...
%                                                 [ 1:3 ndims(a)+1 5:ndims(a) 4]);
% filter1=permute(filters{1},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x len_filter
% da1 = bsxfun(@times,neighb1, filter1); % m x d x d x 2 x My x len_filter
% da1= sum(da1,ndims(da1)); % m x d x d x 2 x My
% da(:,:,:,[1 end],:)= 1/dx * da1;
%
% % From 2 to 6 pixels far from the boundaries
% for k=2:min(6,Mx/2)
%     % Get neighborhoods
%     len_filter = 2*(k-1)+1;
%     neighb1 = zeros([ m d d 2 My len_filter]);
%     neighb1(:,:,:,1,:,:) = permute( a(:,:,:, 1:len_filter ,:) , ...
%                                                 [ 1:3 ndims(a)+1 5:ndims(a) 4]);
%     neighb1(:,:,:,2,:,:) = permute( a(:,:,:, (end-len_filter+1):end ,:) , ...
%                                                 [ 1:3 ndims(a)+1 5:ndims(a) 4]);
%     filter1=permute(filters{k},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x len_filter
%     da1 = bsxfun(@times,neighb1, filter1); % m x d x d x 2 x My x len_filter
%     da1= sum(da1,ndims(da1)); % m x d x d x 2 x My
%     da(:,:,:,[k (Mx-(k-1)) ],:)= 1/dx * da1;
% end
%
% if Mx >= 13
%     % Interior
%     % Get neighborhoods
%     len_filter = 13;
%     neighb1 = get_neighborhood_mat_1D( a, len_filter); % m x d x d x Mx x My x len_filter
%     neighb1 = neighb1(:,:,:, 7:(end-6),:,:); % m x d x d x (Mx-12) x My x len_filter
%     filter1=permute(filters{7},[2:ndims(neighb1) 1 ]); % 1 x 1 x 1 x 1 x len_filter
%     da1 = bsxfun(@times,neighb1, filter1); % m x d x d x (Mx-12) x My x len_filter
%     da1= sum(da1,ndims(da1));  % m x d x d x (Mx-12) x My
%     da(:,:,:,7:(end-6),:)= 1/dx * da1 ;
% end
%
% da = permute(da, [1:3 ndims(da)+1 4:ndims(da)] ); % m x d x d x 1 x Mx x My
% da = permute(da,[(1:ndims(da)-d) ndims(da)-d+1+mod((1-l):(d-l),d)]);
% % Put the first coordinate in first place again
%
% end
%
%
% function neighb = get_neighborhood_mat_1D(A,len)
% % Get the len^2-1 neighborhood (in 2D) of each points
% % A must have the size :  ? x ? x ? x Mx x My (2D!!!)
% %
%
% s=size(A);
% neighb = zeros([s len]);% ... x Mx x My x len
%
% le=ceil((len-1)/2);
%
% for i=-le:0
%     neighb(:,:,:,1:end+i,:,le+1+i)=A(:,:,:,1-i:end,:);
% end
% for i =1:le
%     neighb(:,:,:,(1+i):end,:,le+1+i)=A(:,:,:,1:end-i,:);
% end
%
% neighb=neighb(:,:,:,:,:,end:-1:1);
%
% end
