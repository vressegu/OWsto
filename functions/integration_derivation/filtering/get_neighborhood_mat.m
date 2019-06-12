function neighb = get_neighborhood_mat(A,len)
% Get the len^2-1 neighborhood (in 2D) of each points and, then, symetrize
% A must have the size :  ? x ? x Mx x My (2D!!!)
%

% Neighborhood in first dimension
s=size(A);
neighb = zeros([s len]);% ... x Mx x My x len

l=ceil((len-1)/2);

for i=-l:0
    neighb(:,:,1:end+i,:,l+1+i)=A(:,:,1-i:end,:);% ,<-- partie a changer a chque fois selon la taille de A
end
for i =1:l
    neighb(:,:,(1+i):end,:,l+1+i)=A(:,:,1:end-i,:);% ,<-- partie a changer a chque fois selon la taille de A
end
% neighb = neighb(:,:,:,:,[3 2 1]);

A=neighb;

% Neighborhood in second dimension
s=size(A);
neighb = zeros([s 3]);% ... x Mx x My x len x len

for i=-l:0
    neighb(:,:,:,1:end+i,:,l+1+i)=A(:,:,:,1-i:end,:);% ,<-- partie a changer a chque fois selon la taille de A
end
for i =1:l
    neighb(:,:,:,(1+i):end,:,l+1+i)=A(:,:,:,1:end-i,:);% ,<-- partie a changer a chque fois selon la taille de A
end
% neighb = neighb(:,:,:,:,:,[3 2 1]);
end