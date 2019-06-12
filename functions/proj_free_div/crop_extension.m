function [f,MX]=crop_extension(f,MX,n)
% Crop the function to reduce the domain
% f : N x M
%

if nargin <3
    n=10;
end

n_add= 2*n;
N=size(f,1);
d=length(MX);

f=reshape(f,[N MX]);

for k=1:d
    MX_ = MX;
    MX_(k)=[];
    M_=prod(MX_);
    M_=max(M_,1);
    f = permute(f,[1 ndims(f)-d+1+mod((k-1):(k-2+d),d)]);
    % Put the k coordinate in first place
    f=permute(f,[1 3:ndims(f) 2]);
    % Put the k coordinate in last place
    f=reshape(f,[N*M_ MX(k)]);
    
    % Crop
    f=f(:,n+1:end-n);
    MX(k)=MX(k)-n_add;
    
    f=reshape(f,[N MX_ MX(k)]);
    f=permute(f,[1 ndims(f) 2:ndims(f)-1]);
    % Put the k coordinate in the first place
    f = permute(f,[1 ndims(f)-d+1+mod((1-k):(d-k),d)]);
    % Put the k coordinate in the right place
end
M=prod(MX);
f=reshape(f,[N M]);

end