function da = diff_rob_l(sigma,a,l,dX)
% a must have the size :  m x d1 x d x 1 x Mx x My (x Mz)
% Compute da(:,:,:,1,x_1,x_2(,x_3))= d_{x_l} a(:,:,:,1,x_1,x_2(,x_3))
%

if sigma == 0
    da = diff_l(a,l,dX);
else
    size_filter = max([ceil(2*sigma) 1]);
    idx=-size_filter:size_filter;
    gaussian = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx.^2);
    gaussian = 1/sum(gaussian) * gaussian;
    d_gaussian = -1/(sigma^2) * idx .* gaussian;
    filter = d_gaussian' * gaussian;
    
    siz = size(a);
    d=siz(3);
    if d == 3
        pgaussian = permute(gaussian, [1 3 2]);
        filter = bsxfun(@times, filter, pgaussian); % 5 x 5 x 5
    end
    if ~(length(siz)==4+d && siz(4) == 1 )
        error('wrong size');
    end
    dx=dX(l);
    
    % Get neighborhoods
    s_nei = 2*size_filter+1;% diameter of the neighborhood
    if d==3
        error('pas encore code');
    end
    a=permute(a,[1:3 5:ndims(a) 4]); % m x d1 x d x Mx x My
    neighb = get_neighborhood_mat(a,s_nei); % m x d1 x d x Mx x My x 3 x 3
    
    % Robust derivation using the derivation of gaussian
    if l==1
        filter=permute(filter,[d+(1:(3+d)) 1 2]); % 1 x 1 x 1 x 1 x 1 x 3 x 3
    elseif l== 2
        filter=permute(filter,[d+(1:(3+d)) 2 1]); % 1 x 1 x 1 x 1 x 1 x 3 x 3
    end
    da = bsxfun(@times,neighb, filter); % m x d1 x d x Mx x My x 3 x 3
    da=sum(da,ndims(da)); % m x d1 x d x Mx x My x 3
    da= sum(da,ndims(da)); % m x d1 x d x Mx x My
    da= 1/dx * da;
    
    
    da = permute(da, [1:3 ndims(da)+1 4:ndims(da)]);
    
end

    function neighb = get_neighborhood_mat(A,len)
        % Get the len^2-1 neighborhood (in 2D) of each points and, then, symetrize
        % A must have the size :  ? x ? x ? x Mx x My (2D!!!)
        %
        
        % Neighborhood in first dimension
        s=size(A);
        neighb = zeros([s len]);% ... x Mx x My x len
        
        le=ceil((len-1)/2);
        
        for i=-le:0
            neighb(:,:,:,1:end+i,:,le+1+i)=A(:,:,:,1-i:end,:);
        end
        for i =1:le
            neighb(:,:,:,(1+i):end,:,le+1+i)=A(:,:,:,1:end-i,:);
        end
        
        A=neighb;
        
        % Neighborhood in second dimension
        s=size(A);
        neighb = zeros([s 3]);% ... x Mx x My x len x len
        
        for i=-le:0
            neighb(:,:,:,:,1:end+i,:,le+1+i)=A(:,:,:,:,1-i:end,:);
        end
        for i =1:le
            neighb(:,:,:,:,(1+i):end,:,le+1+i)=A(:,:,:,:,1:end-i,:);
        end
        
    end
end
