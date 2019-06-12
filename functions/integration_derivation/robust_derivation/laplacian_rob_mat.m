function d2f = laplacian_rob_mat(sigma,f,dX)
% Compute the laplacian of a matrix field f
% We must have size(f) = [d1 d d_{2+1} ... d_{2+d}]
% where d_{2+1} ... d_{2+d} are the space dimensions of R^d.
% The result will be of the same size
%
if sigma == 0
    d2f = laplacian_mat(f,dX);
elseif sigma == -1
    d2f = laplacian_acc_mat(f,dX,'standard');
else
    siz = size(f);
    d=siz(2);
    if ~(length(siz)==2+d)
        error('wrong size');
    end
    
    % Definition of the filter
    size_filter = max([ceil(4*sigma) 1]);
    idx=(-size_filter:1:size_filter)';
    l_idx=length(idx);
    gaussian = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx.^2);
    gaussian = 1/sum(gaussian) * gaussian;
    filter =zeros(l_idx);
    for k=1:d
        d2_gaussian = -1/(dX(k)*sigma)^2 * ( 1 - (1/(sigma) * idx).^2 ) ;
        if k ~=1
            d2_gaussian = permute(d2_gaussian,[2:k 1]);
        end
        vec_rep = l_idx*ones(1,d);
        vec_rep(k) = 1;
        d2_gaussian = repmat(d2_gaussian,vec_rep);
        filter = filter + d2_gaussian;
    end
    for k=1:d
        gaussiank = gaussian;
        if k ~=1
            gaussiank = permute(gaussiank,[2:k 1]);
        end
        filter = bsxfun(@times,filter,gaussiank);
    end
    
    % Get neighborhoods
    s_nei = 2*size_filter+1;% diameter of the neighborhood
    if d==3
        error('pas encore code');
    end
    neighb = get_neighborhood_mat(f,s_nei); % d1 x d x Mx x My x 3 x 3
    
    % Robust derivation using the laplacian of gaussian
    
    filter=permute(filter,[d+(1:(2+d)) 1 2]); % 1 x 1 x 1 x 1 x 3 x 3
    d2f = bsxfun(@times,neighb, filter); % d1 x d x Mx x My x 3 x 3
    for k = 1:d
        d2f=sum(d2f,ndims(d2f)); % d1 x d x Mx x My x 3
    end
    
end

    function neighb = get_neighborhood_mat(A,len)
        % Get the len^2-1 neighborhood (in 2D) of each points and, then, symetrize
        % A must have the size :  ? x ? x Mx x My (2D!!!)
        %
        
        % Neighborhood in first dimension
        s=size(A);
        neighb = zeros([s len]);% ... x Mx x My x len
        
        l=ceil((len-1)/2);
        
        for i=-l:0
            neighb(:,:,1:end+i,:,l+1+i)=A(:,:,1-i:end,:);
        end
        for i =1:l
            neighb(:,:,(1+i):end,:,l+1+i)=A(:,:,1:end-i,:);
        end
        
        A=neighb;
        
        % Neighborhood in second dimension
        s=size(A);
        neighb = zeros([s 3]);% ... x Mx x My x len x len
        
        for i=-l:0
            neighb(:,:,:,1:end+i,:,l+1+i)=A(:,:,:,1-i:end,:);
        end
        for i =1:l
            neighb(:,:,:,(1+i):end,:,l+1+i)=A(:,:,:,1:end-i,:);
        end
    end
end