function df = gradient_rob_mat2(sigma,f,dX)
% Compute the gradient of a matrix field f
% We must have size(f) = [d1 d d_{2+1} ... d_{2+d}]
% where d_{2+1} ... d_{2+d} are the space dimensions of R^d.
% The result will be of size [ d1 d2 d_{2+1} ... d_{2+d} d]
%

if sigma ==0
    df = gradient_mat(f,dX);
else
    
    siz = size(f);
    d=siz(2);
    if ~(length(siz)==2+d)
        error('wrong size');
    end
    % Filter for interrior
    size_filter = max([ceil(2*sigma) 1]);
    idx=-size_filter:size_filter;
    gaussian = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx.^2);
    gaussian = 1/sum(gaussian) * gaussian;
    d_gaussian = -1/(sigma^2) * idx .* gaussian;
    filter = d_gaussian' * gaussian;
    
    % Filters for boundaries
    % Left boundary
    idx_left = idx + 1/2;
    gaussian_left = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx_left.^2);
    gaussian_left = 1/sum(gaussian_left) * gaussian_left;
    d_gaussian_left = -1/(sigma^2) * idx_left .* gaussian_left;
    filter_left = d_gaussian_left' * gaussian;
    % Right boudary
    idx_right = idx - 1/2;
    gaussian_right = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx_right.^2);
    gaussian_right = 1/sum(gaussian_right) * gaussian_right;
    d_gaussian_right = -1/(sigma^2) * idx_right .* gaussian_right;
    filter_right = d_gaussian_right' * gaussian;
    
    if d == 3
        pgaussian = permute(gaussian, [1 3 2]);
        filter = bsxfun(@times, filter, pgaussian); % 3 x 3 x 3
    end
    
    % Get neighborhoods
    s_nei = 2*size_filter+1;% diameter of the neighborhood
    if d==3
        error('pas encore code');
    end
    neighb = get_neighborhood_mat(f,s_nei); % d1 x d x Mx x My x 3 x 3
    
    % Robust derivation using the derivation of gaussian
    df = zeros([size(f) d]);
    
    % First dimension
    % Interior
    filter1=permute(filter,[d+(1:(2+d)) 1 2]); % 1 x 1 x 1 x 1 x 3 x 3
    df1 = bsxfun(@times,neighb(:,:,2:end-1,:,:,:), filter1); % d1 x d x Mx-2 x My x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x Mx-2 x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x Mx-2 x My
    df(:,:,2:end-1,:,1)= 1/dX(1) * df1;
    % Left
    filter1=permute(filter_left,[d+(1:(2+d)) 1 2]); % 1 x 1 x 1 x 1 x 3 x 3
    df1 = bsxfun(@times,neighb(:,:,1,:,:,:), filter1); % d1 x d x 1 x My x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x 1 x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x 1 x My
    df(:,:,1,:,1)= 1/dX(1) * df1;
    % Right
    filter1=permute(filter_right,[d+(1:(2+d)) 1 2]); % 1 x 1 x 1 x 1 x 3 x 3
    df1 = bsxfun(@times,neighb(:,:,end,:,:,:), filter1); % d1 x d x 1 x My x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x 1 x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x 1 x My
    df(:,:,end,:,1)= 1/dX(1) * df1;
    
    % Second dimension
    % Interior
    filter2=permute(filter,[d+(1:(2+d)) 2 1]); % 1 x 1 x 1 x 1 x 3 x 3
    df2 = bsxfun(@times,neighb(:,:,:,2:end-1,:,:), filter2); % d1 x d x Mx x My-2 x 3 x 3
    df2=sum(df2,ndims(df2)); % d1 x d x Mx x My-2 x 3
    df2= sum(df2,ndims(df2)); % d1 x d x Mx x My-2
    df(:,:,:,2:end-1,2)= 1/dX(2) * df2;
    % Left
    filter1=permute(filter_left,[d+(1:(2+d)) 2 1]); % 1 x 1 x 1 x 1 x 3 x 3
    df1 = bsxfun(@times,neighb(:,:,:,1,:,:), filter1); % d1 x d x Mx x 1 x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x 1 x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x 1 x My
    df(:,:,:,1,1)= 1/dX(2) * df1;
    % Right
    filter1=permute(filter_right,[d+(1:(2+d)) 2 1]); % 1 x 1 x 1 x 1 x 3 x 3
    df1 = bsxfun(@times,neighb(:,:,:,end,:,:), filter1); % d1 x d x Mx x 1 x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x 1 x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x 1 x My
    df(:,:,:,end,1)= 1/dX(2) * df1;
    
end

    function neighb = get_neighborhood_mat(A,len)
        % Get the len^2-1 neighborhood (in 2D) of each points and, then, symetrize
        % A must have the size :  ? x ? x Mx x My (2D!!!)
        %
        
        % Neighborhood in first dimension
        s=size(A);
        neighb = zeros([s len]);% ... x Mx x My x len
        
        l=ceil((len-1)/2);
        
        for k=-l:0
            neighb(:,:,1:end+k,:,l+1+k)=A(:,:,1-k:end,:);
        end
        for k =1:l
            neighb(:,:,(1+k):end,:,l+1+k)=A(:,:,1:end-k,:);
        end
        
        A=neighb;
        
        % Neighborhood in second dimension
        s=size(A);
        neighb = zeros([s 3]);% ... x Mx x My x len x len
        
        for k=-l:0
            neighb(:,:,:,1:end+k,:,l+1+k)=A(:,:,:,1-k:end,:);
        end
        for k =1:l
            neighb(:,:,:,(1+k):end,:,l+1+k)=A(:,:,:,1:end-k,:);
        end
    end

end