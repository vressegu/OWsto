function df = gradient_rob_mat(sigma,f,dX)
% Compute the gradient of a matrix field f
% We must have size(f) = [d1 d d_{2+1} ... d_{2+d}]
% where d_{2+1} ... d_{2+d} are the space dimensions of R^d.
% The result will be of size [ d1 d2 d_{2+1} ... d_{2+d} d]
%

if sigma ==0
    df = gradient_mat(f,dX);
elseif sigma == -1
    df = gradient_acc_mat(f,dX,'standard');
else
    
    siz = size(f);
    d=siz(2);
    if ~(length(siz)==2+d)
        error('wrong size');
    end
    size_filter = max([ceil(2*sigma) 1]);
    idx=-size_filter:size_filter;
    gaussian = 1/(sqrt(2*pi)*sigma) * exp(-1/(2*sigma^2) * idx.^2);
    gaussian = 1/sum(gaussian) * gaussian;
    d_gaussian = -1/(sigma^2) * idx .* gaussian;
    filter = d_gaussian' * gaussian;
    
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
    df1 = bsxfun(@times,neighb, filter1); % d1 x d x Mx x My x 3 x 3
    df1=sum(df1,ndims(df1)); % d1 x d x Mx x My x 3
    df1= sum(df1,ndims(df1)); % d1 x d x Mx x My
    df(:,:,:,:,1)= 1/dX(1) * df1;
    
    % Second dimension
    filter2=permute(filter,[d+(1:(2+d)) 2 1]); % 1 x 1 x 1 x 1 x 3 x 3
    df2 = bsxfun(@times,neighb, filter2); % d1 x d x Mx x My x 3 x 3
    df2=sum(df2,ndims(df2)); % d1 x d x Mx x My x 3
    df2= sum(df2,ndims(df2)); % d1 x d x Mx x My
    df(:,:,:,:,2)= 1/dX(2) * df2;
    
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


% param_CL={ {'periodic' 'periodic'} , ...
%     {[0 V0] [0 0]} };
% 
% d=size(f,2);
% %% Limits conditions
% for k=1:d
%     CL_temp = param_CL{k}; % CL of the axis k
%     CL_xmin = CL_temp{1};
%     CL_xmax = CL_temp{2};
%     
%     idx='';
%     for k_dim=1:d
%         idx = [idx ':,'];
%     end
%     
% idx(1+2*(k-1))=1;
% 
%     switch CL_xmin
%         case periodic
%             df(:,:,:,:,l)= 1/dX(1) * df1;
%         eval(['df(:,:,' idx 'l)=diff_mat(f,k_dim,dX,type);']);
%     end
% end
% 
