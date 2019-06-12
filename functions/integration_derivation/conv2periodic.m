function c = conv2periodic(a,b)
% compute the 2D convotution assuming that a and b are periodic images
%

% warning('on pourrait integrer seulement sur un quadrant');
if any(size(a)~=size(b))
    error('wrong size');
end

[na,ma]=size(a);
c=zeros(size(a));
for n=1:na
    for m=1:ma
%         c(n:end,m:end) = c(n:end,m:end) +  b(n,m) *...
%             a( 1:(end+1-n),1:(end+1-m)) ;
% %             a( max(1,1-n):min(end,end-n),max(1,1-m):min(end,end-m)) ;
        c = c +  b(n,m) *...
            a( mod((1-n):(end-n),na)+1 ,mod((1-m):(end-m),ma) +1) ;
    end
end
% 
% % %% Sub function
% %     function c = convperiodic(a,b)
% %         % compute the 1D convotution assuming that a and b are periodic
% %         % signals
% %         %
% %         [na,~]=size(a);
% 
% %         [nb,mb]=size(b);
% 
% % c=zeros(size(a));
% % for p=1:na
% %     for q=1:ma
% %         for n=1:na
% %             for m=1:ma
% % %         for n=1:p
% % %             for m=1:q
% %                 c(p,q) = c(p,q) +  b(n,m) *...
% %                     a( mod(p-n,na)+1,mod(q-m,ma)+1) ;
% % %                 c(p,q) = c(p,q) +  b(n,m) *...
% % %                     a( p+1-n,q+1-m) ;
% % % %                 c(p,q) = c(p,q) +  b(n,m) *...
% % % %                     a( p+1-n,q+1-m) ;
% % % %                     a( p-n,q-m) ;
% %             end
% %         end
% % %         for n=p+1:na
% % %             for m=q+1:ma
% % %                 c(p,q) = c(p,q) +  b(n,m) *...
% % %                     a( na+p+1-n,ma+q+1-m) ;
% % % %                     a( p-n,q-m) ;
% % %             end
% % %         end
% %     end
% % end
% 
% % for n=1:na
% %     for m=1:ma
% % %         c(n:end,m:end) = c(n:end,m:end) +  b(n,m) *...
% % %             a( 1:(end+1-n),1:(end+1-m)) ;
% % % %             a( max(1,1-n):min(end,end-n),max(1,1-m):min(end,end-m)) ;
% %         c = c +  b(n,m) *...
% %             a( mod((1-n):(end-n),na)+1 ,mod((1-m):(end-m),ma) +1) ;
% %     end
% % end
% % %         for k=1:na
% % %             c(k,:) = sum a(k-n,:) b(n,:)
% % %         end
% % 
% % %     end
% % % end