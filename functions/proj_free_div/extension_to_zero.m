function [f,MX]=extension_to_zero(f,MX,n)
% Extrapole the function to zero outside the domain
% f : N x M
%
if nargin <3
    n=10;
end

n_add= 2*n;
n=n-2;
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
    
    % Extension to right
    f=extention_right(f,n);
    % Extension to left
    f=f(:,end:-1:1);
    f=extention_right(f,n);
    f=f(:,end:-1:1);
    MX(k)=MX(k)+n_add;
    
    f=reshape(f,[N MX_ MX(k)]);
    f=permute(f,[1 ndims(f) 2:ndims(f)-1]);
    % Put the k coordinate in the first place
    f = permute(f,[1 ndims(f)-d+1+mod((1-k):(d-k),d)]);
    % Put the k coordinate in the right place
end
M=prod(MX);
f=reshape(f,[N M]);

    function f=extention_right(f,n)
        % Extrapolation to right
        %
        
        f0=f(:,end);
        df0=f0-f(:,end-1);
        add = add_f(f0,df0,n);
        Mx = size(f,2);
        f=[f add];
        
        function [add] = add_f(f,df,n)
            % Extrapolation to add on the right
            %
            
            [mu, inv_sigma_2, A] = extension_parameter(f,df,n);
            N0=length(f);
            add =nan(N0,n);
            idx_n=(inv_sigma_2<=1e-10 |  A./f >10 );
            idx_n=(idx_n | isnan(A) | isnan(inv_sigma_2) | isnan(mu) );
            idx_a=(mu>0 & A>10*f);
            for i=1:n
                add(:,i)=A.*exp(-1/2*(i-mu).^2 .* inv_sigma_2 );
                add(find(idx_n),i)=f(find(idx_n))*exp(-3*i/n);
                add(find(idx_a),i)=f(find(idx_a))*(1+sin(3*pi/(2*n)*i));
            end
            
            add = [add zeros(N0,2)];
            
            function [mu, inv_sigma_2, A] = extension_parameter(f,df,n)
                % Parameters of the extrapolation
                %

                df_on_f = df./f;
                epsilon = ones(length(f),1)*norm(f)*1e-3/length(f);
                idx_e=(epsilon>abs(f));
                epsilon(idx_e)=abs(f(idx_e))*1e-1;
                
                inv_sigma_2 = 2/n^2 * ( log(abs(f./epsilon)) + n*df_on_f );
                mu = df_on_f./inv_sigma_2 ;
                A=f.*exp(df_on_f.*mu/2);
                
                idx_d0= (abs(df)<1e-10);
                mu(idx_d0)=0;
                A(idx_d0)=f(idx_d0);
                inv_sigma_2(idx_d0) = 2/n^2 * ( log(abs(f(idx_d0)./epsilon(idx_d0))) );
                
                idx_0= (abs(f)<1e-10) & idx_d0;
                mu(idx_0)=0;
                A(idx_0)=0;
                inv_sigma_2(idx_0)=0;
                
            end
            
        end
    end


end