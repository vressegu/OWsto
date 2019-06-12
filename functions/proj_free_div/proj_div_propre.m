function S = proj_div_propre(S,MX)
% proj_div_propre = (Id - proj_free_div)
% with proj_free_div forcing the constraint div(S) = 0
% Size(S) = M x N x d
% M is the number of point in the space
% N is arbitrary (for instance, number of time step)
% MX = [Mx My] or [ Mx My Mz] is the number of point in respectively the
% axis Ox, Oy and Oz
% d = 2 or 3 is the dimension
%

[M,N,d]=size(S);
bool_real = all(all(all(isreal(S))));

%% Extension
S=reshape(S,[M N*d])';
% Extrapole the function to zero outside the domain in order to be able to
% use Dirichlet conditions
[S,MX]=extension_to_zero(S,MX);
M=prod(MX);
S=reshape(S',[M N d]);

%% Fourrier transform
S = reshape(S,[MX N d]); % Mx x My (x Mz) x N x d
for i=1:d
    S = fft(S,[],i);
end
MX=size(S);
MX=MX(1:d);

%% Projection on the subspace of free divergence matrix

% Spacial frequency k
k=zeros([MX d]);
p=(MX-1)/2;
f_p_tot = floor(p);
for i = 1:d
    f_p=f_p_tot(i);
    MX_rep = MX;
    MX_rep(i)=1;
    if d==2
        k(:,:,i)= repmat( ...
            permute( ...
            [ 0:1/MX(i):f_p/MX(i) -1+( (f_p+1)/MX(i):1/MX(i):(1-1/MX(i)) ) ]...
            ,[ 3:(i+1) 2 1]) ...
                                    , MX_rep);
    else
        k(:,:,:,i)= repmat( permute( ...
            [ 0:1/MX(i):f_p/MX(i) -1+( (f_p+1)/MX(i):1/MX(i):(1-1/MX(i)) ) ]...
            ,[ 3:(i+1) 2 1]) ...
            , MX_rep);
    end
end

% Projection operator A
% A_{i,j} = \delta_{i,j} - k_i k_j / ||k||_2^2
% Mx x My x Mz x d x d
A = bsxfun(@times, k, permute(k,[(1:d) d+2 d+1]) ); % Mx x My x Mz x d x d
norm_k_2 = sum( k.^2, d+1); % Mx x My x Mz x 1
A = bsxfun(@times, 1./norm_k_2, A); % Mx x My x Mz x d x d

% The operator do not modify the constant value
if d==2
    A(1,1,:,:)=zeros(d);
else
    A(1,1,1,:,:)=zeros(d);
end

% Application of the operator
A = permute(A, [1:d d+3 d+1 d+2]); % Mx x My x Mz x 1 x d x d
S = bsxfun( @times, S, A); % Mx x My x Mz x N x d x d
S = squeeze(sum(S, d+2)); % Mx x My x Mz x N x d

%% Inverse Fourrier transform
for i=1:d
    S = ifft(S,[],i);
end
MX=size(S);
MX=MX(1:d);

if bool_real
    S = real(S);
end
S = reshape(S,[M N d]); % M x N x d

%% Crop to remove the extrapolation
S=reshape(S,[M N*d])'; % N*d x M
[S,MX]=crop_extension(S,MX);
M=prod(MX);
S=reshape(S',[M N d]);

