function f = integration_mat(f,dX,big_data)
% Integration of the function f over the space
% Using a Newton Cotes integration of degree 4 (Booles's rule)
% - f ( n m Mx My (Mz) )is a matrix ( n m ) field sampled on Mx x My x Mz
% points
% - dX is a vector containing the space step along x, y (and z)

if nargin <= 2
    big_data=false;
end
% Get size
sizef = size(f);
% [ n m ] = sizef(1:2);
MX = sizef(3:end);
d = length(MX);
if d ~= size(dX)
    error('wrong size');
end

% zero padding
MX4 = mod(4 - mod(MX-1,4),4);
if d == 2
    f(:,:,MX(1)+(1:MX4(1)), MX(2)+(1:MX4(2)))=0;
else
    f(:,:,MX(1)+(1:MX4(1)), MX(2)+(1:MX4(2)), MX(3)+(1:MX4(3)))=0;
end
MX = MX + MX4;

if any(MX<9)
    error('there are not enought space step');
end
MX = (MX-1)/4-2;

coef_int_begin = 1/90 * [ 7 32 12 32 14];
coef_int = 1/90 * [32 12 32 14];
coef_int_end = 1/90 * [32 12 32 7];

idx='';
for k=1:d
    idx= [ idx ',:' ];
end

N=sizef(1);
for k=d:-1:1
    % for k=1:d
    coef_intk = repmat(coef_int,[1,MX(k)]);
    coef_intk = [coef_int_begin coef_intk coef_int_end];
    vec_permute = [1 3:ndims(f)];
    vec_permute = [vec_permute(1:1+k) 2 vec_permute(2+k:end)];
    coef_intk = permute(coef_intk, vec_permute);
    if big_data
        sizef = size(f);
        res=nan(sizef(1:end-1));
        %         if d == 2
        for t=1:N
            eval(['f(t,:' idx ')=bsxfun(@times,f(t,:' idx '),coef_intk);']);
            eval(['res(t' idx ')=4*dX(k)*sum(f(t,:' idx '),2+k);']);
        end
        %         else
        %             parfor t=1:N
        %                 f(t,:,:,:,:)=bsxfun(@times,f(t,:,:,:,:),coef_intk);
        %                 res(t,:)=4*dX(k)*sum(f(t,:,:,:,:),2+k);
        %             end
        %         end
        f=res;
        clear res coef_intk;
        idx([end-1 end])=[];
    else
        f=bsxfun(@times,f,coef_intk);
        clear coef_intk;
        f=4*dX(k)*sum(f,2+k);
    end
end
