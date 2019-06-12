function [f_filtered, var_residu, n_day, residu] = fct_filtered_several_day(f,n_day,dt)
% Temporally filterd the data
%

%% Set parameter
n_dt = n_day*(24*3600)/dt;
% if nargin < 2
%     n_dt = 7;
% end
[Mx,My,Mz,N]=size(f);
MX=[Mx My Mz];
M=prod(MX);

%% Reshape
f = reshape(f,[M N]);
f_filtered = fct_filter(f,n_dt*dt,dt);

%% Subsample
N_d_n_dt = mod(N,n_dt);
if N_d_n_dt > 0
    warning([ num2str(N_d_n_dt) ' time step are removed at the end to filter correctly'])
   f(:,(end + 1 - N_d_n_dt):end)=[];
   f_filtered(:,(end + 1 - N_d_n_dt):end)=[];
   N = size(f,2);
end
f_filtered = reshape(f_filtered,[M n_dt N/n_dt]);
residu = bsxfun(@plus, f, - reshape(f_filtered, [M N]));
f_filtered = f_filtered(:,ceil(n_dt/2),:); 

var_residu = residu.^2;
n_dt_res = max( 4 *3600*24, n_dt*dt);
var_residu = fct_filter(var_residu,n_dt_res,dt);
warning('the residual is averaged on at least 4 days');
var_residu = var_residu(:,ceil(n_dt/2):n_dt:end);

f_filtered = reshape(f_filtered, [MX N/n_dt]);
var_residu = reshape(var_residu, [MX N/n_dt]);
if nargout > 3
    residu = reshape(residu, [MX N]);
end

