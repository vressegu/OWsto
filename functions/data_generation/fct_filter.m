function f = fct_filter(f,sigma,dt)
% Filter f along second dimension of f 
%

N_ini = size(f,2);
sigma = 1/8 * sigma;

%% Mirror (since the data are non-periodic)
f = [f(:,end:-1:1) f];
N = size(f,2);

% Ensure that the number of time step is even
N_d_n_dt = mod(N,2);
if N_d_n_dt > 0
    warning([ num2str(N_d_n_dt) ' time step are removed at the end to filter correctly'])
   f(:,(end + 1 - N_d_n_dt):end)=[];
   N = size(f,2);
end
f = fft(f,[],2);

%% Filter in temporal Fourier space
P=N/2;
freq =1/N*[ 0:(P-1) 0 (1-P):-1] ;
freq=2*pi/dt*freq;
freq(P+1)=0;

filter = exp( - 1/2 * (2*pi*sigma)^2 * freq.^2 );
filter(P+1)=0;

f = bsxfun(@times, f , filter);
f = real(ifft(f,[],2));

%% Remove irror (since the data are non-periodic)
f = f(:,(N_ini+1):(2*N_ini));

