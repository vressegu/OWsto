function filter_exp = design_filter_3(model)
% Design a Gaussian filter for smoothing between forward and backward
% advection
%

%% Get parameters
sigma_filter=model.advection.sigma_filter;

%% Grid
x= model.grid.x;
x = x -mean(x);
y= model.grid.y;
y = y - mean(y);
[x,y]=ndgrid(x,y);

%% Filter
filter_exp = exp( - 1/(2*sigma_filter^2) * (x.^2 + y.^2) );
filter_exp = 1/sum(filter_exp(:)) * filter_exp ;