% STK_EXPCOV_ISO computes the isotropic Exponential covariance function
%
% CALL: K = stk_expcov_iso (PARAM, X, Y)
%
%   computes  the covariance matrix K  between  the sets of locations  X  and Y, 
%   using the isotropic Exponential covariance function with parameters PARAM.  The
%   output matrix K has size NX x NY, where NX is the number of rows in X and NY
%   the number of rows in Y. The vector of parameters must have two elements:
%
%     * PARAM(1) = log (SIGMA ^ 2), where SIGMA is the standard deviation,
%
%     * PARAM(2) = - log (RHO), where RHO is the range parameter.
%
% CALL: dK = stk_expcov_iso (PARAM, X, Y, DIFF)
%
%   computes the derivative of the covariance matrix with respect to PARAM(DIFF)
%   if DIFF is equal to 1 or 2, or the covariance matrix itself if DIFF is equal
%   to -1 (in which case this is equivalent to stk_gausscov_iso (PARAM, X, Y)).
%
% CALL: K = stk_expcov_iso (PARAM, X, Y, DIFF, PAIRWISE)
%
%   computes the covariance vector  (or a derivative of it if DIFF > 0)  between
%   the sets of locations X and Y.  The output K is a vector of length N,  where
%   N is the common number of rows of X and Y.
%
% See also: stk_gausscov_iso


function k = stk_expcov_iso (param, x, y, diff, pairwise)

if nargin > 5,
    stk_error ('Too many input arguments.', 'TooManyInputArgs');
end

persistent x0 y0 param0 pairwise0 D

% process input arguments
x = double (x);
y = double (y);
if nargin < 4, diff = -1; end
if nargin < 5, pairwise = false; end

% extract parameters from the "param" vector
Sigma2 = exp (param(1));
invRho = exp (param(2));

% check parameter values
if ~ (Sigma2 > 0) || ~ (invRho > 0),
    error ('Incorrect parameter value.');
end

% check if all input arguments are the same as before
% (or if this is the first call to the function)
if isempty (x0) || isempty (y0) || isempty (param0) || ...
        ~ isequal ({x, y, param}, {x0, y0, param0}) || ...
        ~ isequal (pairwise, pairwise0)
    % compute the distance matrix
    D  = invRho * stk_dist (x, y, pairwise);
    % save arguments for the nex call
    x0 = x;  y0 = y;  param0 = param;  pairwise0 = pairwise;
end

if diff == -1,
    % compute the value (not a derivative)
    k = Sigma2 * stk_sf_expcorr (D, -1);
elseif diff == 1,
    % diff wrt param(1) = log(Sigma2)
    k = Sigma2 * stk_sf_expcorr (D, -1);
elseif diff == 2,
    % diff wrt param(3) = log(invRho)
    k = D .* (Sigma2 * stk_sf_expcorr (D, 1));
else
    error('there must be something wrong here !');
end

end % function stk_expcov_iso


