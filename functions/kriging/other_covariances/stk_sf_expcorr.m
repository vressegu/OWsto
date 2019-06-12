% STK_SF_EXPCORR computes the Exp correlation function
%
% CALL: K = stk_sf_expcorr (H)
%
%    computes the value of the Gaussian correlation function at distance H.
%
% CALL: K = stk_sf_expcorr (H, DIFF)
%
%    computes the derivative  of the Gaussian correlation function  with respect
%    to the distance H  if DIFF is equal to 1.  If DIFF is equal to -1,  this is
%    the same as K = stk_sf_expcorr (H).
%
% NOTES:
%
% See also: stk_sf_gausscor


function k = stk_sf_expcorr (h, diff)

if nargin > 2,
    stk_error ('Too many input arguments.', 'TooManyInputArgs');
end

% default: compute the value (not a derivative)
if nargin < 2,
    diff = -1;
end

switch diff,
    
    case -1, % value of the covariance function
        
        k = exp (- h );
        
    case 1, % derivative wrt h
        
        k = - exp (- h );
        
    otherwise
        
        error ('incorrect value for diff.');
        
end

end % function stk_sf_expcorr

