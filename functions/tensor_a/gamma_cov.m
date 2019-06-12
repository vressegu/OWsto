function res = gamma_cov(d_2,model,deriv)
% Compute the prior covariance gamma and its derivatives gamma' and gamma''
%

if nargin <3
    deriv=0;
end

cov_type = model.kriging.covariance_type;
param = model.kriging.param;
Sigma2 = exp (param(1));
invRho = exp (param(2));

switch deriv
    case 0
        switch cov_type
            case 'stk_expcov_iso'
                d = sqrt(d_2); clear d_2
                res = Sigma2 * exp ( - invRho * d );
            case 'stk_gausscov_iso'
                res = Sigma2 * exp ( - invRho^2 * d_2 );
            otherwise
                error('not coded yet');
        end
    case 1
        switch cov_type
            case 'stk_expcov_iso'
                warning('not reread');
                l_pixel = model.l_pixel;
                d = sqrt(d_2); clear d_2
                res = - invRho/2 * Sigma2 * ( 1./(d + l_pixel) ) .* exp ( - invRho * d );
                %         res = - invRho/2 * Sigma2 * d_2.^(-1/2) .* exp ( - invRho * sqrt(d_2) );
            case 'stk_gausscov_iso'
                res = - invRho^2 * Sigma2 * exp ( - invRho^2 * d_2 );
            otherwise
                error('not coded yet');
        end
    case 2
        switch cov_type
            case 'stk_expcov_iso'
                 error('not coded yet');
            case 'stk_gausscov_iso'
                res = invRho^4 * Sigma2 * exp ( - invRho^2 * d_2 );
            otherwise
                error('not coded yet');
        end
    case 'fft'
        switch cov_type
            case 'stk_expcov_iso'
                error('not coded yet');
            case 'stk_gausscov_iso'
%                 keyboard;
%                 res = Sigma2 * exp ( log(pi) - 2*log(invRho) - (pi/invRho)^2 * d_2 );
                res = (Sigma2 * pi / invRho^2 ) * exp ( - (pi/invRho)^2 * d_2 );
            otherwise
                error('not coded yet');
        end
    otherwise
        error('the third derivative of the covariance should not be needed');
end


