function c_k_M = fct_c_k_M(a0,k_m,k_M,beta)
% Compute the factor c_{kappa_M} which characterize the statististics of 
% grad (sigma) dBt

d = size(a0,1);
a0 = trace(a0)/d;

if beta ~= 3
    c_inf =  a0 * k_m^2 / (beta-3);
    c_k_M = 1 - (beta^2-1)/8 * (k_M/k_m).^(3-beta);
    c_k_M = c_inf * c_k_M ;
else
    c_k_M =  a0 * k_m^2 * ...
        ( log(k_M/k_m) - 3/4 );
end