function [u_res,err,tau,RHS] = inv_J_tau_FD_homogene_old(model,XP)
% Inverse the operator (f J - 1/(\rho) \tau ) with the RHS g \nabla SSH
%

M=size(XP,1);
normalized_eq=false;

%% Defines coefficients of the PDE

if normalized_eq
    %     coef_a =  model.dt_a ...
    %         * (model.physical_constant.g/f_coriolis(model.coriolis)).^2 ...
    %         * ( exp(model.kriging.lognoisevariance) / (model.l_pixel_nugget/2)^2 ...
    %         - 2 * gamma_cov (0, model,1) );
    invRho_= exp(model.kriging.param(2));
    coef_a = invRho_^2 / f_coriolis(model.coriolis);
    coef_f = 1/f_coriolis(model.coriolis);
    coef_h = 1 / f_coriolis(model.coriolis);
    coef_x = invRho_;
else
    coef_a = 1;
    coef_f = 1;
    coef_h = 1;
    coef_x = 1;
end

% tau = coef_a * matrix_tau_periodic(model);
tau = coef_a * matrix_tau_homogene(model);
scalar_A = coef_f * f_coriolis(model.coriolis) * 1i;
RHS = ( coef_h * model.physical_constant.g ) ...
    * fct_unity_approx2(XP,model.kriging.BOX)' ...
    .* multiprod( [1 1i] , - fct_grad_h(XP,model) ).';

tau = scalar_A*eye(M) - tau;

%% Solving equation
u_res=linsolve(tau,RHS);

%% Error estimation
err = abs((tau * u_res - RHS)./RHS);

%% Get back in ( L^2 ( R^2 ) )^2
u_res(:,2)=imag(u_res);
u_res(:,1)=real(u_res(:,1));
