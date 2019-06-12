    function u = myufun_inv_J_tau_reg(problem,region,state)
%     function u = myufun_inv_J_tau_reg_model(problem,region,model)
        
        l_pixel=model.l_pixel_nugget;
        nugget = exp(model.kriging.lognoisevariance);
        lambda = model.dt_a ...
            * (model.physical_constant.g/f_coriolis(model.coriolis)).^2 ...
            * ( nugget / (l_pixel/2)^2 ...
            - 2 * gamma_cov (0, model,1) );
        invRho=exp(model.kriging.param(2));
        lambda = lambda * invRho^2 / f_coriolis(model.coriolis);
        
        % Solution of min || f x w - g grad(h) ||^2_2
        u = myufun_inv_J_tau_model(problem,region,model);
        
        % Solution of min || f x w - g grad(h) ||^2_2 + lambda * || w ||^2_2
        u = 1/(1+lambda) * u;