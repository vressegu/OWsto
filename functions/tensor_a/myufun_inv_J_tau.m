    function u = myufun_inv_J_tau(problem,region,state)
%     function u = myufun_inv_J_tau_model(problem,region,model)
        % u must 
        x=region.x;
        y=region.y;
        pts=stk_dataframe ([x(:) y(:)]);
        
        % Gradient of ssh in Complex form
        u = multiprod( [1 1i] , fct_grad_h(pts,model) ); 
        
        % Deterministic geostrophic balance
        u = 1i * (model.physical_constant.g/f_coriolis(model.coriolis)) * u;
    end