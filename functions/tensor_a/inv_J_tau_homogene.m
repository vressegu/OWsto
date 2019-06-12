function [u_res,p,t,e] = inv_J_tau_homogene(model_input)
% Inverse the operator (f J - 1/(\rho) \tau ) with the RHS g \nabla SSH
%

% Adaptative mesh ?
adaptative_mesh=true;
% Type of boundary conditions
n_cas_bound = 4;
% Normalisation of the equation
normalized_eq = false;

%% Check if a is already computed
% persistent model
%
% if  ( isempty(model)  && isempty( u))
%     update_pers = true;
% else
%     update_pers = ~ ( comp_struct(model_input,model) );
% end
%
% if ~ update_pers
%     u_res = u;
%     return
% end
model=model_input;clear model_input;

%% Defines coefficients of the PDE
% fct_grid_mesh(P,T,U,TIME);
if normalized_eq
    %     coef_a =  model.dt_a ...
    %         * (model.physical_constant.g/f_coriolis(model.coriolis)).^2 ...
    %         * ( exp(model.kriging.lognoisevariance) / (model.l_pixel_nugget/2)^2 ...
    %         - 2 * gamma_cov (0, model,1) );
    invRho_=exp(model.kriging.param(2));
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

% fct_a_pde = @(P,T,u,time) ( coef_a * bsxfun(@times, fct_unity_approx2(fct_grid_mesh(1/coef_x*P,T),model.kriging.BOX) , ...
%     fct_tensor_a_pde_homogene(fct_grid_mesh(1/coef_x*P,T),model) ) );
fct_a_pde = @(P,T,u,time) ( coef_a * fct_tensor_a_pde_homogene(fct_grid_mesh(1/coef_x*P,T),model) );
% % fct_a_pde = @(P,T,U,TIME) ( fct_tensor_a_pde(fct_grid_mesh(P,T,U,TIME),model) );

scalar_A = coef_f * 2 * f_coriolis(model.coriolis) * 1i;

RHS = @(P,T,u,time) ( coef_h * 2 * model.physical_constant.g ...
    * fct_unity_approx2(fct_grid_mesh(1/coef_x*P,T),model.kriging.BOX) ...
    .* multiprod( [1 1i] , - fct_grad_h(fct_grid_mesh(1/coef_x*P,T),model) ) );
% RHS = @(P,T,U,TIME) ( 2* model.physical_constant.rho * model.physical_constant.g ...
%             * multiprod( [1 1i] , fct_grad_h(fct_grid_mesh(P,T,U,TIME),model) ) );

%% Geometry of the domain

% Geometry Description matrix (of the domain) (cf decsg)
gd=[ 3 ; ...                     % Rectangle
    4 ; ...                      % number of line segments in the boundary
    model.kriging.BOX(1,1) ; ...
    model.kriging.BOX(1,1) ; ...
    model.kriging.BOX(2,1) ; ...
    model.kriging.BOX(2,1) ; ... % x-coordinates of the starting points
    %       of the edges
    model.kriging.BOX(:,2) ;
    model.kriging.BOX([2 1],2) ];    % y-coordinates of the starting points
%       of the edges
gd(3:end)=coef_x*gd(3:end);

if csgchk(gd)
    error('wrong Geometry Description matrix');
end

% Decomposed Geometry Matrix (of the domain)
dg=decsg(gd);

% Geometry container
pg = pdeGeometryFromEdges(dg);

%% Test
% if model.plot_PDE
%     % visu
%     figure;
%     pdegplot(dg,'edgeLabels','on')
%     %     pdegplot(dg)
%     %     keyboard;
%     %     warning(['Check que les labels des cotes sont coherents' ...
%     %         'avec les conditions de bords et note ces labels']);
% end

%% Boundary conditions
% The best way would be (on the whole system) periodic boundary
% ( apparently not coded in matlab solver ...) on right
% and left boundaries, and free slip (in complex..) on top and bottom
% boundaries

switch n_cas_bound
    case 1 % 1st try : neuman =0
        % Top and bottom : neuman =0
        es = pg.Edges([1,3]);
        % bc1 = pdeBoundaryConditions(es,'u',C0);% Dirichlet
        Q0=0;
        G0=0;
        bc1 = pdeBoundaryConditions(es,'q',Q0,'g',G0);
        es = pg.Edges([2,4]);
        bc2 = pdeBoundaryConditions(es,'q',Q0,'g',G0);
        
        problem = pde(1); % N is the number of equation in the system
        problem.BoundaryConditions = [bc1,bc2];
        
    case 2 % 2nd try : Dirichlet with deterministic geostrophic balance
        es = pg.Edges([1,2,3,4]);
        bc1 = pdeBoundaryConditions(es,'u',@myufun_inv_J_tau,'Vectorized','on' );% Dirichlet
        
        problem = pde(1); % N is the number of equation in the system
        problem.BoundaryConditions = bc1;
        
    case 3 % 3rd try : Dirichlet with deterministic geostrophic balance and regularization
        es = pg.Edges([1,2,3,4]);
        % myufun_inv_J_tau_reg = @(problem,region,state) (myufun_inv_J_tau_reg_model(problem,region,model));
        bc1 = pdeBoundaryConditions(es,'u',@myufun_inv_J_tau_reg,'Vectorized','on' );% Dirichlet
        % bc1 = pdeBoundaryConditions(e1,'u',@myufun_inv_J_tau_reg,'Vectorized','on' );% Dirichlet
        
        
        problem = pde(1); % N is the number of equation in the system
        problem.BoundaryConditions = bc1;
        
    case 4 % Dirichlet with zero velocity
        es = pg.Edges([1,2,3,4]);
        bc1 = pdeBoundaryConditions(es,'u',0,'Vectorized','on' );% Dirichlet
        
        problem = pde(1); % N is the number of equation in the system
        problem.BoundaryConditions = bc1;
    otherwise
        error('boundary condition not defined');
end
% bound_cond=0; % Boundary Condition Matrix
% warning('boundary condition not defined');

%% Mesh
if ~adaptative_mesh
    % Create a first mesh
    [p,e,t] = initmesh(dg);
    % if model.plot_PDE
    %     % visu
    %     figure;
    %     pdemesh(p,e,t)
    % %     keyboard;
    % end
    
    % Refine the mesh
    for k=1:model.nb_mesh_refin
        [p,e,t] = refinemesh(dg,p,e,t); % refine mesh
    end
    if model.plot_PDE
        % visu
        figure;
        pdemesh(p,e,t)
        %     keyboard;
    end
end

%% Solve the PDE

if adaptative_mesh
    % solve equation in ( L^2 ( R^2 ) + i * L^2 ( R^2 ) )
    [u_res,p,e,t] = adaptmesh_perso(dg,problem ,fct_a_pde , scalar_A , RHS , ...
        'Tripick','pdeadgsc','maxt',5e5,'Ngen',inf,'par',4e-7,'Rmethod','longest');
    %     [u_res,p,e,t] = adaptmesh_perso(dg,problem ,fct_a_pde , scalar_A , RHS , ...
    %         'Tripick','pdeadgsc','maxt',4e5,'Ngen',inf,'par',5e-7,'Rmethod','longest');
    %     [u_res,p,e,t] = adaptmesh_perso(dg,problem ,fct_a_pde , scalar_A , RHS , ...
    %         'Tripick','pdeadworst','maxt',5e5,'Ngen',inf,'par',0.9,'Rmethod','regular');
    %     [u_res,p,e,t] = ...
%                 adaptmesh(dg,problem ,fct_a_pde , scalar_A , RHS ,'maxt',5e5,'Ngen',inf,'par',0.5,'Rmethod','longest');
    %     [u_res,p,e,t] = ...
%                  adaptmesh(dg,problem ,fct_a_pde , scalar_A , RHS ,'maxt',5e5,'Ngen',inf,'par',0.9,'Rmethod','regular');
    %     [u_res,p,e,t] = ...
%                  adaptmesh(dg,problem ,fct_a_pde , scalar_A , RHS ,'maxt',5e5,'Ngen',inf,'par',1e-1,'Rmethod','regular');
    %     [u_res,p,e,t] = ...
%                  adaptmesh(dg,problem ,fct_a_pde , scalar_A , RHS ,'maxt',5e5,'Ngen',20,'par',1e-1,'Rmethod','longest');
else
    u_res = assempde(problem ,p,e,t, fct_a_pde , scalar_A , RHS );
    % U=ASSEMPDE(B,P,E,T,C,A,F)
    %       assembles and solves the (elliptic) PDE problem
    %       -div(c*grad(u))+a*u=f, on a mesh described by P, E, and T,
    %       with boundary conditions given by the function name B.
    %       It eliminates the Dirichlet boundary conditions from the
    %       system of linear equations when solving for U.
end

%% Get back in ( L^2 ( R^2 ) )^2
u_res_r=real(u_res);
u_res_i=imag(u_res);

%% Plot
if model.plot_PDE
    %     figure;
    %     pdesurf(p,t,u_res_r)
    % %     pdesurf(p,t,u_res_r')
    %     keyboard;
    
    if adaptative_mesh
        % visu
        figure;
        pdemesh(p,e,t)
        %     keyboard;
    end
    
    figure;
    pdeplot(p,e,t,'xydata',u_res_r)
    %     keyboard;
    
    figure;
    pdeplot(p,e,t,'xydata',u_res_i)
    %     keyboard;
end

p=1/coef_x*p;

% if model.plot_PDE
%     figure;
%     pdeplot(p,e,t,'xydata',u_res_i)
%     keyboard;
% end

%% Get back in ( L^2 ( R^2 ) )^2
% % u_res(:,2)=imag(u_res);
% % u_res(:,1)=real(u_res(:,1));
% u_res(:,2)=u_res_i; clear u_res_i
% u_res(:,1)=u_res_r; clear u_res_r

%% Express in a grig
% uxy = tri2grid(p,t,u_res,x,y)

%% Save as persistent variable
% u=u_res;

%% Subfonctions
    function u = myufun_inv_J_tau(problem,region,state)
        %     function u = myufun_inv_J_tau_model(problem,region,model)
        % u must
        x=region.x;
        y=region.y;
        pts= stk_dataframe (1/coef_x *[x(:) y(:)]);
        
        % Gradient of ssh in Complex form
        u = multiprod( [1 1i] , coef_h * fct_grad_h(pts,model) );
        
        % Deterministic geostrophic balance
        u = 1i * 1/coef_f * (model.physical_constant.g/f_coriolis(model.coriolis)) * u;
    end

    function u = myufun_inv_J_tau_reg(problem,region,state)
        %     function u = myufun_inv_J_tau_reg_model(problem,region,model)
        
        l_pixel=model.l_pixel_nugget;
        nugget = exp(model.kriging.lognoisevariance);
        lambda = model.dt_a ...
            * (model.physical_constant.g/f_coriolis(model.coriolis)).^2 ...
            * ( nugget / (l_pixel/2)^2 ...
            - 2 * gamma_cov (0, model,1) );
        invRho=exp(model.kriging.param(2));
        lambda = lambda/2 * invRho^2 / f_coriolis(model.coriolis);
        lambda = lambda * coef_a/coef_f;
        
        % Solution of min || f x w - g grad(h) ||^2_2
        x=region.x;
        y=region.y;
        pts= stk_dataframe (1/coef_x *[x(:) y(:)]);
        % Gradient of ssh in Complex form
        u = multiprod( [1 1i] , coef_h * fct_grad_h(pts,model) );
        % Deterministic geostrophic balance
        u = 1i * 1/coef_f * (model.physical_constant.g/f_coriolis(model.coriolis)) * u;
        
        % Solution of min || f x w - g grad(h) ||^2_2 + lambda * || w ||^2_2
        u = 1/(1+lambda) * u;
    end


end
