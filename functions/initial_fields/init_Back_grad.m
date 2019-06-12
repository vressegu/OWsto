function b_S = init_Back_grad(model,X,Y)
% Generate a 2D buoyancy field with two cold cyclones and two warm
% anticyclones with a slight perturbation in the initial condition
%

Ly = model.grid.dX(2) * model.grid.MX(2);
odg_b = model.odg_b;

b_S = 2/Ly * (Y-mean(Y(:)));
        
% Specify the amplitude of the buoyancy
b_S = odg_b * b_S;