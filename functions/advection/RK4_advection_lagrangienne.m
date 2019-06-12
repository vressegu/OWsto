function [X,dX] = RK4_advection_lagrangienne(model, X, w, X0)
% Advection of T with the speed w, in Fourrier space
%
dt=model.advection.dt_adv;
Lx=model.grid.dX(1)*model.grid.MX(1);
Ly=model.grid.dX(2)*model.grid.MX(2);

k1 = deriv_advection_lagrangienne(X, w);
k2 = deriv_advection_lagrangienne( X + k1*dt/2, w);
k3 = deriv_advection_lagrangienne( X + k2*dt/2, w);
k4 = deriv_advection_lagrangienne( X + k3*dt, w);

dX = (dt/3)*(k1/2 + k2 + k3 + k4/2);
X = X + dX;
% X = X + (dt/3)*(k1/2 + k2 + k3 + k4/2);

% Period
X(:,1)= mod(X(:,1),Lx);
X(:,2)= mod(X(:,2),Ly);
% X((X(:,2)<0),2)= 0;
% X((X(:,2)>Ly),2)= Ly;
        
%% Sub function
    function dX_=deriv_advection_lagrangienne(Xt,wt)
        
        dX_(:,1)=interp2(X0{1},X0{2},wt(:,:,1)',Xt(:,1),Xt(:,2));
        dX_(:,2)=interp2(X0{1},X0{2},wt(:,:,2)',Xt(:,1),Xt(:,2));
%         dX_(:,1)=interp2(X0(:,1),X0(:,2),w(:,1),Xt(:,1),Xt(:,2));
%         dX_(:,2)=interp2(X0(:,1),X0(:,2),w(:,2),Xt(:,1),Xt(:,2));
        
%         % Period
%         Xt(:,1)= mod(Xt(:,1),Lx);
%         Xt((Xt(:,2)<0),2)= 0;
%         Xt((Xt(:,2)>Ly),2)= Ly;
    end
end