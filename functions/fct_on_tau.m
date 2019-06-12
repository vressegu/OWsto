function on_tau_local = fct_on_tau(model,nabla_phi,w_eul,w_Lag)
% Compute the inverse of the local time scale of differntial advection
%

MX=model.grid.MX;
x=model.grid.dX(1)* (0:(MX(1)-1)) ;
x=x-mean(x);
y=model.grid.dX(2)* (0:(MX(2)-1)) ;
y=y-mean(y);
[x,y]=ndgrid(x,y);
r=sqrt(x.^2+y.^2);

MX=model.grid.MX;
w_Lag=permute(w_Lag,[1 2 4 3]);
nabla_phi=permute(nabla_phi,[1 2 4 3]);

% figure;
% norm_nabla_phi = sum(sum(nabla_phi.^2,4),3);
% imagesc(norm_nabla_phi);
% % imagesc(1./norm_nabla_phi);

% Normalized orthogonal Lagrangian velocity
v_tilde_ortho = fct_ortho(fct_normalize(w_Lag));

% % figure;
% % norm = sum(v_tilde_ortho.^2,3);
% % % norm = sum(v_tilde_ortho.^2,3);
% % imagesc(norm);
% figure;
% quiver(v_tilde_ortho(:,:,1)',v_tilde_ortho(:,:,2)'); 
% axis xy; axis equal

% Norm of the Eulerian velocity
n_v = sqrt(sum(w_eul.^2,3));
% n_v = n_v ./ r;

% figure;
% imagesc(n_v)

% Curvature
w_eul_n = fct_normalize(w_eul);
w_eul_n_ortho = fct_ortho(w_eul_n);
% Gradient of the norm of the Eulerian velocity
grad_w_eul_n = permute( w_eul_n , [1 2 4 3]);
grad_w_eul_n = gradient_perso(model.grid, grad_w_eul_n );
grad_w_eul_n = bsxfun(@times, w_eul_n, grad_w_eul_n);
grad_w_eul_n = sum(grad_w_eul_n,3);
grad_w_eul_n = permute( grad_w_eul_n , [1 2 4 3]);
curv = bsxfun(@times, w_eul_n_ortho, grad_w_eul_n);
curv = sum(curv,3);

% figure;
% % imagesc(curv);
% % imagesc(curv.*r);
% mask = ( n_v > 1e-3 * max(n_v(:)) );
% % imagesc(mask);
% imagesc(curv.*r.*mask);

% Using local angular velocity instead of velocity
n_v = n_v .* curv;

% Gradient of the norm of the Eulerian velocity
grad_n_v = gradient_perso(model.grid, n_v );

% figure;
% norm = sum(grad_n_v.^2,3);
% imagesc(norm);

% Remove boundary pixels
v_tilde_ortho=v_tilde_ortho(2:MX(1)-1,2:MX(2)-1,:);
grad_n_v=grad_n_v(2:MX(1)-1,2:MX(2)-1,:);
MX=MX-2;

% Inversion of nabla_phi
on_tau_local = nan([2 MX]);
% on_tau_local = nan([MX 2]);
grad_n_v=permute(grad_n_v,[3 1 2]);
nabla_phi=permute(nabla_phi,[3 4 1 2]);
for i=1:MX(1)
    for j=1:MX(2)
        on_tau_local(:,i,j) = nabla_phi(:,:,i,j) \ grad_n_v(:,i,j);
%         on_tau_local(i,j,:) = nabla_phi(i,j,:,:) \ grad_n_v(i,j,:);
    end
end
on_tau_local=permute(on_tau_local,[2 3 1]);

% Derivation at the initial time othogonal to the stream direction
on_tau_local = sum( v_tilde_ortho .* on_tau_local , 3) ;

figure(22)
imagesc(on_tau_local.^2)

on_tau2_global = 1/prod(model.grid.MX) * sum(on_tau_local(:).^2);
tau_global = 1/sqrt(on_tau2_global) /(3600*24)
% keyboard;

end

function f_ortho = fct_ortho(f)
% Compute the orthogonal vector in each point of the space
%
f_ortho(:,:,1)= - f(:,:,2);
f_ortho(:,:,2)= + f(:,:,1);
end

function g = fct_normalize(g)
% Compute the orthogonal vector in each point of the space
%
ng = sqrt(sum(g.^2,3));
g = bsxfun( @times, 1./ng , g);
end