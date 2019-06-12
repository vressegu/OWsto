function on_tau_local = fct_on_tau_forward(model,nabla_phi,w_eul,w_Lag)
% Compute the inverse of the local time scale of differntial advection
%

MX=model.grid.MX;

nabla_phi=permute(nabla_phi,[1 2 4 3]);

% Normalized orthogonal Lagrangian velocity
v_tilde_ortho = fct_ortho(fct_normalize(w_Lag));

% Norm of the Eulerian velocity
n_v = sqrt(sum(w_eul.^2,3));

% Gradient of the norm of the Eulerian velocity
grad_n_v = gradient_perso(model.grid, n_v );

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

on_tau2_global = sum(on_tau_local(:).^2);
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