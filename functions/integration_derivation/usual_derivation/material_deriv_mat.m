function D_f = material_deriv_mat(f,w,dt,dX)
% Compute the material derivative of the vector field f
% We must have size(f) = [Mx My n1 N_t]
% where Mx My are the space dimensions of R^2.
% The result will be of size [ Mx My n1 N_t]
%

%%

f=permute(f,[5 3 4 1 2]); % 1 n1 N_t Mx My 
nabla_f = gradient_mat(f,[dt dX]); % 1 n1 N_t Mx My 3
nabla_f=permute(nabla_f,[4 5 6 2 3 1]); % Mx My 3 n1 N_t 

% Derivative along time
deriv_time = nabla_f(:,:,1,:,:); % Mx My 1 n1 N_t 
deriv_time = permute(deriv_time,[1 2 4 5 3]); % Mx My n1 N_t 

% Spatial gradient
nabla_f = nabla_f(:,:,2:3,:,:); % Mx My 2 n1 N_t 
nabla_f = permute(nabla_f,[1 2 3 5 4]); % Mx My 2 N_t n1
nabla_f = bsxfun(@times, w, nabla_f); % Mx My 2 N_t n1
nabla_f = sum(nabla_f,3); % Mx My 1 N_t n1
nabla_f = permute(nabla_f,[1 2 5 4 3]); % Mx My n1 N_t

% Material derivative
D_f = deriv_time + nabla_f;

end