function tau_f = tau_periodic_sparse_by_mat(model,f)
% Compute tau(f) on the grid
% with periodic boudaries conditions on null conditions on y=0
%

model.grid.MX=size(f);
f(:,1)=[];
tau_f = matrix_tau_periodic_sparse2(model) * f(:); clear f
tau_f = reshape(tau_f, model.grid.MX - [0 1]);
tau_f = [ zeros(model.grid.MX(1),1) tau_f];