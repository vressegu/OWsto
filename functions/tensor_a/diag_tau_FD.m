function [D,V] = diag_tau_FD(model)
% Diagonalize the operator  \tau 
%
m_tau = - matrix_tau(model); % symmetric definite positive
m_tr_tau = trace(m_tau);
M=size(m_tau,1);
% nb_eig=ceil(M/1);
% opts=struct();
% opts.issym=true;
% opts.isreal=true;
% sigma='lm';%largest magnitude
% [V,D] = eigs(m_tau,nb_eig,sigma,opts);% eigenvalues and eigenvectors of  (- tau)
% [V,D] = eig(m_tau,'vector');% eigenvalues and eigenvectors of  (- tau)
[V,D,~] = svd(m_tau);% eigenvalues and eigenvectors of  (- tau)
clear m_tau
D=diag(D); % eigenvalues of (-tau)
D=-D; % eigenvalues of tau


%% test
% rate_nrj = -D/m_tr_tau;
rate_nrj = -sum(D)/m_tr_tau

