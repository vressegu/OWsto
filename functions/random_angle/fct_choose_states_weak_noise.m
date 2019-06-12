function state = fct_choose_states_weak_noise(r)

r =abs(r);

inf_hyp = 0.5;
inf_extrem = 0.1;
% inf_sto_hyp = 0.5;
sup_ell = 0.1;
% sup_ell = 0.01;

% state = 0 * ...
%     ( mean_sin_v_modif < sup_ell & m4_sin_v < inf_extrem  ) ...
%     ... % deterministic ellipticity
%       + 1 * ...
%       (  mean_sin_v_modif < sup_ell & m4_sin_v >= inf_extrem ) ...
%       ... % ellipticity due to stochastic elliptic bursts
%       + 2 * ...
%       ( mean_sin_v_modif >= sup_ell & mean_sin_v_modif <= inf_hyp & m4_sin_v < inf_extrem  ) ...
%       ... % weak stochastic hyperbolity
%       + 3 * ...
%       (  mean_sin_v_modif >= sup_ell & mean_sin_v_modif <= inf_hyp & m4_sin_v >= inf_extrem ) ...
%       ... % weaker hyperbolity due to stochastic elliptic bursts
%       + 4 * ...
%       ( mean_sin_v_modif >= inf_hyp  & m4_sin_v >= inf_extrem & 10.^(alpha2)/12 <= m_sin_zeta) ...
%       ... % deterministic hyperbolicity with stochastic elliptic bursts
%       + 5 * ...
%       ( mean_sin_v_modif >= inf_hyp  & m4_sin_v < inf_extrem & 10.^(alpha2)/12 <= m_sin_zeta) ...
%       ... % deterministic hyperbolicity
%       + 6 * ...
%       ( mean_sin_v_modif >= inf_hyp  & 10.^(alpha2)/12 > m_sin_zeta); 
%           % stochastic hyperbolicity



state = 0 * ...
    ( r >= 1 ) ...
    ... % deterministic ellipticity
      + 4 * ...
    ( r < 1 ) ;
        % deterministic hyperbolicity
        