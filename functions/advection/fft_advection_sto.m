function fft_T_adv = fft_advection_sto(model, fft_T, w, sigma_dBT_dt)
% Advection of T with the speed w, in Fourrier space
% Euler?Maruyama method
%

dt=model.advection.dt_adv;

model.advection.step='finite_variation';
vf = deriv_fft_advection2(model, fft_T, w +  sigma_dBT_dt );
fft_T_adv = fft_T + vf*dt;

% model.advection.step='finite_variation';
% vf = deriv_fft_advection2(model, fft_T, w );
% model.advection.step='martingale';
% mart = deriv_fft_advection2(model, fft_T, sigma_dBT_dt);
% fft_T_adv = fft_T + (vf + mart)*dt;

% w=w+sigma_dBT_dt;
% model.advection.step='finite_variation';
% 
% k1 = deriv_fft_advection2(model, fft_T, w);
% k2 = deriv_fft_advection2(model, fft_T + k1*dt/2, w);
% k3 = deriv_fft_advection2(model, fft_T + k2*dt/2, w);
% k4 = deriv_fft_advection2(model, fft_T + k3*dt, w);
% 
% fft_T_adv = fft_T + (dt/3)*(k1/2 + k2 + k3 + k4/2);