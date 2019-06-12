function plot_dot_obs(model)

hold on,
xi = double(model.obs.x);
n_xi=size(xi,1);

% % max_xi=max(xi);
% % min_xi=min(xi);
% % dx_maille =  100e3; % 100 km
% % marge=ceil(dx_maille/2);
% % idx= (xi > repmat(max_xi - marge,[n_xi,1]) ) & (xi < repmat(min_xi + marge,[n_xi,1]) );
% 
% idx= (xi > repmat(model.kriging.BOX(2,:),[n_xi,1]) ) ...
%     | (xi < repmat(model.kriging.BOX(1,:),[n_xi,1]) );
% idx=any(idx,2);
% xi(idx,:)=[];

DOT_STYLE = {'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4};
plot (xi(:, 1), xi(:, 2), DOT_STYLE{:});
hold off
axis equal