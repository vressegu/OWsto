function nabla_phi = fct_nabla_phi(model,X)
% This function compute the gradient of the flow
%

Xplot=reshape(X,[model.grid.MX,2]);
% Mx My 2

nabla_phi_x = nan([model.grid.MX,2]);
nabla_phi_y = nan([model.grid.MX,2]);

% Interior
nabla_phi_x(2:end-1,:,:) = 1/(2*model.grid.dX(1)) ...
    *(Xplot(3:end,:,:)-Xplot(1:end-2,:,:));
nabla_phi_y(:,2:end-1,:) = 1/(2*model.grid.dX(2)) ...
    *(Xplot(:,3:end,:)-Xplot(:,1:end-2,:));
% Mx-2 My-2 2 2

% Periodic bourndaries
nabla_phi_x(1,:,:)= 1/(2*model.grid.dX(1)) ...
    *(Xplot(2,:,:)-Xplot(end,:,:));
nabla_phi_x(end,:,:)= 1/(2*model.grid.dX(1)) ...
    *(Xplot(1,:,:)-Xplot(end-1,:,:));
nabla_phi_y(:,1,:) = 1/(2*model.grid.dX(2)) ...
    *(Xplot(:,2,:)-Xplot(:,end,:));
nabla_phi_y(:,end,:) = 1/(2*model.grid.dX(2)) ...
    *(Xplot(:,1,:)-Xplot(:,end-1,:));

% Without boundaries
% nabla_phi_x = 1/(2*model.grid.dX(1)) ...
%     *(Xplot(3:end,2:end-1,:)-Xplot(1:end-2,2:end-1,:));
% nabla_phi_y = 1/(2*model.grid.dX(2)) ...
%     *(Xplot(2:end-1,3:end,:)-Xplot(2:end-1,1:end-2,:));
% % Mx-2 My-2 2 2

nabla_phi(:,:,:,1) = nabla_phi_x;
nabla_phi(:,:,:,2) = nabla_phi_y;

nabla_phi_x = nabla_phi(:,:,1,:);
nabla_phi_y = nabla_phi(:,:,2,:);

iiix = nabla_phi_x >= model.grid.MX(1)/4;
%         figure;imagesc(iiix(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiix(:,:,:,2)');axis xy; axis equal; drawnow;
iiix = iiix(:);
if any(iiix)
    sX = size(nabla_phi_x);
    nabla_phi_x = nabla_phi_x(:);
    nabla_phi_x(iiix) = nabla_phi_x(iiix) - model.grid.MX(1)/2;
    nabla_phi_x = reshape(nabla_phi_x,sX);
end

iiix = nabla_phi_x < -model.grid.MX(1)/4;
%         figure;imagesc(iiix(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiix(:,:,:,2)');axis xy; axis equal; drawnow;
iiix = iiix(:);
if any(iiix)
    sX = size(nabla_phi_x);
    nabla_phi_x = nabla_phi_x(:);
    nabla_phi_x(iiix) = nabla_phi_x(iiix) + model.grid.MX(1)/2;
    nabla_phi_x = reshape(nabla_phi_x,sX);
end

iiiy = nabla_phi_y >= model.grid.MX(2)/4;
%         figure;imagesc(iiiy(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiiy(:,:,:,2)');axis xy; axis equal; drawnow;
iiiy = iiiy(:);
if any(iiiy)
    sX = size(nabla_phi_y);
    nabla_phi_y = nabla_phi_y(:);
    nabla_phi_y(iiiy) = nabla_phi_y(iiiy) - model.grid.MX(2)/2;
    nabla_phi_y = reshape(nabla_phi_y,sX);
end

iiiy = nabla_phi_y < -model.grid.MX(2)/4;
%         figure;imagesc(iiiy(:,:,:,1)');axis xy; axis equal; drawnow;
%         figure;imagesc(iiiy(:,:,:,2)');axis xy; axis equal; drawnow;
iiiy = iiiy(:);
if any(iiiy)
    sX = size(nabla_phi_y);
    nabla_phi_y = nabla_phi_y(:);
    nabla_phi_y(iiiy) = nabla_phi_y(iiiy) + model.grid.MX(2)/2;
    nabla_phi_y = reshape(nabla_phi_y,sX);
end

nabla_phi(:,:,1,:) = nabla_phi_x;
nabla_phi(:,:,2,:) = nabla_phi_y;


% figure;
% subplot(2,2,1);imagesc(nabla_phi(:,:,1,1)');axis xy; axis equal; drawnow;
% subplot(2,2,2);imagesc(nabla_phi(:,:,2,1)');axis xy; axis equal; drawnow;
% subplot(2,2,3);imagesc(nabla_phi(:,:,1,2)');axis xy; axis equal; drawnow;
% subplot(2,2,4);imagesc(nabla_phi(:,:,2,2)');axis xy; axis equal; drawnow;
% keyboard;
