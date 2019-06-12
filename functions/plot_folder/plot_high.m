%% Plots high

load('data/simu2008_30.mat','ssh','dx','dy');

[Mx,My,N]=size(ssh);

dX=[dx dy];
x=dx*(1:Mx);
y=dy*((-My/2+1/2):(+My/2-1/2));
% [x, y] = ndgrid(x,y);

ssh=ssh-mean(mean(mean(ssh)));

L=[min(x) max(x) ; min(y) max(y) ] ;

% Xn = fct_grid_satellite(L,0,200e3,1000e3);
Xn = fct_grid_satellite(L,pi/8);
% Xn = fct_grid(L,dX);

for t=1:N
%     surf(x,y,100*ssh(:,:,t)')
    image(x,y,60*ssh(:,:,t)+30)
    hold on
    plot(Xn(1,:),Xn(2,:),'.');
    hold off
    shading flat;
    axis equal;
    drawnow;
%     pause;
end