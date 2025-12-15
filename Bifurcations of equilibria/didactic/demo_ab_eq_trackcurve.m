%% Solve ODE for AB reaction
clear;
addpath(fullfile('didactic_tools/'));
beta=14;
gamma=2;
ab=@(x,alpha)[...     % ODE r.h.s.
    -x(1)+alpha*(1-x(1))*exp(x(2));...
    -x(2)+beta*alpha*(1-x(1))*exp(x(2))-gamma*x(2)];
%% track curve of equilibria
f=@(y)ab(y(1:2),y(3));
y0=[0;0;0];
ytan=[0;0;1];
ylist=TrackCurve(f,y0,ytan,'nmax',100,'maxstep',0.05,'print',2,'stop',@(y)y(3)>0.2);
%%
clf;
plot(ylist(3,:),ylist(1,:),'.-')
set(gca,'fontweight','bold','ylim',[0,1]);
xlabel('alpha');
ylabel('concentration');
grid on