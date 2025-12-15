%% Demo AB reaction
clear;
format compact
addpath(fullfile('didactic_tools/'));
beta=14;
gamma=2;
n=2;
ab=@(x,alpha)[...     % ODE r.h.s.
    -x(1)+alpha*(1-x(1))*exp(x(2));...
    -x(2)+beta*alpha*(1-x(1))*exp(x(2))-gamma*x(2)];
%% Track curve of equilibria
f=@(y)ab(y(1:2),y(3));
y0=[0;0;0];
ytan=[0;0;1];
ylist=TrackCurve(f,y0,ytan,'nmax',100,'maxstep',0.05,'print',2,'stop',@(y)y(end)>0.2);
xlist=ylist(1:end-1,:);
plist=ylist(end,:);
%% Plot family of equilibria
clf;
plot(plist(1,:),xlist(1,:),'k-','LineWidth',2,'DisplayName','equilibria')
set(gca,'fontsize',18,'ylim',[0,1]);
xlabel('$\alpha$','Interpreter','latex');
ylabel('concentration $x_1$','Interpreter','latex');
grid on
%% Calculate & plot stability
np=length(plist);
evs=NaN(n,np);
hdev=1e-6;
for k=1:np
    J=Jacobian(@(x)ab(x,plist(k)),xlist(:,k),hdev);
    evs(:,k)=eig(J);
end
nunst=sum(real(evs)>0);
sink=nunst==0;
saddle=nunst==1;
source=nunst==2;
%% Add stability information to plot
hold on
clr=lines();lw3={'LineWidth',3};
plot(plist(sink),xlist(1,sink),'o','Color',clr(1,:),'DisplayName','sink',lw3{:});...
plot(plist(saddle),xlist(1,saddle),'+','Color',clr(2,:),'DisplayName','saddle',lw3{:});
plot(plist(source),xlist(1,source),'d','Color',clr(3,:),'DisplayName','source',lw3{:});
legend('location','best')
