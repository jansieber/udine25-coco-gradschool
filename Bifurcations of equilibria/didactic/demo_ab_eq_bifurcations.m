%% Equilibrium bifurcations for AB reaction
clear;
demo_ab_eq_stability;
%% find folds
indf=[find(nunst==1,1,'first');find(nunst==2,1,'first')];
yf=ylist(:,indf);
evf=NaN(2,length(indf));
dfdx=@(y)Jacobian(@(x)ab(x,y(3)),y(1:2),hdev);
ffold=@(y)[f(y);det(dfdx(y))];
for i=1:length(indf)
    yf(:,i)=Solve(ffold,yf(:,i));
    evf(:,i)=eig(dfdx(yf(:,i)));
end
format short g
disp('eigenvalues at saddle-node bifurcations');
disp([(1:length(indf))',evf])
xf=yf(1:2,:);
pf=yf(3,:);
plot(pf,xf(1,:),'ko','linewidth',1,'MarkerFaceColor',clr(6,:),'MarkerSize',10,'DisplayName','Fold');
%% Find Hopf bifurcation in AB reaction
indh=find(nunst==2,1,'last');
yh=ylist(:,indh);
fhopf=@(y)[f(y);...
    trace(dfdx(y))];
yh=Solve(fhopf,yh);
xh=yh(1:2,:);
ph=yh(3,:);
plot(ph,xh(1),'kd','linewidth',1,'MarkerFaceColor',clr(5,:),'MarkerSize',10,'DisplayName','Hopf');
evh=eig(dfdx(yh));
disp('eigenvalues at Hopf bifurcation');
disp(evh)
