%% Find Hopf bifurcation in AB reaction
indh=find(nunst==2,1,'last');
yh=ylist(:,indh);
fhopf=@(y)[f(y);...
    trace(dfdx(y))];
yh=Solve(fhopf,yh);
xh=yh(1:2,:);
ph=yh(3,:);
plhopf=plot(ph,xh(1),'k+','linewidth',2,'markerface','k');
legend([plhopf;plfold;pl],[{'Hopf','fold'},eqlabels],'location','best')
