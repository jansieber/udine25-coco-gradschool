%% plot growth histories along periodic orbit
%
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities']);
format compact
load('results_daphnia.mat')
%%
figure(1);clf;ax1=gca;hold(ax1,'on');
bval=@(x)par0(ip.sigma).*x(ix.r,:)./(1+par0(ip.sigma).*x(ix.r,:)).*x(ix.bd,:); %unscaled birth rate

Plot2dBranch(pos_eqs_fine,'funcs',funcs,'ax',ax1,'y',@(p)bval(p.x));
Plot2dBranch(psolbr,'funcs',funcs,'ax',ax1,'y',@(p)max(bval(p.profile),[],2));
xlabel(ax1,'mortality $\mu$','Interpreter','latex');
ylabel(ax1,'(max) birth rate $b$','Interpreter','latex');
legend(ax1,'Location','EastOutSide');
ylim(ax1,[0,4e-3])
%%
figure(2);clf;ax2=gca;
Plot2dBranch(hopf_fine,'lgname','first Hopf bifurcation');
hold(ax2,'on')
Plot2dBranch(transcr,'lgname','Tanscritical bifurcation');
Plot2dBranch(mbranch,'parameter',[ip.mu,ip.cap],'lgname',...
    sprintf('Singularity\n(max of maturation age=%4.1f)\n',max_r0*par0(ip.amax)));
Plot2dBranch(pofold,'lgname','Fold of periodic orbits','funcs',C1info(1).funcs);
mutr=getpar(transcr,ip.mu);
xlabel(ax2,'mortality $\mu$','Interpreter','latex');
ylabel(ax2,'resource carrying capacity $C$','Interpreter','latex');
legend(ax2,'location','southeast');
%%
pt=p0;
[yarr,istruc]=dde_chain_arrays(tab_dist,'chain_delay',pt);
size_arr=yarr.profile.sd.int+exp(-yarr.mesh.tau*pt.parameter(ip.growth));
%% plotting
fn={'fontsize',14};
lw={'LineWidth',2};
figure(3);tiledlayout(2,2);nexttile;ax1=gca;
clr=colormap('default');
contourf(ax1,yarr.mesh.time,yarr.mesh.tau,size_arr);
colorbar(ax1);
hold(ax1,'on');
szA=pt.parameter(ip.szA);
contour(ax1,yarr.mesh.time,yarr.mesh.tau,size_arr,[szA,szA],'m','linewidth',2)
title(ax1,'Individuals'' size by age over time');
set(ax1,fn{:},lw{:});
nexttile;ax2=gca;
contourf(ax2,yarr.mesh.time,yarr.mesh.tau,yarr.profile.ceff.int);
colorbar(ax2);
set(ax2,fn{:},lw{:});
title(ax2,'cumulative population effect');
nexttile;
[ax3,h1,h2]=plotyy(pt.mesh*pt.period,pt.profile([ix.r,ix.raA],:),...
    pt.mesh*pt.period,pt.profile(ix.bd,:));axis(ax3,'tight');
%hold(ax3(1),'on')
%plot(ax3(1),pt.parameter(ipe.t0)*pt.period,e_raA'*dde_coll_eva(pt,pt.parameter(ipe.t0)),'ko',...
%    'markerfacecolor','k','DisplayName','max. neg. slope')
set(h1(1),'LineWidth',2,'DisplayName','$r(t)$');
set(h1(2),'LineWidth',2,'DisplayName','$a_A(t)/a_{\max}$');
set(h2,'DisplayName','birth rate','LineWidth',2);
clr=lines();
h2.Color=clr(4,:);
set(ax3(2),'YColor',clr(4,:),fn{:},lw{:},'ylim',[0,max(pt.profile(ix.bd,:))*1.1]);
set(ax3(1),fn{:},lw{:},'ylim',[0,1]);
legend(ax3(2),'Interpreter','latex');
nexttile;ax4=gca;
contourf(ax4,yarr.mesh.time,yarr.mesh.tau,log10(yarr.profile.ceff.id));
colorbar(ax4);
title(ax4,'log10 population effect density');
set(ax4,fn{:},lw{:});
axis(ax4,'equal')