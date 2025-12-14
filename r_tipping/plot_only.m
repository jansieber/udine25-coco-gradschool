%% connecting orbit representing rate-induced tipping
% reference Thoraya Alharti PhD thesis
% <https://ore.exeter.ac.uk/ndownloader/files/56835545>
%%
%#ok<*NASGU>
%% Define r.h.s of ODE and its derivatives using symbolic toolbox 
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
vnames={'x','p','q'};
[ip,npars]=structind_from_names(pnames);
[iv,dim]=structind_from_names(vnames);
id_pars=@(name1,name2,free)struct('match1',name1,'match2',name2,'free',free);
get1=@(x)reshape(x(1,:),[1,size(x,2:ndims(x))]);
% each subsystem (u_gamma, u_plus and u_minus) has a full set of parameters
% their names will be prepended by ug, up and um
name_prep=@(prep,names)cellfun(@(s)[prep,'.',s],names,'UniformOutput',false);
%% plot u_- after re-reading bifurcation diagram
bd_incTm=coco_bd_read('incTm');
incTmlabs=coco_bd_labs(bd_incTm,'EP');
sol_incTm=bvp_read_solution('u_minus','incTm',incTmlabs(2));
txt={'Fontsize',20,'FontName','Courier','FontWeight','bold'};
ltx={'Interpreter','latex'};
lw2={'LineWidth',2};
figure(1);clf;ax=gca;
pminus=plot(ax,sol_incTm{1}.tbp-sol_incTm{1}.T,sol_incTm{1}.xbp,lw2{:});
labminus={'$x_-$','$p_-$','$q_-$'};
legend(ax,pminus,labminus,'Location','best',txt{:},ltx{:});
xlabel('unscaled time $t$',txt{:},ltx{:});
set(ax,txt{:})
drawnow
%% insert resulting Gamma^u into plot
bdg=coco_bd_read('fixUg');
epglabs=coco_bd_labs(bdg,'EP');
solg=bvp_read_solution('u_gamma','fixUg',epglabs(end));
hold(ax,'on');
pgam=plot(ax,sol_incTm{1}.T+solg{1}.tbp,solg{1}.xbp,lw2{:});
labgam={'$x_g$','$p_g$','$q_g$'};
legend(ax,[pminus;pgam],[labminus,labgam],'Location','best');
hold(ax,'off')
drawnow
%% Plot composite solution with gap
bd_incTm=coco_bd_read('incTm');
bd_incTp=coco_bd_read('incTp');
incTplabs=coco_bd_labs(bd_incTp,'EP');
incTmlabs=coco_bd_labs(bd_incTm,'EP');
seglist={'u_gamma','u_plus','u_minus'};
runlist={'incTp','incTp','incTm'};
lablist={incTplabs(end),incTplabs(end),incTmlabs(end)};
figure(1);clf;ax=gca;
plotsol(ax,runlist,lablist,seglist,iv,ip)
drawnow
%% plot solution with gap closed
bd_closegap=coco_bd_read('closegap');
gap0lab=coco_bd_labs(bd_closegap,'GAP');
seglist={'u_gamma','u_plus','u_minus'};
runlist='closegap';
lablist=gap0lab;
include_L=true;
figure(2);clf;ax2=gca;
plotsol(ax2,runlist,lablist,seglist,iv,ip,include_L)
drawnow
%% plot lin gap as function of phi
bd_eta_of_phi=coco_bd_read('eta_of_phi');
phivals=coco_bd_col(bd_eta_of_phi,'um.phi');
etavals=coco_bd_col(bd_eta_of_phi,'eta');
WuPsvals=coco_bd_col(bd_eta_of_phi,'WuPs');
WsGuvals=coco_bd_col(bd_eta_of_phi,'WsGu');
figure(3);clf;
lw={'linewidth',2,'MarkerSize',4};
tiledlayout(2,1,'TileSpacing','tight');
nexttile; %subplot(2,1,1);
ax3=gca;hold(ax3,'on');
peta=plot(ax3,phivals/pi,etavals,'DisplayName','$\eta(\phi)$',lw{:});
plot(ax3,[-1,1],[0,0],'k-');
getval=@(s)get1(coco_bd_col(bd_eta_of_phi,['up.',s])');
title(ax3,sprintf('$r=%6.4f$, $a=%4.2g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    getval('r'),getval('a'),getval('cp'),getval('beta')),'Interpreter','latex')
legend(ax3,peta,'Interpreter','latex');
grid(ax3,'on');
set(ax3,'FontSize',18,'box','on','LineWidth',1,'XTick',[],'YLim',[min(etavals),max(etavals)]*1.1);
nexttile;% subplot(2,1,2);
ax3=gca;hold(ax3,'on');
pwups=plot(ax3,phivals/pi,WuPsvals,'DisplayName','$W^u(P^s,\phi)$',lw{:});
pwsgu=plot(ax3,phivals/pi,WsGuvals,'DisplayName','$W^s(\Gamma^u,\phi)$',lw{:});
xlabel(ax3,'$\phi/\pi$','Interpreter','latex');
legend(ax3,[pwups,pwsgu],'Interpreter','latex','Location','best');
grid(ax3,'on');
set(ax3,'FontSize',18,'box','on','LineWidth',1);
drawnow
%% plot critical rate depending on phi
bd_r_of_phi=coco_bd_table('r_of_phi');
phivals=bd_r_of_phi.('um.phi');
rcritvals=bd_r_of_phi.('um.r');
figure(4);clf;
ax3=gca;hold(ax3,'on');
prcrit=plot(ax3,phivals/pi,rcritvals,'DisplayName','$r_\mathrm{crit}(\phi)$',lw{:});
xlabel(ax3,'$\phi/\pi$','Interpreter','latex');
title(ax3,sprintf('$a=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    bd_r_of_phi{1,'um.a'},bd_r_of_phi{1,'um.cp'},bd_r_of_phi{1,'um.beta'}),ltx{:})
legend(ax3,'Interpreter','latex','Location','best');
grid(ax3,'on');
set(ax3,'FontSize',18,'box','on','LineWidth',1);
drawnow
%% plot critical rate and its approximate derivative depending on phi
bd_r_of_phi_wcopy=coco_bd_table('r_of_phi_wcopy');
phivals=bd_r_of_phi_wcopy.('um.phi');
rvals=bd_r_of_phi_wcopy.('um.r');
drvals=bd_r_of_phi_wcopy.('dr');
figure(4);clf;tiledlayout(2,1,'TileSpacing','tight');
ax4a=nexttile(1);
prcrit=plot(ax4a,phivals/pi,rvals,'DisplayName','$r_\mathrm{crit}(\phi)$',lw{:});
title(ax4a,sprintf('$a=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    bd_r_of_phi_wcopy{1,'um.a'},bd_r_of_phi_wcopy{1,'um.cp'},bd_r_of_phi_wcopy{1,'um.beta'}),ltx{:})
legend(ax4a,'Interpreter','latex','Location','best');
grid(ax4a,'on');
set(ax4a,'FontSize',18,'box','on','LineWidth',1,'xlim',[-1,1],'XTickLabel',{});
ax4b=nexttile(2);
prcrit=plot(ax4b,phivals/pi,drvals,'DisplayName','$r_\mathrm{crit}''(\phi)$',lw{:});
xlabel(ax4b,'$\phi/\pi$','Interpreter','latex');
yline(ax4b,0,'k',lw2{:});
legend(ax4b,'Interpreter','latex','Location','best');
grid(ax4b,'on');
set(ax4b,'FontSize',18,'box','on','LineWidth',1,'xlim',[-1,1]);
drawnow
%% plot critical rates at tangency for phi in a-r plane
clear rvals avals bd_a_r
for k=2:-1:1
    runid=sprintf('a_r_tangency_%d',k);
    bd_a_r{k}=coco_bd_table(runid);
    rvals{k}=bd_a_r{k}{:,'um.r'};
    avals{k}=bd_a_r{k}{:,'um.a'};
end
bd_a_r{3}=coco_bd_table('a_r_phi=0');
rvals{3}=bd_a_r{3}{:,'um.r'};
avals{3}=bd_a_r{3}{:,'um.a'};
kmax=round(1.5+double(rvals{1}(end)<rvals{2}(end))-0.5);
kmin=3-kmax;
%
clr=lines();
figure(5);clf;ax5=gca;hold(ax5,'on');
plot(ax5,rvals{3},avals{3},'DisplayName','critical rate for phase $\phi=0$',lw2{:},'color',clr(5,:));
plot(ax5,rvals{kmax},avals{kmax},'DisplayName','tangency at lower rate (tipping becomes inevitable',lw2{:},'color',clr(2,:));
plot(ax5,rvals{kmin},avals{kmin},'DisplayName','tangency at lower rate (tipping becomes possible',lw2{:},'color',clr(1,:));
ylabel(ax5,'$a$','Interpreter','latex');
xlabel(ax5,'$r$','Interpreter','latex');
title(ax5,sprintf('tipping boundaries for $c_+=%4.2g$, $\\beta=%4.2g$',...
    bd_a_r{1}{1,'um.cp'},bd_a_r{1}{1,'um.beta'}),'Interpreter','latex')
legend(ax5,'Interpreter','latex','Location','best');
grid(ax5,'on');
set(ax5,'FontSize',18,'box','on','LineWidth',2);
drawnow
%% plot solution for smallest r, largest r, r at phi=0
figure(6);clf;
rolab=coco_bd_labs('a_r_phi=0','RO');
for i=1:length(rolab)
    clf
    [plot_L,plot_driver]=deal(true,false);
    plotsol(gca,'a_r_phi=0',rolab(i),seglist,iv,ip,plot_L,plot_driver);
    drawnow
    pause(0.1)
end
