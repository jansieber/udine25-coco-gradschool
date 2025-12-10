%% plots of outcomes from computations for connecting orbit 
% representing rate-induced tipping
% reference Thoraya Alharti PhD thesis
% <https://ore.exeter.ac.uk/ndownloader/files/56835545>
%%
%#ok<*NASGU>
clear
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
[ip,npars]=structind_from_names(pnames);
get1=@(x)reshape(x(1,:),[1,size(x,2:ndims(x))]);
%% plot initial u_- after re-reading bifurcation diagram
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
%% insert Gamma^u into initial plot
bdg=coco_bd_read('fixUg');
epglabs=coco_bd_labs(bdg,'EP');
solg=bvp_read_solution('u_gamma','fixUg',epglabs(end));
hold(ax,'on');
pgam=plot(ax,sol_incTm{1}.T+solg{1}.tbp,solg{1}.xbp);
labgam={'$x_g$','$p_g$','$q_g$'};
legend(ax,[pminus;pgam],[labminus,labgam],'Location','best');
hold(ax,'off')
drawnow
%% Plot initial composite solution with gap
bd_incTm=coco_bd_read('incTm');
bd_incTp=coco_bd_read('incTp');
incTplabs=coco_bd_labs(bd_incTp,'EP');
incTmlabs=coco_bd_labs(bd_incTm,'EP');
seglist={'u_gamma','u_plus','u_minus'};
runlist={'incTp','incTp','incTm'};
lablist={incTplabs(end),incTplabs(end),incTmlabs(end)};
figure(1);clf;ax=gca;
plotsol(ax,runlist,lablist,seglist,ip)
drawnow
%% plot solution with gap closed
bd_closegap=coco_bd_read('closegap');
gap0lab=coco_bd_labs(bd_closegap,'GAP');
seglist={'u_gamma','u_plus','u_minus'};
runlist='closegap';
lablist=gap0lab;
figure(2);clf;ax2
plotsol(ax2,runlist,lablist,seglist,ip)
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
bd_r_of_phi=coco_bd_read('r_of_phi');
phivals=coco_bd_col(bd_r_of_phi,'um.phi');
rcritvals=coco_bd_col(bd_r_of_phi,'um.r');
figure(4);clf;
ax3=gca;hold(ax3,'on');
prcrit=plot(ax3,phivals/pi,rcritvals,'DisplayName','$r_\mathrm{crit}(\phi)$',lw{:});
xlabel(ax3,'$\phi/\pi$','Interpreter','latex');
getval=@(s)get1(coco_bd_col(bd_r_of_phi,['up.',s])');
title(ax3,sprintf('$a=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    getval('a'),getval('cp'),getval('beta')),'Interpreter','latex')
legend(ax3,'Interpreter','latex','Location','best');
grid(ax3,'on');
set(ax3,'FontSize',18,'box','on','LineWidth',1);
drawnow
%% plot critical rate at phi=0 for different a
bd_a_r=coco_bd_read('a_r_plane');
rvals=coco_bd_col(bd_a_r,'um.r');
avals=coco_bd_col(bd_a_r,'um.a');
figure(5);clf;ax4=gca;hold(ax4,'on');
plot(ax4,avals,rvals,'DisplayName','$\eta=0$ for $\phi=0$',lw{:});
xlabel(ax4,'$a$','Interpreter','latex');
ylabel(ax4,'$r$','Interpreter','latex');
getval=@(s)get1(coco_bd_col(bd_a_r,['um.',s])');
title(ax4,sprintf('$\\phi=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    getval('phi'),getval('cp'),getval('beta')),'Interpreter','latex')
legend(ax4,'Interpreter','latex','Location','best');
grid(ax4,'on');
set(ax4,'FontSize',18,'box','on','LineWidth',2);
drawnow
%% plot solution for smallest r
figure(6);clf;
uzarlab=coco_bd_labs(bd_a_r,'EP');
plotsol(gca,'a_r_plane',uzarlab(1),seglist,ip)
drawnow
%% plot composite solution
function plotsol(ax,run,lab,seglist,ip)
vnames={'x','p','q'};
segnames={'\gamma','+','-'};
if ~iscell(run)
    run=repmat({run},1,length(seglist));
    lab=repmat({lab},1,length(seglist));
end
for i=length(seglist):-1:1
    sol(i)=bvp_read_solution(seglist{i},run{i},lab{i});
end
t_offset=[sol{2}.T,0,-sol{3}.T];
lw={'linewidth',2,'MarkerSize',4};
lstyle={':','-','-'};
ish=ishold(ax);
hold(ax,'on');
for i=1:length(seglist)
    for k=1:length(vnames)
        plot(ax,sol{i}.tbp+t_offset(i),sol{i}.xbp(:,k),lstyle{i},...
            'DisplayName',['$',vnames{k},'_',segnames{i},'$'],lw{:});
    end
end
yrg=[min(sol{1}.xbp(:)),max(sol{1}.xbp(:))]*1.1;
plot(ax,[0,0],yrg,'k-','DisplayName','Lin section',lw{:});
plot(ax,t_offset(1)+[0,0],yrg,'k:','DisplayName','$T_+$',lw{:});
plot(ax,t_offset(3)+[0,0],yrg,'k--','DisplayName','$T_-$',lw{:});
legend(ax,'Interpreter','LaTeX','Location','eastoutside');
title(ax,sprintf('$r=%6.4g$, $a=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    sol{1}.p(ip.r),sol{1}.p(ip.a),sol{1}.p(ip.cp),sol{1}.p(ip.beta)),'Interpreter','latex')
set(ax,'FontSize',18,'box','on','LineWidth',2);
xlabel(ax,'time $t$','Interpreter','latex');
ax.XLim(1)=t_offset(3)*1.05;
ylim(yrg);
if ~ish
    hold(ax,'off');
end
end