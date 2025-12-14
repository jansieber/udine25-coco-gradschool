function [yarr,istruc]=contour_daphnia(fg,tab,ip,ix,pt)
[yarr,istruc]=dde_chain_arrays(tab,'chain_delay',pt);
size_arr=yarr.profile.sd.int+exp(-yarr.mesh.tau*pt.parameter(ip.growth));
%% plotting
fn={'fontsize',14};
lw={'LineWidth',2};
figure(fg);tiledlayout(2,2);nexttile;ax1=gca;
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
tfine=linspace(0,1,10000);
yfine=dde_coll_eva(pt,tfine);
[ax3,h1,h2]=plotyy(tfine*pt.period,...
    [yfine(ix.r,:)/pt.parameter(ip.cap);yfine(ix.raA,:)],...
    tfine*pt.period,[yfine(ix.bd,:);yfine(ix.cd,:)]);axis(ax3,'tight');
%hold(ax3(1),'on')
%plot(ax3(1),pt.parameter(ipe.t0)*pt.period,e_raA'*dde_coll_eva(pt,pt.parameter(ipe.t0)),'ko',...
%    'markerfacecolor','k','DisplayName','max. neg. slope')
set(h1(1),'LineWidth',2,'DisplayName','$r(t)/C$');
set(h1(2),'LineWidth',2,'DisplayName','$a_A(t)/a_{\max}$');
set(h2(1),'DisplayName','birth rate','LineWidth',2);
set(h2(2),'LineWidth',2,'DisplayName','$c_d$');
clr=lines();
h2(1).Color=clr(4,:);
h2(2).Color=clr(5,:);
set(ax3(2),'YColor',clr(4,:),fn{:},lw{:});
set(ax3(1),fn{:},lw{:},'ylim',[0,1]);
legend(ax3(2),'Interpreter','latex','Location','eastoutside');
nexttile;ax4=gca;
contourf(ax4,yarr.mesh.time,yarr.mesh.tau,log10(abs(yarr.profile.ceff.id)));
colorbar(ax4);
title(ax4,'log10 population effect density');
set(ax4,fn{:},lw{:});
axis(ax4,'equal')
end