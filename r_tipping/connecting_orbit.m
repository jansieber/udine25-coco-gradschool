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
%% r.h.s. and its derivatives
frhs={@(y,par)r_tipping_rhs(iv,ip,y,par)};
%% each subsystem (u_gamma, u_plus and u_minus) has  a full set of parameters
% their names will be prepended by ug, up and um
name_prep=@(prep,names)cellfun(@(s)[prep,'.',s],names,'UniformOutput',false);
%% Setup BVP for u_-
% b.c.:
bc_um=@(d,T,x0,x1,p)[...
    (x0(1)+p(ip.cm))^2-p(ip.a);...
    x1(2)-cos(p(ip.phi))/sqrt(2);...
    x1(3)-sin(p(ip.phi))/sqrt(2)];
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 300],'norm',inf);
% initial parameters
um_pars([ip.r, ip.cm, ip.beta, ip.omega, ip.cp, ip.a, ip.phi])=...
    [      1,    0,     0.5,       2*pi,     4,     1,    0];
um_pars=um_pars(:);
% initialize (x,p,q) such that the satisfy all boundary conditions
[                   a0,            phi0,          cm0]=...
    deal(um_pars(ip.a),um_pars(ip.phi),um_pars(ip.cm));
um_x=[-sqrt(a0)-cm0;cos(phi0)/sqrt(2);sin(phi0)/sqrt(2)];
% start from a short time segment (T_minus<<1)
Tm0=1e-3;
t0ini=linspace(0,Tm0,20)';
umxini=repmat(um_x,1,length(t0ini))';
umnames=name_prep('um',pnames);
%%
prob=ode_isol2bvp(prob,'u_minus',frhs{:},t0ini,umxini,um_pars,umnames,bc_um);
% find index of period and add as parameter
[fdata,uidx]=coco_get_func_data(prob,'u_minus.bvp.seg1.coll','data','uidx');
prob=coco_add_pars(prob,'Tu_minus',uidx(fdata.coll_seg.maps.T_idx),'Tu_minus');
bd_incTm=coco(prob,'incTm',[],1,'Tu_minus',[0,4]); 
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
%% setup BVP for u_gamma
% we integrate backwards in time for 20 periods in the future limit system,
% then take one period of the result
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 20);
ug_pars=um_pars;
omega0=ug_pars(ip.omega);
f=@(t,x)get1(feval(frhs{1},[x;cos(omega0*t);sin(omega0*t)],ug_pars));
[~,xg0]=ode45(f,[0,-20*pi/omega0],5);
[tg0,xg0]=ode45(f,[0,-2*pi/omega0],xg0(end));
treverse=tg0(end:-1:1)-tg0(end);
xreverse=xg0(end:-1:1);
ug0=[xreverse,cos(omega0*tg0),sin(omega0*tg0)];
bc_gam=@(~,T,x0,x1,p)[x0-x1;x1(2)*sin(p(ip.phi))-x1(3)*cos(p(ip.phi))];
ugnames=name_prep('ug',pnames);
prob=ode_isol2bvp(prob,'u_gamma',frhs{:},treverse,ug0,ug_pars,ugnames,bc_gam);
[fdata,uidx]=coco_get_func_data(prob,'u_gamma.bvp.seg1.coll','data','uidx');
prob=coco_add_pars(prob,'Tu_gamma',uidx(fdata.coll_seg.maps.T_idx),'Tgamma');
coco(prob,'fixUg',[],0,'Tgamma',[0,4]);

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
%% setup BVP for u_+
% First re-read u_gamma. This re-reading happens quite often, so there is a
% separate function below. We again start from a short segment
% (Tu_plus<<1). The phase of the end point varies and the phase of Gamma^u
% varies with it as it is glued to the end point of u_plus in condition
% pglue.
clear fdata
bdg=coco_bd_read('fixUg');
epglabs=coco_bd_labs(bdg,'EP');
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [0 600],'norm',inf);
[prob,uidx,u0,maps]=reread_sols(prob,{'u_gamma'},{'fixUg'},{epglabs(end)},ip);
% add u_+ problem
% only two b.c. are purely on the u_+ segment
bc_up=@(d,T,x0,x1,p)[...
   x0(2)-cos(p(ip.phi))/sqrt(2);...
   x0(3)-sin(p(ip.phi))/sqrt(2)];
up_pars=um_pars;
Tp0=1e-3;
t0ini=linspace(0,Tp0,20)';
up_x=[u0.u_gamma(uidx.u_gamma(maps.u_gamma.x1_idx(1)));cos(phi0)/sqrt(2);sin(phi0)/sqrt(2)];
upxini=repmat(up_x,1,length(t0ini))';
upnames=name_prep('up',pnames);
prob=ode_isol2bvp(prob,'u_plus',frhs{:},t0ini,upxini,up_pars,upnames,bc_up);
% find index of period and add as parameter
[fdata.u_plus,uidx.u_plus]=coco_get_func_data(prob,'u_plus.bvp.seg1.coll','data','uidx');
maps.u_plus=fdata.u_plus.coll_seg.maps;
prob=coco_add_pars(prob,'Tu_plus',uidx.u_plus(maps.u_plus.T_idx),'Tu_plus');
% glue u_g(1) and u_+(1) together (x value and phase of (p,q)), collected in
% extra function below
prob=match_plus_gamma(prob,'pglue',uidx,maps); 
bd_incTp=coco(prob,'incTp',[],1,{'Tu_plus','ug.phi','Tu_gamma'},[0,4]);
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
%% combine u_-, u_+ and u_gamma
bd_incTm=coco_bd_read('incTm');
bd_incTp=coco_bd_read('incTp');
incTplabs=coco_bd_labs(bd_incTp,'EP');
incTmlabs=coco_bd_labs(bd_incTm,'EP');
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', 100,'norm',inf);
seglist={'u_gamma','u_plus','u_minus'};
runlist={'incTp','incTp','incTm'};
lablist={incTplabs,incTplabs,incTmlabs};
glue_segments=[id_pars('u_minus','u_plus',[]),id_pars('u_gamma','u_plus',ip.phi)];
[prob,uidx,u0,maps]=reread_sols(prob,seglist,runlist,lablist,ip,...
    'match_plus_gamma','pglue','add_gap_monitor','eta','gap_parname','eta',...
    'identify_parameters',glue_segments); %#ok<ASGLU>
nonphi=setdiff(1:npars,ip.phi);
freepars=[{'eta','um.r','Tu_gamma','ug.phi'},ugnames(nonphi),upnames];
bd_closegap=coco(prob,'closegap',[],1,freepars,[-0.2,0.2]);
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
%% continue Lin gap in phi
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', 100,'norm',inf);
[prob,uidx,u0,maps]=reread_sols(prob,seglist,'closegap',gap0lab,ip,...
    'match_plus_gamma','pglue','add_gap_monitor','gap','gap_parname','eta','identify_parameters',glue_segments); %#ok<ASGLU>);
nonphi=setdiff(1:npars,ip.phi);
prob=coco_add_pars(prob,'W_sigma',...
    [uidx.u_minus(maps.u_minus.x1_idx(1));uidx.u_plus(maps.u_plus.x0_idx(1))],{'WuPs','WsGu'});
freepars=[{'um.phi','eta','WuPs','WsGu','Tu_gamma','ug.phi'},ugnames(nonphi),upnames];
bd_eta_of_phi=coco(prob,'eta_of_phi',[],1,freepars,[-pi,pi]);
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
%xlabel(ax3,'$\phi/\pi$','Interpreter','latex');
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
%% continue critical rate in phi for fixed a
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', 100,'norm',inf,'h0',5e-3);
seglist={'u_gamma','u_plus','u_minus'};
[prob,uidx,u0,maps]=reread_sols(prob,seglist,'closegap',gap0lab,ip,...
    'match_plus_gamma','pglue','add_gap_monitor','gap','gap_parname','eta',...
    'identify_parameters',glue_segments); %#ok<ASGLU>);
%
nonphi=setdiff(1:npars,ip.phi);
freepars=[{'um.r','um.phi','Tu_gamma','ug.phi'},ugnames(nonphi),upnames];
h_dev=1e-3;
prob=coco_add_event(prob,'UZ','um.phi',h_dev*0.5*[-1,1]);
bd_r_of_phi=coco(prob,'r_of_phi',[],1,freepars,{[],[-pi,pi]});
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
%% continue critical rate in phi for fixed a adding copy
UZ=coco_bd_labs('r_of_phi','UZ');
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', 100,'norm',inf);
[prob,uidx,u0,maps]=reread_sols_copy(prob,'r_of_phi',UZ(1:2),ip,'h_dev',h_dev,'init',true);
prob=coco_add_event(prob,'TAN','dr',0);
%
freepars=[{'dr','um.phi','um.r','Tu_gamma','ug.phi','copy.Tu_gamma'},ugnames(nonphi),upnames];
coco(prob,'r_of_phi_wcopy',[],1,freepars,{[],[-pi,pi]});
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
%% continue tangencies in r,a,phi gradually increasing T_+, T_-
bd_r_of_phi_wcopy=coco_bd_table('r_of_phi_wcopy');
TAN=coco_bd_labs('r_of_phi_wcopy','TAN');
for k=1:2
    prob=coco_prob();
    prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false,'NTSTMX',300);
    prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [200,300],'norm',inf,'NPR',10,'MaxRes',100);
    prob=reread_sols_copy(prob,'r_of_phi_wcopy',TAN(k),ip,'h_dev',h_dev,'fix_r_T','rtfix','init',false);
    freepars=[{'um.r','um.a','um.phi','Tu_gamma','Tu_plus','Tu_minus','copy.Tu_plus','copy.Tu_minus','ug.phi','copy.Tu_gamma'},ugnames(nonphi),upnames];
    runid=sprintf('a_r_tangency_%d',k);
    coco(prob,runid,[],1,freepars,[1e-1,3]);
end
%% continue critical rate for phi=0 in a-r plane
prob=coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 30,'MXCL',false);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'PtMX', [4000,200],'norm',inf,'NPR',100);
prob=reread_sols(prob,seglist,'closegap',gap0lab,ip,...
    'match_plus_gamma','pglue','add_gap_monitor','gap','gap_parname','eta',...
    'identify_parameters',glue_segments,'fix_r_T','rtfix');
freepars=[{'um.r','um.a','Tu_gamma','Tu_plus','Tu_minus','ug.phi'},ugnames(nonphi),upnames];
coco(prob,'a_r_phi=0',[],1,freepars,{[1e-1,3]});
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
%% r.h.s. for R-tipping problem
function f=r_tipping_rhs(iv,ip,y,par)
% assign names to parameters and variables
[       r,          cm,          beta,          omega,          cp,          a]=deal(...
 par(ip.r,:),par(ip.cm,:),par(ip.beta,:),par(ip.omega,:),par(ip.cp,:),par(ip.a,:)); % parameter names
[        x,        p,        q]=deal(...
    y(iv.x,:),y(iv.p,:),y(iv.q,:));
%% r.h.s
rho=p.^2+q.^2;
srho=sqrt(rho);
L=cm+rho.*(cp-cm)+beta.*q.*srho;
nlin=r-r.*rho;
f=[(x+L).^2-a;   p.*nlin+omega.*q;  -omega.*p+q.*nlin];
end
