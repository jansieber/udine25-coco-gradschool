%% Test Daphnia demo
% see |daphnia_math.tex| for mathematical formulas
%% Load DDE-Biftool and extension into Path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities']);
s=load('index_daphnia.mat');
[     ip,  ix,  ip_sd,  ix_sd,  ip_ceff,  ix_ceff]=deal(...
    s.ip,s.ix,s.ip_sd,s.ix_sd,s.ip_ceff,s.ix_ceff);
%% initial values for parameters and states
par0=NaN(1,length(fieldnames(ip)));
par0([ip.szB,ip.szA,ip.szmax,ip.growth,ip.sigma,ip.feed])=...
      [  0.8,   2.5,    6,      0.15,      7,      1.8];
par0([ip.rep,ip.mu,ip.cap,ip.flow,ip.amax,ip.amaxrel])=...
      [  0.1,  0.25,  0.5,   0.5,    70,      1];
%% Equations for equilibria
% ststfun is function of r and parameters that equals 1 in equilibrium
% ststcomp returns other components of equilibrium for given r and parameters
ststfun=dde_sym2fun(@sym_ststeq,'rhs');
ststcomp=dde_sym2fun(@sym_stst,'rhs');
nx=length(fieldnames(ix));
u0=NaN(nx,1);
u0(ix.r)=fzero(@(x)ststfun(x,par0)-1,[0,2]);
u0(2:end)=ststcomp(u0(ix.r),par0);
%% define rhs DDAE part
lhsmat=zeros(sym_daphnia_dde('nf'),nx);
lhsmat(1,1)=1;
f_funcs=set_symfuncs(@sym_daphnia_dde,...
    'lhs_matrix',lhsmat);
tab_de=dde_add_funcs([],'rhs',f_funcs,'x',u0,'par',par0);
%% Define chain of integrals
ix_chain=struct('r',1,'bd',2);
taugrid=[0,0.05,0.1,0.15,0.2:0.1:1]; % fine mesh
% taugrid=[0,0.1,0.2,0.5,1]; % coarse mesh
idd=dde_create_chain_delay(ip.amax,[ix.r;ix.bd],...
    'grid',taugrid,'degree',4);
idd=dde_add_chain_delay(idd,'sd',@sym_sd_int,...
    'igx',ix_chain.r,'igp',ip_sd);
idd=dde_add_chain_delay(idd,'ceff',@sym_ceff_int,...
    'igx',{...
           'x' ,[ix_chain.r,ix_chain.bd],'id';...
           'sd',  1,                     'int' ...
           },...
    'igp',ip_ceff);
tab_dist=dde_close_chain_delay(tab_de,idd,...
    'ix',{'sd',1,'int';'ceff',1,'int'},...
    'value',[ix.bd;ix.cd;ix.sdA],...
    'M',{[0,-1; 0,0; 1,0], [0,1; 0,1; 0,0]}, ...%[0,0,-1,1; 0,0,0,1; 1,0,0,0],...
    'it',[ix.raA,ip.amaxrel],'t_type',{'x','parameter'});
funcs=dde_combined_funcs(tab_dist);
getpar=@(br,ind)arrayfun(@(p)p.parameter(ind),br.point);
bounds={'min_bound',[ip.mu, 0.05;ip.cap,0.05],...
    'max_bound',[ip.mu 1;ip.cap,2]};
%%
[pos_eqs,suc] = SetupStst(funcs,'x',u0,'parameter',par0,...
    'contpar',ip.mu,'step',0.01,'max_step',[ip.mu, 0.01;0,0.1],...
    bounds{:},...
    'print_residual_info',1);
%%
pos_eqs=br_add_stop(pos_eqs,'name','minimalsize','state','corrector',...
    'online',@(p)p(end).x(ix.bd)<0);
figure(1);clf; ax1=gca;hold(ax1,'on');set(ax1,'FontSize',20);
xlabel(ax1,'mortality $\mu$','Interpreter','latex');
ylabel(ax1,'resource $r$ max','Interpreter','latex');
pos_eqs=br_contn(funcs,pos_eqs,40,'ax',ax1);
pos_eqs=br_rvers(pos_eqs);
pos_eqs=br_contn(funcs,pos_eqs,40,'ax',ax1);
%% compare to analytical results
mustst=getpar(pos_eqs,ip.mu);
xstst=cell2mat(arrayfun(@(p){p.x},pos_eqs.point));
for i=length(mustst):-1:1
    xtest(:,i)=get_daphnia_stst(par0,ix,ip,'mu',mustst(i),xstst(ix.r,i));
end
%% check stability: several Hopf bifurcations are suggested
[pos_eqs,nunst_eqs,dom_eqs,err_eqs]=br_stabl(funcs,pos_eqs,0,1);
%%
[pos_eqs_fine,~,ind_bif,stst_bifs]=LocateSpecialPoints(funcs,pos_eqs,...
    'print_residual_info',1,'distance',1e-3,'min_iterations',5,'debug',true);
%% Continue transcritical bifurcation
ptrans0=pos_eqs_fine.point(ind_bif(1));
e_bd=zeros(nx,1);
e_bd(ix.bd)=1;
b0_cond=@(p,pref)dde_stst_lincond(p,1,'x','stateproj',e_bd','condprojmat',1,'trafo',0);
[trfuncs,transcr_ini,suc]=SetupFold(funcs,pos_eqs_fine,ind_bif(1),'print_residual_info',1,...
    'contpar',[ip.mu,ip.cap],'dir',ip.cap,'step',0.01,'max_step',[0,0.5],...
    'extracolumns',1,'usercond',{b0_cond},'outputfuncs',true,'extra_condition',true);
figure(2);clf;ax2=gca;set(ax2,'FontSize',20);
xlabel(ax2,'mortality $\mu$','Interpreter','latex');
ylabel(ax2,'carrying capacity of resource $C$','Interpreter','latex');
transcr_ini.method.continuation.stops=[];
transcr=br_contn(trfuncs,transcr_ini,100,'ax',ax2);
transcr=br_rvers(transcr);
transcr=br_contn(trfuncs,transcr,100,'ax',ax2);
%%
[transcr,nunst_tr]=br_stabl(trfuncs,transcr);
%% continue Hopf bifurcation in two parameters (ignoring tau3 here)
[hopf_ini,suc]=SetupHopf(funcs,pos_eqs_fine,ind_bif(2),'print_residual_info',1,...
    'contpar',[ip.mu,ip.cap],'dir',ip.cap,'step',0.01,'max_step',[0,0.05]);
figure(2);
hopf=br_contn(funcs,hopf_ini,100,'ax',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,100,'ax',ax2);
%%
[hopf_fine,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hopf,...
    'exclude',@(p,z)abs(imag(z))>0.6);
%% branch off at Hopf bifurcation
[psolbr_ini,suc]=SetupPsol(funcs,pos_eqs_fine,ind_bif(2),'print_residual_info',1,...
    'plot_measure',{@(p)p.parameter(ip.mu),@(p)max(p.profile(ix.r,:))},...
    'max_step',[ip.mu, 0.01;0,0.5],'degree',6,'intervals',20,'radius',1e-3);%,'use_tangent',1,'minimal_angle',0.5);
%% stop when maximum of raA becomes large
e_raA=zeros(nx,1);
e_raA(ix.raA)=1;
maxsel=@(i,v)max(abs(v(i,:)));
raAslope=@(p)maxsel(ix.raA,dde_coll_eva(p,p.mesh,'diff',1));
slopemax=10;
psolbr_ini=br_add_stop(psolbr_ini,'name','maxraAslope','state','corrector',...
    'online',@(p)raAslope(p(end))>slopemax);
%%
figure(1);
psolbr=br_contn(funcs,psolbr_ini,50,'ax',ax1);
%%
figure(3);clf;ax3=gca;hold(ax3,'on');
set(ax3,'FontSize',20);
xlabel(ax3,'time $t$','Interpreter','latex');
ylabel(ax3,'maturation age $a_A(t)$','Interpreter','latex');
title(ax3,'maturation age along branch of periodic orbits');
tfine=linspace(0,1,10000);
for i=1:length(psolbr.point)
    pt=psolbr.point(i);
    y=dde_coll_eva(pt,tfine);
    plot(tfine*pt.period,y(ix.raA,:)*pt.parameter(ip.amax),'-',...
        pt.mesh*pt.period,pt.profile(ix.raA,:)*pt.parameter(ip.amax),'k.');
    title(sprintf('i=%d of %d, mu=%g',i,length(psolbr.point),pt.parameter(ip.mu)));
    pause(0.1);
    drawnow;
end
%%
[psolbr,nunst_psol,dom_psol,err_psol]=br_stabl(funcs,psolbr,0,1);
%% Plot bifurcation diagram with birth rate as y-axis
figure(1);clf;ax1=gca;hold(ax1,'on');
bval=@(x)par0(ip.sigma).*x(ix.r,:)./(1+par0(ip.sigma).*x(ix.r,:)).*x(ix.bd,:); %unscaled birth rate

Plot2dBranch(pos_eqs_fine,'funcs',funcs,'ax',ax1,'y',@(p)bval(p.x));
Plot2dBranch(psolbr,'funcs',funcs,'ax',ax1,'y',@(p)max(bval(p.profile),[],2));
xlabel(ax1,'mortality $\mu$','Interpreter','latex');
ylabel(ax1,'(max) birth rate $b$','Interpreter','latex');
legend(ax1,'Location','EastOutSide');
ylim(ax1,[0,4e-3])
set(ax1,'FontSize',20);
%% From branch of psol detect point with maximal relative maturation age~0.6
t_extreme=arrayfun(@(p){dde_coll_roots(p,e_raA','diff',1)},psolbr.point(2:end));
max_r=cellfun(@(p,t)max(e_raA'*dde_coll_eva(p,t)),num2cell(psolbr.point(2:end)),t_extreme);
[~,ind_r]=min(abs(max_r-0.6));
p0=p_remesh(psolbr.point(ind_r));
max_t=dde_coll_roots(p0,e_raA','diff',1);
[max_r0,indt]=max(e_raA'*dde_coll_eva(p0,max_t));
t0=max_t(indt);
contour_daphnia(4,tab_dist,ip,ix,p0);
%% There is a peak of rel. a_A rising to above 0.6
% we aim to track solutions with peak of this height in two parameters
ipe=ip;
ipe.t0=length(fieldnames(ip))+1;
ipe.val=ipe.t0+1;
maxr_cond=@(p,pref)dde_extreme_cond(p,e_raA',ipe.val,ipe.t0);
p0.parameter([ipe.t0,ipe.val])=[t0,max_r0];
%% track periodic orbits with fixed maximal negative slope as approximation for singularity
% we switch off phase condition and instead fix t0 such that the inflection
% point always lies inside a collocation interval. We also switch off mesh
% adaptation to avoid that the inflection point crosses a collocation
% interval boundary.
mbranch=setfield(psolbr,'point',p0);
[mfuncs,mbranch,suc]=ChangeBranchParameters(funcs,mbranch,1,...
    'contpar',[ipe.mu,ipe.cap],'extra_condition',true,...
    'phase_condition',false,...
    'usercond',{maxr_cond},'outputfuncs',true,...
    'print_residual_info',1,'plot_measure',[],...
    bounds{:},...
    'adapt_mesh_after_correct',false,'step',0.1,...
    'newton_max_iterations',6);
mbranch.method.continuation.stops=[];
%%
figure(2);ax2=gca;
mbranch=br_contn(mfuncs,mbranch,100,'plotaxis',ax2);
mbranch=br_rvers(mbranch);
mbranch=br_contn(mfuncs,mbranch,100,'plotaxis',ax2);
%%
[mbranch,nunst_m,dom_m,triv_m]=br_stabl(mfuncs,mbranch,'recompute',1,...
    'max_number_of_eigenvalues',4);
%%
figure(3);clf;ax3=gca;hold(ax3,'on');
set(ax3,'FontSize',20);
xlabel(ax3,'time $t$','Interpreter','latex');
ylabel(ax3,'maturation age $a_A(t)$','Interpreter','latex');
title(ax3,sprintf('maturation age along branch of periodic orbits\nwhere its maximum is fixed large'));
for i=1:length(mbranch.point)
    pt=mbranch.point(i);
    y=dde_coll_eva(pt,tfine);
    plot(ax3,pt.mesh*pt.period,pt.profile(ix.raA,:)*pt.parameter(ip.amax),'k.',...
        tfine*pt.period,y(ix.raA,:)*pt.parameter(ip.amax),'-');
    title(ax3,sprintf('i=%d of %d, mu=%g',i,length(mbranch.point),pt.parameter(ip.mu)));
    pause(0.1);
    drawnow;
end
%% Track fold of periodic orbits (works for fine mesh)
t_args={'print_residual_info',1,'adapt_mesh_after_correct',false,'step',1e-2,'plot_measure',[]};
ipf=find(nunst_psol==1,1,'last');
[pfuncs,pofold,suc]=SetupPOfold(funcs,psolbr,ipf,'contpar',[ip.mu,ip.cap],'dir',[ip.mu],'step',0.1,t_args{:});
pofold=br_contn(pfuncs,pofold,30,'ax',ax2);
pofold=br_rvers(pofold);
pofold=br_contn(pfuncs,pofold,30,'ax',ax2);
[pofold,pfnunst]=br_stabl(pfuncs,pofold);
%%
figure(2);clf;ax2=gca;
Plot2dBranch(hopf_fine,'lgname','first Hopf bifurcation');
hold(ax2,'on')
Plot2dBranch(transcr,'lgname','Tanscritical bifurcation');
Plot2dBranch(pofold,'lgname','Fold p.o.','funcs',pfuncs);
Plot2dBranch(mbranch,'parameter',[ip.mu,ip.cap],'lgname',sprintf('Singularity\n(max of rel a_A=%4.2f)\n',max_r0));
mutr=getpar(transcr,ip.mu);
xlabel(ax2,'mortality $\mu$','Interpreter','latex');
ylabel(ax2,'carrying capacity of resource $C$','Interpreter','latex');
set(ax2,'FontSize',20);
legend(ax2,'Location','eastoutside');

%%
save('results_daphnia_fine.mat')
