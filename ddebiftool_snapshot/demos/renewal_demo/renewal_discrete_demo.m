%% compare basic renewal equation example to its equivalent DDE with discrete delay
%
% original equation, shown in the demo renewal_demo:
% $$x(t)=\int_{-\tau_1-\tau_2}^{-\tau_2}
% \frac{\gamma}{2}x(t+s)*1-x(t+s))\mathrm{d}s$$
%
% defining $\tau_3=\tau_2+\tau_1$, the integral is from $-\tau_2$ to
% $\tau_3$ with length of integral $\tau_1$.
% (from [23] D. Breda, O. Diekmann, D. Liessi, and F. Scarabel. Numerical
% bifurcation analysis of a class of nonlinear renewal equations. Electron.
% J. Qual. Theory Differ. Equ., 65:1–24, 2016.)
%
% This equation can be reformulated as a DDAE with discrete delays and
% translational symmetry, shown for comparison is this demo:
%
% $$\dot u(t)= x(t)(1-x(t))-\beta$$
%
% $$0=\frac{\gamma}{2}[u(t-tau_2)-u(t-\tau_1-\tau_2)+\beta\tau_1]-x(t)$$
%
% where we impose $u=0$ or $u(0)=0$. The demo illustrate the treatment of
% translational symmetry in $u$ by introducing artificial parameter
% $\beta$, breaking/compensating this smmetry.
%% Load DDE-Biftool and extension into Path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_psol'],...
    [base,'ddebiftool_extra_rotsym/'],...
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities']);
%% Set number of delays and create parameter names as strings
% The demo has the parameters |gamma|, |tau1|, |tau2|, |tau3| and (for the DDE version) |beta|.
varnames={'x','u'};
cind=[varnames;num2cell(1:length(varnames))];
ix=struct(cind{:});
parnames={'gamma','tau1','tau2','tau3','beta'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
%% initial values for parameters and states
par0([ip.gamma, ip.tau1,ip.tau2,ip.tau3,ip.beta])=...
      [   1,        8,      1,      9,      0];
%% integral kernel
g=@(s,x,p)x.*(1-x);
dg={@(s,x,p,ds,dx,dp)dx-2*x.*dx,...
    @(s,x,p,ds,dx,dp)-2*dx.^2};
%% renewal equation & integral: 
%  0=gamma/2 * (u(t-tau2)-u(t-tau1-tau2))+beta*tau1 - x(t)
%  u'(t)=x(t)*(1-x(t))
f=@(x,u2,u3,gam,beta,tau1)...
    [(gam/2).*(u2-u3+beta.*tau1)-x;...
     g(0,x,[])-beta];
%% derivatives
df=@(x,u2,u3,gam,beta,tau1,  dx,du2,du3,dgam,dbeta,dtau1)...
    [(dgam/2).*(u2-u3+beta.*tau1)+gam/2.*(du2-du3+dbeta.*tau1+beta.*dtau1)-dx;...
     dg{1}(0,x,[], 0,dx,[])-dbeta];
d2f=@(x,u2,u3,gam,beta,tau1, dx,du2,du3,dgam,dbeta,dtau1)...
    [dgam.*(du2-du3)+(dgam/2).*(beta.*dtau1+dbeta.*tau1)+gam.*dbeta.*dtau1;...
     dg{2}(0,x,[], 0,dx,[])];
fc=@(xx,p)f(xx(1,1,:),xx(2,3,:),xx(2,4,:),p(1,ip.gamma,:),p(1,ip.beta,:),p(1,ip.tau1,:));
dfc={@(xx,p,dxx,dp)df(xx(1,1,:),xx(2,3,:),xx(2,4,:),p(1,ip.gamma,:),p(1,ip.beta,:),p(1,ip.tau1,:),...
    dxx(1,1,:),dxx(2,3,:),dxx(2,4,:),dp(1,ip.gamma,:),dp(1,ip.beta,:),dp(1,ip.tau1,:)),...
    @(xx,p,dxx,dp)d2f(xx(1,1,:),xx(2,3,:),xx(2,4,:),p(1,ip.gamma,:),p(1,ip.beta,:),p(1,ip.tau1,:),...
    dxx(1,1,:),dxx(2,3,:),dxx(2,4,:),dp(1,ip.gamma,:),dp(1,ip.beta,:),dp(1,ip.tau1,:))};
%% condition relating delays
taucond=dde_sys_cond_create('name','tau3fix','fun',...
    @(tau)tau(1,3,:)-tau(1,2,:)-tau(1,1,:),'args',{'parameter',{1,[ip.tau1,ip.tau2,ip.tau3]}});
u0ststcond=dde_sys_cond_create('name','u0fixstst','fun',...
    @(x)x,'args',{'x',{ix.u}});
u0psolcond=dde_sys_cond_create('name','u0fixpsol','fun',...
    @(x)x,'args',{'profile',{ix.u,1}});
%% basic function structure (without conditions)
derivs={'sys_dirderi',dfc};
funcs=set_funcs('sys_rhs',fc,derivs{:},'sys_tau',@()[ip.tau1,ip.tau2,ip.tau3],'lhs_matrix',diag([0,1]),...
    'x_vectorized',true,'p_vectorized',true);
%% initial [x;u], beta
u0(ix.x)=1-2/(par0(ip.gamma)*par0(ip.tau1));
u0(ix.u)=0;
u0=u0(:);
par0(ip.beta)=u0(ix.x)*(1-u0(ix.x));
%%
bd={'min_bound',[ip.gamma, 0;ip.tau1,0],'max_bound',[ip.gamma 5;ip.tau1,8]};
[ststfuncs,pos_eqs,suc] = SetupStst(funcs,'x',u0,'parameter',par0,...
    'contpar',[ip.gamma,ip.beta,ip.tau3],'step',0.1,'max_step',[ip.gamma, 0.1],...
    bd{:},...
    'usercond',[taucond,u0ststcond],'outputfuncs',true,...
    'print_residual_info',1,'extra_condition',true,...
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)p.x(ix.x)});
%%
figure(1);clf; ax1=gca;hold(ax1,'on');
[pos_eqs] = br_contn(ststfuncs,pos_eqs,40,'ax',ax1);
%% detect bifurcations (exclude 0 eigenvalue)
[pos_eqs_fine,~,stst_bifs,ind_bif]=MonitorChange(ststfuncs,pos_eqs,'distance',1e-3,'min_iterations',5,...
    'exclude_trivial',true,'locate_trivial',@(p)0) 
%% branch off at Hopf bifurcation
[psolfuncs,psolbr_ini,suc]=SetupPsol(funcs,pos_eqs_fine,ind_bif(1),...
    'excludefreqs',0,'SetupHopf.usercond',[taucond,u0ststcond],...
    'SetupHopf.extra_condition',true,...
    'usercond',[taucond,u0psolcond],'outputfuncs',true,...
    'print_residual_info',1,'extra_condition',true,...
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)max(p.profile(ix.x,:))});
%%
figure(1);
psolbr=br_contn(psolfuncs,psolbr_ini,100,'ax',ax1);
%% detect bifurcations (exclude 1 eigenvalue)
[psol_fine,~,psol_bifs,ind_psol]=MonitorChange(psolfuncs,psolbr,'distance',1e-3,'min_iterations',5,...
    'exclude_trivial',true,'locate_trivial',@(p)[1;1],'range',2:length(psolbr.point));
%% find dominant eigenvalue
[nunst_psolbif,dom_psolbif,err_psolbif]=GetStability(psol_bifs,...
    'exclude_trivial',true,'locate_trivial',@(p)[1;1]);
%% plot branches with stability
figure(1);clf;ax1=gca;
ltx={'Interpreter','LaTeX'};
Plot2dBranch(pos_eqs_fine,'funcs',ststfuncs,'locate_trivial',@(p)0,'y',@(p)p.x(ix.x),'ax',ax1);
hold(ax1,'on');
Plot2dBranch(psol_fine,'funcs',psolfuncs,'locate_trivial',@(p)[1;1],'y',@(p)max(p.profile(ix.x,:)),'ax',ax1);
%% Add bifurcations
% which are not automatically labelled for systems with continuous symmetry
hopfdeco={'ko','LineWidth',3,'DisplayName','Hopf bifurcation'};
lpdeco={'kd','LineWidth',1,'MarkerFaceColor','r','MarkerSize',12,'DisplayName','fold of limit cycles'};
pddeco={'ks','LineWidth',1,'MarkerFaceColor','b','MarkerSize',12,'DisplayName','period doubling'};
hopf_gam=arrayfun(@(p)p.parameter(ip.gamma),stst_bifs);
hopf_x=arrayfun(@(p)p.x(ix.x),stst_bifs);
po_gam=arrayfun(@(p)p.parameter(ip.gamma),psol_bifs);
po_x=arrayfun(@(p)max(p.profile(ix.x,:)),psol_bifs);
po_ev=round(dom_psolbif);
plot(ax1,hopf_gam,hopf_x,hopfdeco{:});
plot(ax1,po_gam(po_ev==-1),po_x(po_ev==-1),pddeco{:});
plot(ax1,po_gam(po_ev==1),po_x(po_ev==1),lpdeco{:});
xlabel(ax1,'$\gamma$',ltx{:});
ylabel(ax1,'$x$, $\max_tx$',ltx{:});
set(ax1,'FontSize',14)
%% Continue Hopf bifurcation in two parameters
[hopffuncs,hopf_ini,suc]=SetupHopf(funcs,pos_eqs_fine,ind_bif(1),'print_residual_info',1,'extra_condition',true,...
    'contpar',[ip.gamma,ip.tau1,ip.beta,ip.tau3],'dir',ip.tau1,'step',0.1,...
    'usercond',[taucond,u0ststcond],'outputfuncs',true,'excludefreqs',0,...
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)p.parameter(ip.tau1)});
%%
figure(2);clf;ax2=gca;
hopf=br_contn(hopffuncs,hopf_ini,100,'ax',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(hopffuncs,hopf,100,'ax',ax2);
%% Continue fold of periodic orbits
fprintf('Continuing fold of periodic orbits\n');
%% create extra condition for variational problem
% which replicates the condition of the original periodic problem: du(0)=0
xdim=length(u0);
du0psolcond=dde_sys_cond_create('name','du0fixpsol','fun',...
    @(x)x,'args',{'profile',{ix.u+xdim,1}});
figure(2);
[pofoldfuncs,pofold,suc]=SetupPOfold(funcs,psol_fine,ind_psol(po_ev==1),...
    'nullparind',ip.beta,...
    'contpar',[ip.gamma,ip.tau1,ip.beta,ip.tau3],'dir',ip.gamma,'step',0.1,...
    'usercond',[taucond,u0psolcond,du0psolcond],'outputfuncs',true,...
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)p.parameter(ip.tau1)});

pofold=br_contn(pofoldfuncs,pofold,30,'plotaxis',ax2);
%% Determine stability
[pofold_wstab,pofold_nunst]=br_stabl(pofoldfuncs,pofold,'exclude_trivial',true,'locate_trivial',@(p)[1;1;1]);
%% Continue period doubling
[pdfuncs,pdbranch,suc]=SetupPeriodDoubling(funcs,psol_fine,ind_psol(po_ev==-1),...
    'contpar',[ip.gamma,ip.tau1,ip.beta,ip.tau3],'dir',ip.gamma,'step',0.1,...
    'usercond',[taucond,u0psolcond],'outputfuncs',true,...
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)p.parameter(ip.tau1)});
figure(2);
pdbranch=br_contn(pdfuncs,pdbranch,40,'plotaxis',ax2);
pdbranch=br_rvers(pdbranch);
pdbranch=br_contn(pdfuncs,pdbranch,40,'plotaxis',ax2);
%% Determine Stability
[pdbranch_wstab,pd_nunst]=br_stabl(pdfuncs,pdbranch,'exclude_trivial',true,'locate_trivial',@(p)[1;1;-1]);
%% 
figure(3);clf;ax3=gca;
Plot2dBranch(hopf,'ax',ax3,'funcs',hopffuncs,'exclude_trivial',true,...
    'locate_trivial',@(p)[1;1i*p.omega;-1i*p.omega]);hold(ax3,'on');
Plot2dBranch(pofold_wstab,'ax',ax3,'funcs',pofoldfuncs,...
    'exclude_trivial',true,'locate_trivial',@(p)[1;1;1]);
Plot2dBranch(pdbranch_wstab,'ax',ax3,'funcs',pdfuncs,...
    'exclude_trivial',true,'locate_trivial',@(p)[1;1;-1]);
hold(ax3,'off');
xlabel('$\gamma$','Interpreter','latex');
ylabel('$\tau_1$','Interpreter','latex');
set(ax3,'FontSize',18);
