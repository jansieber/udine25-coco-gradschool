%% Test basic renewal equation example
%
% $$x(t)=\int_{\tau_2}^{\tau_1+\tau_2}
% \frac{\gamma}{2}x(t-s)\cdot(1-x(t-s))\mathrm{d}s$$
%
% from [23] D. Breda, O. Diekmann, D. Liessi, and F. Scarabel. Numerical
% bifurcation analysis of a class of nonlinear renewal equations. Electron.
% J. Qual. Theory Differ. Equ., 65:1–24, 2016.
%
% In DDE-Biftool, this equation is defined as
%
% $$0=x(t)-\frac{\gamma}{2}y(t-\tau_2)$$
%
% $$0=\int_0^{\tau_1}x(t-s)(1-x(t-s)) \mathrm{d}s-y(t)$$
%
% This equation can be reformulated as a DDAE with discrete delays and
% translational symmetry, which is shown in the demo
% |renewal_discrete_demo|:
%
% $$\dot u(t)= x(t)(1-x(t))-\beta$$
%
% $$0=\frac{\gamma}{2}[u(t-\tau_2)-u(t-\tau_2-\tau_1)+\beta\tau_1]-x(t)$$
%
% where we impose $u=0$ or $u(0)=0$.
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
%% renewal equation: x(t)=tau1*gamma/2 * y(t-tau2)
% enable derivatives by uncommenting for more accurate normal form
% computation and improved convergence along bifurcation curves
[ix,iy,ig,itau1,itau2]=deal(1,2,1,2,3);
f_funcs=set_funcs('sys_rhs',...
    @(x,p)x(ix,1,:)-p(1,ig,:).*x(iy,3,:)/2,...
    'sys_tau',@()[itau1,itau2],...
    'lhs_matrix',[0,0],...
    'x_vectorized',true,'p_vectorized',true);
%% initial values for parameters and states
par0([ig, itau1,itau2])=...
      [2,     2,    1];
u0(ix,1)=1-2/(par0(ig)*par0(itau1));
u0(iy,1)=u0(ix)*2/par0(ig);
tab_de=dde_add_funcs([],'rhs',f_funcs,'x',u0,'par',par0);
%% integral kernel
% enable derivatives by uncommenting for more accurate normal form
% computation and improved convergence along bifurcation curves
g=@(s,x,p)x.*(1-x);
tab_dist=dde_add_dist_delay(tab_de,g,...
    'value',iy,'bound',itau1,'x',ix,'ipar',[],'int',4,'degree',3);
funcs=dde_combined_funcs(tab_dist);
bd={'min_bound',[ig, 0;itau1,0],'max_bound',[ig 5;itau1,8]};
%% Exclude high frequencies
% for low-order approximations of distributed delays in renewal equations,
% stability may find spurious high-frequency instability. We exclude these
% eigenvalues from consideration
excl_highfreq={'exclude',@(p,lam)abs(imag(lam))>5};
%% Initial equilibria
[pos_eqs,suc] = SetupStst(funcs,'x',u0,'parameter',par0,...
    'contpar',ig,'step',0.1,'max_step',[ig, 0.1],bd{:});
%%
figure(1);clf; ax1=gca;hold(ax1,'on');
pos_eqs = br_contn(funcs,pos_eqs,40,'ax',ax1);
%% Stability
[pos_eqs,nunst_eqs]=br_stabl(funcs,pos_eqs,'recompute',true,excl_highfreq{:});
%% detect Hopf bifurcations
[pos_eqs_fine,~,ind_bif,stst_bifs]=LocateSpecialPoints(funcs,pos_eqs,...
    'print_residual_info',1,'distance',1e-3,'min_iterations',5,excl_highfreq{:});
%% branch off at Hopf bifurcation
[psolbr_ini,suc]=SetupPsol(funcs,pos_eqs_fine,ind_bif,'print_residual_info',1,...
    'plot_measure',{@(p)p.parameter(ig),@(p)max(p.profile(1,:))});
%%
psolbr=br_contn(funcs,psolbr_ini,100);
%% Periodic orbits may also have spurious eigenvalues
% We exclude them by finding their eigenfunction and checking if its
% discretization error is large, indicating that this is a high-frequency
% problem 
excl_err={'exclude',@(p,mu)p.stability.err>0.05,'geteigenfuncs',true};
%%
[psolbr,nunst_psol,dom_psol,err_psol]=br_stabl(funcs,psolbr,'recompute',true,excl_err{:});
%% compute error
gam=arrayfun(@(p)p.parameter(ig),psolbr.point);
exactsol=@(sigma,A,t,T)sigma+A*sin(pi/2*t*T);
sigma=1/2+0.25*pi./gam;
A=sqrt(2*sigma.*(1-1./gam-sigma));
pts=psolbr.point;
errval=NaN(1,length(pts));
for i=1:length(pts)
    pt=pts(i);
    [~,errval(i)]=fminbnd(@(t)norm(exactsol(sigma(i),A(i),pt.mesh+t,pt.period)-pt.profile(ix,:)),0,1);
end
%% continue Hopf bifurcation in two parameters
[hopf_ini,suc]=SetupHopf(funcs,pos_eqs_fine,ind_bif,'print_residual_info',1,'extra_condition',false,...
    'contpar',[ig,itau1],'dir',itau1,'step',0.1);
figure(2);clf;ax2=gca;
hopf=br_contn(funcs,hopf_ini,100,'ax',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,100,'ax',ax2);
%%
[hopf,nunst_hopf]=br_stabl(funcs,hopf,'recompute',true,excl_highfreq{:});
%%
[hopf_fine,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hopf,...
    'distance',1e-3,'min_iterations',5,excl_highfreq{:});
%% Branch off at generalized Hopf bifurcation
t_args={'print_residual_info',1};
C1info=BranchFromCodimension2(funcs,hopf_fine,bd{:},t_args{:});
%% Continue fold of periodic orbits
fprintf('Continuing %s\n',[C1info(1).C1type]);
figure(2);
pofold=br_contn(C1info(1).funcs,C1info(1).branch,30,'plotaxis',ax2);
%%
[pofold_wstab,pofold_nunst]=br_stabl(C1info(1).funcs,pofold,...
    'exclude_trivial',true,excl_err{:},'recompute',true);
%%
ind_pd=find(diff(nunst_psol),1,'first');
[pdfuncs,pdbranch,suc]=SetupPeriodDoubling(funcs,psolbr,ind_pd,...
    'contpar',[ig,itau1],'dir',itau1,'step',0.1,'plot_measure',[])
%%
figure(2);
pdbranch=br_contn(pdfuncs,pdbranch,40,'plotaxis',ax2);
pdbranch=br_rvers(pdbranch);
pdbranch=br_contn(pdfuncs,pdbranch,40,'plotaxis',ax2);
%%
[pdbranch_wstab,pd_nunst]=br_stabl(pdfuncs,pdbranch,'exclude_trivial',true,...
    excl_err{:});
%% 
figure(3);clf;ax3=gca;
Plot2dBranch(hopf_fine,'ax',ax3,excl_highfreq{:});hold(ax3,'on');
Plot2dBranch(pofold_wstab,'ax',ax3,'funcs',C1info(1).funcs,excl_err{:});
Plot2dBranch(pdbranch_wstab,'ax',ax3,'funcs',pdfuncs,excl_err{:});
hold(ax3,'off');
xlabel('$\gamma$','Interpreter','latex');
ylabel('$\tau_1$','Interpreter','latex');
set(ax3,'FontSize',18);
