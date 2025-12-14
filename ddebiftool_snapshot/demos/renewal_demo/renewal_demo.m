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
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities']);
%% Set number of delays and create parameter names as strings
% The demo has the parameters |gamma|, |tau1|, |tau2|, |tau3| and (for the DDE version) |beta|.
varnames={'x','y'};
cind=[varnames;num2cell(1:length(varnames))];
ix=struct(cind{:});
parnames={'gamma','tau1','tau2','tau3','beta'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
%% initial values for parameters and states
par0([ip.gamma, ip.tau1,ip.tau2,ip.tau3,ip.beta])=...
      [   2,        2,      1,      3,      0];
u0(ix.x)=1-2/(par0(ip.gamma)*par0(ip.tau1));
u0(ix.y)=u0(ix.x)*2/par0(ip.gamma);
u0=u0(:);
%% integral kernel
% enable derivatives by uncommenting for more accurate normal form
% computation and improved convergence along bifurcation curves
g=@(s,x,p)x.*(1-x);
dg={... 
    @(s,x,p,ds,dx,dp)dx-2*x.*dx,...
    @(s,x,p,ds,dx,dp)-2*dx.^2 ...
    };
%% renewal equation: x(t)=tau1*gamma/2 * y(t-tau2)
% enable derivatives by uncommenting for more accurate normal form
% computation and improved convergence along bifurcation curves
f=@(xx,p)xx(1,1,:)-p(1,ip.gamma,:).*xx(2,3,:)/2;
df={... 
    @(xx,p,dxx,dp)dxx(1,1,:)-1/2*(...
    p(1,ip.gamma,:).*dxx(2,3,:)+...
    dp(1,ip.gamma,:).* xx(2,3,:)),...
    @(xx,p,dxx,dp)-dp(1,ip.gamma,:).*dxx(2,3,:) ...
    };
f_funcs=set_funcs('sys_rhs',f,'sys_dirderi',df,'sys_tau',@()[ip.tau1,ip.tau2,ip.tau3],...
    'lhs_matrix',[0,0],'x_vectorized',true,'p_vectorized',true);
tab_de=dde_add_funcs([],'rhs',f_funcs,'x',u0,'par',par0);
tab_dist=dde_add_dist_delay(tab_de,[{g},dg],...
    'value',ix.y,'bound',ip.tau1,'x',ix.x,'ipar',[],'int',4,'degree',3);
funcs=dde_combined_funcs(tab_dist);
contpar = ip.gamma;
bd={'min_bound',[ip.gamma, 0;ip.tau1,0],'max_bound',[ip.gamma 5;ip.tau1,8]};
%% Exclude high frequencies
% for low-order approximations of istributed delays in renewal equations,
% stability may find spurious high-frequency instability. We exclude these
% eigenvalues from consideration
excl_highfreq={'exclude',@(p,lam)abs(imag(lam))>5};
%% Initial equilibria
[pos_eqs,suc] = SetupStst(funcs,'x',u0,'parameter',par0,...
    'contpar',ip.gamma,'step',0.1,'max_step',[contpar, 0.1],...
    bd{:},...
    'print_residual_info',1);
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
    'plot_measure',{@(p)p.parameter(ip.gamma),@(p)max(p.profile(1,:))});
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
gam=arrayfun(@(p)p.parameter(ip.gamma),psolbr.point);
exactsol=@(sigma,A,t,T)sigma+A*sin(pi/2*t*T);
sigma=1/2+0.25*pi./gam;
A=sqrt(2*sigma.*(1-1./gam-sigma));
pts=psolbr.point;
errval=NaN(1,length(pts));
for i=1:length(pts)
    pt=pts(i);
    [~,errval(i)]=fminbnd(@(t)norm(exactsol(sigma(i),A(i),pt.mesh+t,pt.period)-pt.profile(ix.x,:)),0,1);
end
%% continue Hopf bifurcation in two parameters (ignoring tau3 here)
[hopf_ini,suc]=SetupHopf(funcs,pos_eqs_fine,ind_bif,'print_residual_info',1,'extra_condition',false,...
    'contpar',[ip.gamma,ip.tau1],'dir',ip.tau1,'step',0.1);
figure(2);clf;ax2=gca;
hopf=br_contn(funcs,hopf_ini,100,'ax',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,100,'ax',ax2);
%%
[hopf,nunst_hopf]=br_stabl(funcs,hopf,'recompute',true,excl_highfreq{:});
%%
[hopf_fine,hopftests,hc2_indices,hc2_types]=LocateSpecialPoints(funcs,hopf,...
    'distance',1e-3,'min_iterations',10,excl_highfreq{:});
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
    'contpar',[ip.gamma,ip.tau1],'dir',ip.tau1,'step',0.1,'plot_measure',[])
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
