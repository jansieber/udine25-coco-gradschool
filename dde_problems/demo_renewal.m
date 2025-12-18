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
% The demo has the parameters |gamma|, |tau1|, |tau2|, |tau3| and (for the DDE version) |beta|.
%% Load DDE-Biftool and extension into Path
clear
ddebiftool_path(fullfile('..','ddebiftool_snapshot'));
%% Renewal equation: x(t)=tau1*gamma/2 * y(t-tau2)
% enable derivatives by uncommenting for more accurate normal form
% computation and improved convergence along bifurcation curves
[ix,iy,ig,itau1,itau2]=deal(1,2,1,2,3);
f=  @(x,y,g)           x(1,:)- g.* y(3,:)/2;
df={@(x,y,g,dx,dy,dg) dx(1,:)-(g.*dy(3,:)+dg.* y(3,:))/2,...
    @(x,y,g,dx,dy,dg)        -dg.*dy(3,:)};
f_funcs=set_funcs('sys_rhs',{[ix,iy],ig,f},'sys_dirderi',df,'sys_tau',@()[itau1,itau2],...
    'lhs_matrix',[0,0],'x_vectorized',true,'p_vectorized',true);
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
dg={... 
    @(s,x,p,ds,dx,dp)dx-2*x.*dx,...
    @(s,x,p,ds,dx,dp)-2*dx.^2 ...
    };
tab_dist=dde_add_dist_delay(tab_de,[{g},dg],...
    'value',iy,'bound',itau1,'x',ix,'ipar',[],'int',4,'degree',3);
funcs=dde_combined_funcs(tab_dist);
bd={'min_bound',[ig, 0;itau1,0],'max_bound',[ig 5;itau1,7.5]};
%% Exclude high frequencies
% for low-order approximations of distributed delays in renewal equations,
% stability may find spurious high-frequency instability. We exclude these
% eigenvalues from consideration
excl_highfreq={'exclude',@(p,lam)abs(imag(lam))>5};
%% Initial equilibria
[pos_eqs,suc] = SetupStst(funcs,'x',u0,'parameter',par0,...
    'contpar',ig,'step',0.1,'max_step',[ig, 0.1],bd{:});
txt={'fontsize',16};
ltx={'interpreter','latex','fontsize',20};
lw={'linewidth',2};
clr=lines();
figure(1);clf; ax1=gca;hold(ax1,'on');set(ax1,txt{:},lw{:});
xlabel(ax1,'$\gamma$',ltx{:});
ylabel(ax1,'$x$',ltx{:});
%%
pos_eqs = br_contn(funcs,pos_eqs,40,'ax',ax1);
%% Observe high-frequency instability
[pos_eqs_whigh,nunst_eqs_whigh]=br_stabl(funcs,pos_eqs,'recompute',true);
indhigh=find(nunst_eqs_whigh>0,1,'first');
pt=pos_eqs_whigh.point(indhigh);
st=pt.stability;
lam=st.l0;
figure(4);clf; ax4=gca;hold(ax4,'on');set(ax4,txt{:},lw{:});
plot(ax4,real(lam),imag(lam),'o',lw{:},'Color',clr(1,:));
grid(ax4,'on');
xline(ax4,0,lw{:});
yline(ax4,0,lw{:});
xlabel(ax4,'Re',ltx{:});
ylabel(ax4,'Im',ltx{:});
title(ax4,sprintf('Eigenvalues of discretised RE at  $\\gamma=%g$',pt.parameter(ig)),ltx{:});
%% detect Hopf bifurcations
[pos_eqs_fine,~,ind_bif,stst_bifs]=LocateSpecialPoints(funcs,pos_eqs,...
    'distance',1e-3,'min_iterations',5,excl_highfreq{:});
%% branch off at Hopf bifurcation
[psolbr_ini,suc]=SetupPsol(funcs,pos_eqs_fine,ind_bif,...
    'plot_measure',{@(p)p.parameter(ig),@(p)max(p.profile(1,:))});
%%
figure(1);
psolbr=br_contn(funcs,psolbr_ini,100,'ax',ax1);
%% Periodic orbits may also have spurious eigenvalues
% We exclude them by finding their eigenfunction and checking if its
% discretization error is large, indicating that this is a high-frequency
% problem 
excl_err={'exclude',@(p,mu)p.stability.err>0.05,'geteigenfuncs',true};
%% Observe spurious Floquet multipliers with eigenfunctions oscillating at discretization level
[psolbr_whigh,nunst_psol_whigh,~,err_psol_whigh]=br_stabl(funcs,psolbr,'recompute',true,'geteigenfuncs',true);
pt=psolbr_whigh.point(end);
ef=pt.stability.eigenfuncs(3);
[mu,err]=deal(pt.stability.mu(3),pt.stability.err(3));
figure(5);clf; ax5=gca;hold(ax5,'on');set(ax5,txt{:},lw{:});
plot(ax5,ef.mesh,[real(ef.profile);imag(ef.profile)],'o-',lw{:});
grid(ax5,'on');
xlabel(ax5,'$t/T$',ltx{:});
ylabel(ax5,'$\mbox{Re}\,z(t)$, $\mbox{Im}\,z(t)$',ltx{:});
title(sprintf('Eigenfunction at $\\gamma=%g$ for $\\mu=%4.2f+\\mathrm{i}\\cdot%4.2f$, error$=%4.2f$',...
    pt.parameter(ig),real(mu),imag(mu),err),ltx{:});

%% Excluding the spurious Floquet multipliers
[psolbr,nunst_psol,dom_psol,err_psol]=br_stabl(funcs,psolbr,'recompute',true,excl_err{:});
%% compute error (analytical solution is known)
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
[hopf_ini,suc]=SetupHopf(funcs,pos_eqs_fine,ind_bif,'extra_condition',false,...
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
C1info=BranchFromCodimension2(funcs,hopf_fine,'step',0.01,bd{:});
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
%% save results
save('renewal_results.mat');
