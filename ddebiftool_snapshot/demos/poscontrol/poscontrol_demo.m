%% Bifurcation analysis for for position control problem
%
% The problem originates from H.-O. Walther, "Stable periodic motion of a
% system with state- dependent delay," Differential and Integral Equations
% 15, 923– 944 (2002).
%
% A bat or submarine wants to control its position $x(t)$ to be at position
% $x_0$ using feedback control with gain $k$. Inside the feedback loop
% there is a fixed processing delay $\tau_0$. The "current position", as
% measured is called $x_m(t)$:
%
% $$ x'(t)=k(x_0-x_m(t-\tau_0)) $$
%
% The measurement is done by sending a sound wave to the wall at position
% $0$ and measuring how long it takes for the signal to arrive back. The
% quantity $s(t)$ is the traveling time of the signal to the wall and back.
% The factor $c$ is the speed of sound such that the position estimate
% $x_m(t)$ is
%
% $$ x_m(t)=\frac{c}{2}s(t) $$
%
% The sounds true traveling time of the signal arriving at time $t$ depends
% on the sum of the position when the signal was emitted, $x(t-s(t))$ and
% the current position, $x(t)$:
%
% $$ c s(t)=x(t-s(t))+x(t) $$
%
%% Load paths for symbolic r.h.s. generation
clear
format compact
base=[pwd(),'/../../'];
addpath([base,'ddebiftool']);
addpath([base,'ddebiftool_utilities']);
addpath([base,'ddebiftool_extra_nmfm']); % for normal forms of equilibria
addpath([base,'ddebiftool_extra_psol']); % for bifurcations of periodic orbits
%% Define parameter names and indices
% They have to be in the same order as when genreating the r.h.s.
parnames={'x0','tau0','k','c'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
%% Define r.h.s.
% At this point (when using |set_funcs| or the wrapper |set_symfuncs|) one
% passes on the optional argument |'lhs_matrix'|, which is the matrix $M$
% on the left-hand side.
funcs=set_symfuncs(@sym_poscontrol,'lhs_matrix',diag([1,0]));
%% Initial parameter values
parini([ip.tau0,ip.x0,ip.c,ip.k])=...
    [      0.1,    1,    2,  1.8];
xini=parini(ip.x0)*[1;1];
bounds={'min_bound',[ip.tau0,0;ip.k,0],'max_bound',[ip.tau0,2;ip.k,3]};
%% Follow equilibria, find Hopf bifurcation & its criticality
% From now on DDAEs are treated the same way as DDEs from the user
% perspective. Set up first two points along branch of equilibria. Continue
% in parameter |tau| .This will not change the equilibrium but its
% stability.
eqs=SetupStst(funcs,'x',xini,'parameter',parini,'step',0.1,...
    'contpar',ip.tau0,'max_step',[ip.tau0,0.05;ip.k,0.1],bounds{:});
% continue the branch (with online plotting)
figure(1);clf;ax1=gca;xlabel(ax1,'tau0');ylabel(ax1,'x (eq)');
eqs=br_contn(funcs,eqs,100,'plotaxis',ax1);
%% Check Stability along branch
% First output reports number of unstable eigenvalues, final output updates
% |branch.point|.
[eqs,eq_nunst]=br_stabl(funcs,eqs);
%% Show stability
% We note that the equilibrium loses its stability when |tau0| is inreased.
disp(eq_nunst.');
Plot2dBranch(eqs,'ax',ax1);
xlabel(ax1,'tau0');ylabel(ax1,'x (eq)');
%% Continue Hopf bifurcations in two parameter plane
% We continue in two parameters: |tau0| and |x0|, starting from the
% approximate Hopf point found previously.
indhopf=find(eq_nunst==2,1,'first');
[hopf,suc]=SetupHopf(funcs,eqs,indhopf,'contpar',[ip.tau0,ip.k],'dir',ip.k,'step',0.1);
figure(2);clf;ax2=gca;xlabel(ax2,'tau0');ylabel(ax2,'k');
hopf=br_contn(funcs,hopf,30,'ax',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(funcs,hopf,30,'ax',ax2);
[hopf_wbifs,hopffuncs,bif2ind,bif2type]=LocateSpecialPoints(funcs,hopf);
%% plot 2d bif
figure(3);clf;ax3=gca;
Plot2dBranch(hopf_wbifs,'ax',ax3);
%% Branch off at Hopf bifurcation for periodic orbits
% We find the Hopf bifurcation where we take the first unstable point as
% the initial guess. |SetupPsol| initializes first two points of branch near Hopf point
discpars={'degree',4,'intervals',20};
[per,suc]=SetupPsol(funcs,eqs,indhopf,'max_step',[0,0.5],...
    discpars{:},'print_residual_info',1,'matrix','sparse','eigmatrix','sparse',...
    'collocation_parameters','force_smooth');
figure(1);clf;ax1=gca;
per=br_contn(funcs,per,40,'plotaxis',ax1);
%% Check stability
% computation of eigenfunctions permits error estimate
per.method.stability.collocation_parameters='force_smooth';
per.method.stability.geteigenfuncs=true;
[per,pnunst,pdom,perr]=br_stabl(funcs,per);
%% Find time of minimal s for last periodic orbit
per_smin_t=arrayfun(@(p){dde_coll_roots(p,[0,1],'diff',1)'},per.point(2:end));
s_eval=@(p,t)[0,1]*dde_coll_eva(p,t(:)'); % evaluate x1 at t in point p
[dum,itmin]=min(s_eval(per.point(end),per_smin_t{end}));
tmin=per_smin_t{end}(itmin);
%% Continue grazing bifurcation, 
% when periodic motion nearly hits the wall: we add two conditions
% s(t0)=1e-2 and s'(t0)=0. We also add t0 as free parameter and continue in
% two system parameters.
ipe=ip;
ipe.t0=length(fieldnames(ip))+1;
ipe.val=ipe.t0+1;
graze0=per;
graze0.point=per.point(end);
graze0.point.parameter([ipe.t0,ipe.val])=[tmin,1e-2];
min_cond=@(p,pref)dde_extreme_cond(p,[0,1],ipe.val,ipe.t0);
[gfuncs,grazing,suc]=ChangeBranchParameters(funcs,graze0,1,...
    'contpar',[ipe.tau0,ipe.k,ipe.t0],...
    'usercond',{min_cond},'outputfuncs',true,...
    'print_residual_info',1,'plot_measure',[],'extra_condition',true);
figure(2);
grazing=br_contn(gfuncs,grazing,80,'ax',ax2);
grazing=br_rvers(grazing);
grazing=br_contn(gfuncs,grazing,30,'ax',ax2);

%% Derivative of time profiles of solution
% Check rate of change for travel time. This is the time-dependent part of
% the delay and should never have s'(t)>=1
figure(3);clf;ax3=gca;hold(ax3,'on');
for i=1:length(per.point)
    pt=per.point(i);
    pt.profile=dde_coll_eva(pt,pt.mesh,'diff',1);
    plot(ax3,pt.mesh,pt.profile(2,:)/pt.period,'.-','linewidth',2);
end
yline(-1);
yline(1);
xlabel(ax3,'time/period');
ylabel(ax3,'d/dt of travel time s of signal');
title(ax3,'Check that time-dependent delay changes with speed<1');
%% Continue fold of periodic orbits
% The family of periodic orbits showed a fold bifurcation. We continue it
% (note the slightly different output format for |SetupPOfold|). We start
% from the last unstable point. Below initializes the first two points on
% the branch.
ipfold=find(diff(pnunst)==-1,1,'first');
[pfuncs,pofold,suc]=SetupPOfold(funcs,per,ipfold,'contpar',[ip.tau0,ip.k],'dir',ip.k,...
    'max_step',[0,0.2],'step',1e-1,'matrix','sparse','print_progress',1,...
    'print_residual_info',1)
%% Continue in fold of periodic orbits in both directions
figure(2);
pofold=br_contn(pfuncs,pofold,40,'plotaxis',ax2,'print_residual_info',1);
pofold=br_rvers(pofold);
pofold=br_contn(pfuncs,pofold,20,'plotaxis',ax2);
%% Check Stability for PO folds (all stable except neutral direction)
% then replot 2d bifurcation diagram
[nunst_pf,dom_pf,triv_pf,pofold.point]=GetStability(pofold,'funcs',pfuncs,...
    'exclude_trivial',true);
disp(nunst_pf.');
figure(2);clf;ax2=gca;hold(ax2,'on');xlabel(ax2,'tau0');ylabel(ax2,'k');
Plot2dBranch(hopf_wbifs,'ax',ax2);
Plot2dBranch(pofold,'ax',ax2,'funcs',pfuncs);
Plot2dBranch(grazing,'ax',ax2,'funcs',gfuncs,'color',[0,0,1]);
set(ax2,'fontsize',18,'fontweight','bold','fontname','courier','linewidth',2,'box','on');
grid on;
