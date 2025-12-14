%% Minimalistic DDE-Biftool demo Mackey-Glass Equation 
%
% The Mackey-Glass equation is given by
% 
% $$x'(t)=\beta \frac{x(t-\tau)}{1+x(t-\tau)^n}-\gamma x(t)$$
% 
% Parameters are (in this order) |beta|, |n|, |tau| (|gamma| is not part of
% parameter vector).
% See manual for extensive comments.
%% load DDE-Biftool into path
clear
base=[pwd(),'../../']; % set to folder where DDE-Biftool routines are
addpath('../../ddebiftool',...       % load base routines
    '../../ddebiftool_extra_psol',...% for bifurcation of periodic orbits
    '../../ddebiftool_utilities');   % interface routines (almost always needed)
format compact
%% Define problem
[ib,in,itau,ig]=deal(1,2,3,4); % define parameter indices
f=@(x,xd,b,n,g)b.*xd./(1+xd.^n)-g.*x; % r.h.s., all vectorized for speed-up
%% Convert r.h.s. to DDE-Biftool f(x,p) format,
% note that the parameters are row vectors for historic reasons
sys_rhs=@(xx,p)f(xx(1,1,:),xx(1,2,:),p(1,ib,:),p(1,in,:),p(1,ig,:));
funcs=set_funcs('sys_rhs',sys_rhs,'sys_tau',@()itau,...
    'x_vectorized',true,'p_vectorized',true); % set up problem
%% Set initial guess for state and parameter
% set boundaries, max stepsize, and perform initial correction.
[b,n,tau,g]=deal(2, 10, 0,1);   % initial parameters
x0=(b/g-1)^(1/n);               % initial non-trivial equilibrium
bounds={'max_bound',[itau,2;ib,5],'max_step',[0,0.3]};
[eqbr,suc]=SetupStst(funcs,'x',x0,'parameter',[b,n,tau,g],...
    'step',0.1,'contpar',itau,bounds{:})
%% Compute, find stability and bifurcations of non-trivial equilibria 
disp('Nontrivial equilibria');
figure(1);clf;ax1=gca;
eqbr=br_contn(funcs,eqbr,10,'ax',ax1);
[eqbr,nunst_eqs]=br_stabl(funcs,eqbr);
ihopf=find(diff(nunst_eqs));
fprintf('Hopf bifurcation near point %d, tau=%g\n',ihopf,eqbr.point(ihopf).parameter(itau));
%% Continue Hopf bifurcation in two parameters b and tau
% starting from point where stability change was detected
% (eqbr.point(ihopf)).
[hopfbr,suc]=SetupHopf(funcs,eqbr,ihopf,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1)
figure(2);clf;ax2=gca; xlabel('b');ylabel('tau');
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
hopfbr=br_rvers(hopfbr);
hopfbr=br_contn(funcs,hopfbr,30,'ax',ax2);
[hopfbr,nunst_h]=br_stabl(funcs,hopfbr);
%% Re-plot Hopf bifurcation in 2 parameter plane
figure(3);clf;xlabel('b');ylabel('tau');
Plot2dBranch(hopfbr);
%% Branch off toward periodic orbits
[per_orb,suc]=SetupPsol(funcs,eqbr,ihopf,'intervals',20,'degree',4,'max_step',[itau,0.5])
%% Plot initial periodic orbit to illustrate discretization
figure(4);clf;
plot(per_orb.point(2).mesh,per_orb.point(2).profile,'o-');
xlabel('time/period');ylabel('x');
%% Continuation with "live plot" of periodic orbits
% Floquet multipliers are determined in call to br_stabl. The number of
% unstable Floquet multipliers, nunst_per, indicates that stability is lost
% for increasing tau.
figure(1);
per_orb=br_contn(funcs,per_orb,60,'ax',ax1);
[per_orb,nunst_per,dom_per]=br_stabl(funcs,per_orb);
%% Check that first stability change is a period doubling
ipd=find(diff(nunst_per)==1,1,'first')
per_orb.point(ipd+1).stability.mu(1:3)
dom_per(ipd+1)
%% plot all profiles, highlighting first period doubling
figure(4);clf;hold on;
for i=1:length(per_orb.point)
    pt=per_orb.point(i);
    plot(pt.mesh,pt.profile,'r-');
end
pt=per_orb.point(ipd+1);
plot(pt.mesh,pt.profile,'k-','linewidth',3);
xlabel('time/period');ylabel('x(t)');
%% Find period doubling bifurcations in two parameters
[pdfuncs,pdbr1,suc]=SetupPeriodDoubling(funcs,per_orb,ipd,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,bounds{:})
%% Continuation
figure(2);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
pdbr1=br_rvers(pdbr1);
pdbr1=br_contn(pdfuncs,pdbr1,30,'ax',ax2);
[pdbr1,nunst_pd]=br_stabl(pdfuncs,pdbr1);
%% Branch off at period doubling 
% (Solutions at far end get inaccurate.)
[per2,suc]=DoublePsol(funcs,per_orb,ipd);
figure(1);
per2=br_contn(funcs,per2,60,'ax',ax1);
[per2,nunst_p2,dom_p2,triv_defect2]=br_stabl(funcs,per2); 
%% Plot final one-parameter diagram
% Note that stability boundaries are not accurate as we have not refined
% the branch there. Use MonitorChange or LocateSpecialPoints for
% refinement.)
figure(5);clf;
Plot2dBranch({eqbr,per_orb,per2});
legend('Location','southeast')
xlabel('tau');ylabel('x');
%% Find secondary period doubling and track in two parameters
ipd2=find(diff(nunst_p2)==1,1,'first');
[pd2funcs,pdbr2,suc]=SetupPeriodDoubling(funcs,per2,ipd2,...
    'contpar',[ib,itau],'dir',ib,'step',1e-1,bounds{:});
figure(2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
pdbr2=br_rvers(pdbr2);
pdbr2=br_contn(pd2funcs,pdbr2,30,'ax',ax2);
[pdbr2,nunst_pd2]=br_stabl(pd2funcs,pdbr2);
%% Two-parameter bifurcation plot cleaned-up
% pass on pdfuncs for period doublings to help identify type of bifurcation
% (without this hint, pdbr1 and pdbr2 look just like periodic orbits)
figure(3);hold on
Plot2dBranch({pdbr1,pdbr2},'funcs',pd2funcs);
