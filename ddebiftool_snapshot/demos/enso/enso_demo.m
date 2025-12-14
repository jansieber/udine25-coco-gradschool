%% ENSO demo for implementation of periodic forcing and resonances
%
%%  ENSO model by Ghil, Zaliapin, Thompson
%
% M. Ghil, I. Zaliapin, and S. Thompson, "A delay differential model of
% ENSO variability: Parametric instability and the distribution of
% extremes," Nonlinear Processes Geophys. 15, 417–433 (2008). 
% 
% See also review by Keane, A., Krauskopf, B., & Postlethwaite, C. M.
% (2017). Climate models with delay differential equations. Chaos: An
% Interdisciplinary Journal of Nonlinear Science, 27(11).
%
% for bifurcation analysis. Equation as implemented here:
%
% $$ \begin{array}{rl}
% x'(t)&=-a \tanh(\kappa x(t-t\tau))+b y(t),\\
% y'(t)&=\lambda y(t)-\omega z(t)-y(t)(y(t)^2+z(t)^2)),\\
% z'(t)&=\omega y(t)+\lambda z(t)-z(t)(y(t)^2+z(t)^2)).
% \end{array}
% $$
%
% The equations for $y$ and $z$ implement the Hopf bifurcation normal form,
% resulting in solutions of the form $y(t)=\sqrt{\lambda}\cos(\omega t)$,
% if we force $z(0)=0$. Bifurcation analysis will be performed in tau and
% lambda, where $\sqrt{lambda}$ is the forcing amplitude. The forcing
% period is always $1$ year, such that $\omega=2\pi$ is fixed.
%% Load path
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool/'],...
    [base,'ddebiftool_extra_psol/'],...
    [base,'ddebiftool_extra_nmfm/'],...
    [base,'ddebiftool_utilities/'],...
    [pwd(),'/extra_resonance']);
format compact
%% Define variable names
parnames={'a','b','kappa','tau','omega','lambda'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
varnames={'x','y','z'};
cind=[varnames;num2cell(1:length(varnames))];
iv=struct(cind{:});
dim=length(varnames);
ntau=1;
%% Define right-hand side
rhs=@(x,xt,y,z,a,b,kappa,omega,lambda)[...
    -a.*tanh(kappa.*xt)+b.*y;...
    lambda.*y- omega.*z-y.*(y.^2+z.^2);...
     omega.*y+lambda.*z-z.*(y.^2+z.^2)];
sys_rhs=@(u,p)rhs(u(iv.x,1,:),u(iv.x,2,:),u(iv.y,1,:),u(iv.z,1,:),...
    p(1,ip.a,:),p(1,ip.b,:),p(1,ip.kappa,:),p(1,ip.omega,:),p(1,ip.lambda,:));
funcs=set_funcs(...
    'sys_rhs',sys_rhs,'sys_tau',@()ip.tau,'x_vectorized',true,'p_vectorized',true);
%% Initial parameter values and initial solution
par0([ip.a,ip.b,ip.kappa,ip.tau,ip.omega,ip.lambda])=...
     [   1,   1,     11,      0,    2*pi,     -1];
u0=zeros(dim,1);
bounds={'max_bound',[ip.tau,0.6;ip.lambda,1]};
getp=@(br,c)arrayfun(@(p)p.parameter(ip.(c)),br.point);
getT=@(br)arrayfun(@(p)p.period,br.point);
%% Track branch of trivial equilibria
% check stability change
[triv_eqs,suc]=SetupStst(funcs,'parameter',par0,'x',u0,...
    'contpar',ip.tau,'step',0.1,bounds{:});
figure(1);clf;ax1=gca;
triv_eqs=br_contn(funcs,triv_eqs,10,'ax',ax1);
[triv_eqs,nunst_eqs,eqbifs,eqbifind]=MonitorChange(funcs,triv_eqs);
%% branch off for autonomous periodic orbits (without forcing)
indhopf=find(nunst_eqs==2,1,'first');
[aut_per_tau,suc]=SetupPsol(funcs,triv_eqs,indhopf);
aut_per_tau=br_contn(funcs,aut_per_tau,50,'ax',ax1);
[aut_per_tau,nunst_pertau]=br_stabl(funcs,aut_per_tau);
figure(2);clf;
plot(getp(aut_per_tau,'tau'),getT(aut_per_tau),'.-');
grid on
xlabel('\tau');ylabel('period');
%% increase lambda tracking Hopf bifurcation
[hopftau,suc]=SetupHopf(funcs,triv_eqs,indhopf,'contpar',[ip.tau,ip.lambda],...
    'step',0.1,'dir',ip.lambda)
figure(3);clf;ax3=gca;
hopftau=br_contn(funcs,hopftau,60,'ax',ax3);
[hopftau,hopftests,hopfbifind,hopfbiftype]=LocateSpecialPoints(funcs,hopftau);
%% Extra conditrion to ensure that forcing has phase cos(t) for forced orbits
phase_cond=dde_sys_cond_create('name','forcedphase','fun',@(p)p,'args',{'profile',{iv.z,1}});
%% Start off from Hopf-Hopf point
t_args={'print_residual_info',1,...
    'stop_1_2',false,'degree',5,'intervals',31};
C1info=BranchFromCodimension2(funcs,hopftau,t_args{:},bounds{:},'step',1e-2);
%%
clear C1branches testfuncs;
istorusbr=find(strcmp({C1info.C1type},'TorusBifurcation'));
for i=length(C1info):-1:1
    if ismember(i,istorusbr)
        trlam=getp(C1info(i).branch,'lambda');
        %% enforce phase in forcing only if forcing nonzero
        if any(abs(trlam)>1e-7)
            psolphase=strcmp({C1info(i).funcs.sys_cond.name},'psol_phasecondition');
            C1info(i).funcs.sys_cond=C1info(i).funcs.sys_cond(~psolphase);
            C1info(i).funcs=dde_funcs_add_cond(C1info(i).funcs,phase_cond);
            C1info(i).branch.method.point.phase_condition=false;
            C1info(i).branch.parameter.max_step=[ip.tau,1e-2];
            C1info(i).branch=br_add_stop(C1info(i).branch,'name','rotation0',...
                'state','corrector','online',@(p)C1info(i).funcs.get_comp(p,'omega')>2);
        end
    end
    C1branches(i)=br_contn(C1info(i).funcs,C1info(i).branch,50,'ax',ax3,...
        'print_progress',1);
    [C1branches(i),testfuncs{i},indc2{i},ind2ctype{i}]=...
        MonitorChange(C1info(i).funcs,C1branches(i),'printlevel',2,'print_residual_info',0);
    [C1branches(i),C1_nunst{i},C1dom{i},C1defect{i}]=...
        br_stabl(C1info(i).funcs,C1branches(i),'exclude_trivial',true);
end
trbranches=C1branches(istorusbr);
trfuncs=[C1info(istorusbr).funcs];
hopfbranches=C1branches(strcmp({C1info.C1type},'hopf'));
hopflam=hopfbranches(any(getp(hopfbranches(1),'lambda')>1e-7));
is_forcedtr=arrayfun(@(br)any(getp(br,'lambda')>1e-7),trbranches);
forcedtr=trbranches(is_forcedtr);
unforcedtr=trbranches(~is_forcedtr);
%% Along torus bifurcations find resonances
[resonance,farey]=deal(cell(1,3));
max_denom=7;
for i=1:length(trbranches)
    [resonance{i},farey{i}]=dde_locate_resonances(trfuncs(i),trbranches(i),max_denom);
    tau=arrayfun(@(p)p.parameter(ip.tau),resonance{i});
    [dum,imaxtau]=max(tau);
    resonance{i}=resonance{i}(1:imaxtau);
    farey{i}=farey{i}(:,1:imaxtau);
    not12=farey{i}(2,:)>2;
    [resonance{i}(not12),suc]=dde_correct_resonance(trfuncs(i),trbranches(i),resonance{i}(not12),...
        'print_residual_info',1);
end
%%
figure(4);clf;ax4=gca;hold(ax4,'on');
Plot2dBranch({hopftau,C1branches(~istorusbr)},'ax',ax4);
Plot2dBranch(trbranches(1),'funcs',trfuncs(1),'ax',ax4);
Plot2dBranch(trbranches(2),'funcs',trfuncs(2),'ax',ax4);
for i=1:length(trbranches)
    for k=1:length(resonance{i})
        pt=resonance{i}(k);
        plot(ax4,pt.parameter(ip.tau),pt.parameter(ip.lambda),...
        'ko','markerfacecolor','k','DisplayName','');
        text(ax4,pt.parameter(ip.tau),pt.parameter(ip.lambda),...
        sprintf('%d/%d',farey{i}(1,k),farey{i}(2,k)),...
        'VerticalAlignment','baseline','FontSize',14);
    end
end
xlim(ax4,[0.1,0.25]);
ylim(ax4,[0,0.15]);
xlabel(ax4,'tau');
ylabel(ax4,'lambda (square of forcing amplitude)')
legend(ax4,findall(ax4.Children,'UserData','Plot2dBranch'))
%% increase lambda to branch off at lambda=0 for a stable periodic response
% initial run are equilibria to lambda=0, then we branch off at Hopf
[lam_eqs,suc]=ChangeBranchParameters(funcs,triv_eqs,1, ...
    'contpar',ip.lambda,'step',0.05,bounds{:})
figure(5);clf;ax5=gca;
lam_eqs=br_contn(funcs,lam_eqs,22,'ax',ax5);
[lam_eqs,nunst_eqlam]=br_stabl(funcs,lam_eqs);
indhopflam=find(nunst_eqlam==2,1,'first')
[psfuncs,psol_lam,suc]=SetupPsol(funcs,lam_eqs,indhopflam,...
    'print_residual_info',1,'intervals',30,'degree',5,'usercond',phase_cond,...
    'extra_condition',true,'phase_condition',false,...
    'SetupHopf.phase_condition',true,'outputfuncs',true)
psol_lam=br_contn(psfuncs,psol_lam,100,'ax',ax5);
[psol_lam,nunst_pslam]=br_stabl(psfuncs,psol_lam);
%% Scan periodic response for different forcing amplitudes
lamvals=linspace(0,(ceil(max(getp(forcedtr,'lambda'))*10)+1)/10,6);
lamvals=lamvals(2:end);
pslamvals=getp(psol_lam,'lambda');
figure(6);clf;ax6=gca;
xlabel(ax6,'tau');
clear nunst_psol
for i=1:length(lamvals)
    lam=lamvals(i);
    [dum,ind]=min(abs(pslamvals-lam));
    ps0=psol_lam.point(ind);
    ps0.parameter(ip.lambda)=lam;
    br0=setfield(psol_lam,'point',ps0);
    [psol(i),suc]=ChangeBranchParameters(psfuncs,br0,1,...
        'contpar',ip.tau,'max_step',[0,0.02]);
    psol(i)=br_contn(psfuncs,psol(i),100,'ax',ax6);
    psol(i)=MonitorChange(psfuncs,psol(i));
    [psol(i),nunst_psol{i}]=br_stabl(psfuncs,psol(i));
end
%%
figure(7);clf;ax7=gca;hold(ax7,'on');
xlabel(ax7,'tau');
for i=1:length(psol)
    Plot2dBranch(psol(i),'funcs',psfuncs,'ax',ax7,...
    'lgname',sprintf('lambda=%g',lamvals(i)));
end
%% Continue fold bifurcation of periodic responses
ipsolfold=3;
ipfold=find(diff(nunst_psol{ipsolfold}==1),1,'first');
[pfoldfuncs,pfold,suc]=SetupPOfold(psfuncs,psol(ipsolfold),ipfold,...
    'contpar',[ip.tau,ip.lambda],'step',0.1,'dir',ip.tau,'psol_phase_condition',false,...
    'max_step',[0,0.1;ip.lambda,0.005],'use_tangent',true,'minimal_angle',0);
pfold=br_add_stop(pfold,'name','forcing0',...
    'state','predictor','online',@(p)p(end).parameter(ip.lambda)<1e-3);
figure(3);
pfold=br_contn(pfoldfuncs,pfold,80,'ax',ax3);
pfold=br_rvers(pfold);
pfold=br_contn(pfoldfuncs,pfold,50,'ax',ax3);
%% Continue far torus bifurcation of periodic
ipsoltr=length(psol);
itr=find(diff(nunst_psol{ipsoltr}==2),1,'first');
[trfuncs(3),trbranches(3),suc]=SetupTorusBifurcation(psfuncs,psol(ipsoltr),itr,...
    'contpar',[ip.tau,ip.lambda],'step',0.1,'dir',ip.tau,'psol_phase_condition',false,...
    'max_step',[0,0.1;ip.lambda,0.05;ip.tau,0.02]);
trbranches(3)=br_add_stop(trbranches(3),'name','rotation0',...
    'state','corrector',...
    'online',@(p)trfuncs(3).get_comp(p(end),'omega')>2|trfuncs(3).get_comp(p(end),'omega')<0);
trbranches(3)=br_contn(trfuncs(3),trbranches(3),40,'ax',ax3);
trbranches(3)=br_rvers(trbranches(3));
trbranches(3)=br_contn(trfuncs(3),trbranches(3),40,'ax',ax3);
[trbranches(3),nunst_tr]=br_stabl(trfuncs(3),trbranches(3));
%% locate resonances on trps
[resonance{3},farey{3}]=dde_locate_resonances(trfuncs(3),trbranches(3),max_denom);
not12=farey{3}(2,:)>2;
[resonance{3}(not12),suc]=dde_correct_resonance(trfuncs(3),trbranches(3),resonance{3}(not12),...
        'print_residual_info',1);
%% Insert new curves into 2d bifurcation diagram
figure(4);
Plot2dBranch(pfold,'funcs',pfoldfuncs,'ax',ax4);
Plot2dBranch(trbranches(3),'funcs',trfuncs(3),'ax',ax4);
for k=1:length(resonance{3})
    pt=resonance{3}(k);
    plot(ax4,pt.parameter(ip.tau),pt.parameter(ip.lambda),...
        'ko','markerfacecolor','k','DisplayName','');
    text(ax4,pt.parameter(ip.tau),pt.parameter(ip.lambda),...
        sprintf('%d/%d',farey{3}(1,k),farey{3}(2,k)),...
        'VerticalAlignment','baseline','FontSize',14);
end
xlim(ax4,[0.1,0.45]);
ylim(ax4,[-1e-2,1]);
legend(ax4,findall(ax4.Children,'UserData','Plot2dBranch'),'location','northwest')
%%
save('enso_results.mat')