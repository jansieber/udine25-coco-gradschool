%% Dynamics of a inverted pendulum balancing on a cart with delayed PD feedback
%
% This demo illustrates continuation of symmetry-breaking bifurcations for
% equilibria and peridoic orbits.
%
% This file uses the right-hand side generated using symbolic toolbox in
% |gen_sym_balancing|.
%
% Differential equations are 
% 
% $$x''(t)=\sin x(t)-\cos x(t)[ax(t-\tau)+bx'(t-\tau)]$$
%
% after non-dimensionalization. ant 0<=phi<=1, 0.5<=psi<=1.
%
% Reference:
% 
% J. Sieber, B. Krauskopf, Bifurcation analysis of an inverted pendulum
% with delayed feedback control near a triple-zero eigenvalue. Nonlinearity
% 17 (1), pp. 85-104, 2004.
%
%%
clear
addpath([pwd(),'/../../ddebiftool']);
addpath([pwd(),'/../../ddebiftool_utilities'],...
    [pwd(),'/../../ddebiftool_extra_nmfm'],...
    [pwd(),'/../../ddebiftool_extra_psol'])    
parnames={'a','b','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
funcs=set_symfuncs(@sym_balancing,'sys_tau',@()ip.tau);
%% initial parameters
parini([ip.a,ip.b,ip.tau])=...
    [    0.8, 1.6,  1];
xini=[0;0];
xdim=length(xini);
bounds={'max_bound',[ip.a,3;ip.b,4.5],...
    'min_bound',[ip.a,0.7;ip.b,0.8]};
%% Follow equilibria
% We know that the equilibria haev a reflection symmetry with matrix R. So,
% we append the condition Rsym*x-x=0.  The options |'extracolumns','auto'|
% then adds artificial parameters (that will be approximately zero), which
% ensure that the Jacobian is square. Conditions of the form
% |Pc*(R-I*exp(2*pi*p/q))*Ps*z=0| or for equilibria can be created wih
% |dde_stst_lincond|, where |z| may be one of the fields |x| or |v| (in the
% bifurcation types,fold or hopf). The arguments |Pc|, |R|, |Ps|, |p,q| are
% passed on as named arguments |'condprojmat'| (default |eye(xdim)|),
% |'trafo'|(default |eye(xdim)|), |'stateproj'| (default |eye(xdim)|), and
% |'rotation'| (default |[0,1]|). The extra conditions get passed on by
% named argument |'usercond'|, requiring |'outputfuncs'| set to |true|,
% such that the first output argument is the modified set of functions
% (including the new conditions in |sys_cond|). With the extra conditions
% branch points (pitchfork bifurcations) will not occur as singular points
% along the branch. The name-value argument pairs are listed in
% |dde_lincond_struct.m|
Rsym=-eye(xdim);
pini=dde_stst_create('x',xini,'parameter',parini);
ststcond=@(p)dde_stst_lincond(p,'x','trafo',Rsym);
[sfuncs,eqs,suc]=SetupStst(funcs,'point',pini,'step',0.1,...
    'contpar',ip.a,'max_step',[0,0.1],bounds{:},...
    'usercond',ststcond,'outputfuncs',true,'extra_condition',true,...
    'extracolumns','auto') %#ok<*NOPTS>
figure(1);clf;ax1=gca;xlabel(ax1,'a');
eqs=br_contn(sfuncs,eqs,100,'plotaxis',ax1);
%% Find Hopf bifurcation & its criticality and pitchfork
% The detection of the symmetry breaking points, if provided with the
% extended functions |sfuncs| will use the standard determining system for
% folds, but with extra conditions on |x|, which will make the determining
% system regular for the branch point. Normal form computations are not
% safe for systems with symmetry (the fold normal form coefficient is zero
% for pitchforks).
[eqs_wbifs,eqs_tests,ind_bifeqs,bifeqs_types]=LocateSpecialPoints(sfuncs,eqs);
%% Track symmetry breaking bifurcations for stst
ibp=ind_bifeqs(strcmp('fold',bifeqs_types));
bp_v_sym=@(p,pref)dde_stst_lincond(p,'v','trafo',zeros(size(Rsym)),...
    'rotation',[0,1],'condprojmat',eye(size(Rsym)));
bp_scale=@(p)dde_scale_cond(p,p,bp_v_sym,'res',1);
bp_scale_ini=@(p)dde_scale_cond(p,p,bp_v_sym,'res',0);
bpcond={bp_scale,ststcond};
%
[bpfuncs,bpbr,suc]=SetupFold(sfuncs,eqs_wbifs,ibp(1),'contpar',[ip.a,ip.b],...
    bounds{:},'max_step',[0,0.2],'outputfuncs',true,'extracolumns','auto',...
    'norm',false,'recompute',true,...
    'dir',ip.b,'step',0.02,'usercond',bpcond,'v_scal',bp_scale_ini);
% continue
figure(2);clf;ax2=gca;xlabel(ax2,'a');ylabel(ax2,'b');
bpbr=br_contn(bpfuncs,bpbr,100,'plotaxis',ax2);
bpbr=br_rvers(bpbr);
bpbr=br_contn(bpfuncs,bpbr,100,'plotaxis',ax2);
% Stability
[bp_wbifs,bp_nunst,bpbifs,bpbifind]=MonitorChange(bpfuncs,bpbr,...
    'printlevel',1,'print_residual_info',0);
%% Continue Hopf bifurcations
% Symmetry of Hopf eigenvectors is |Rsym*v-exp(1i*pi)*v=0| for our
% reflection symmetry, which enforces a symmetry condition with identically
% zero coefficient matrix in this case. We till add this redundant pair of
% conditions as |hopfcond| to demonstrate how to add more symmetry
% conditions to a bifurcation tracking problem. For the current problem the
% standard phase and norm condition of the Hopf bifurcation system are
% suitable. However, they can be replaced by setting |'phase_condition'|
% and |'norm'| to 0 in the |point.method| structure and adding other
% conditions instead. The function |dde_scale_cond| takes an arbitrary
% |sys_cond| |c(p)| and converts it into the condition |c(p)^T*c(p)-res=0|.
% In this way, one may enforce the breaking of one symmetry in the
% eigenvector.
ihopf=ind_bifeqs(strcmp(bifeqs_types,'hopf'));
hopfcond=@(p)dde_stst_lincond(p,'v','trafo',Rsym,'rotation',[1,2]);
hopf_v_id=@(p,pref)dde_stst_lincond(p,'v','trafo',zeros(size(Rsym)),...
    'rotation',[0,1],'condprojmat',eye(size(Rsym)));
hopf_scale=@(p,pref)dde_scale_cond(p,pref,hopf_v_id,'res',1);
hopf_phase=@(p,pref)dde_scale_cond(p,pref,hopf_v_id,'res',0,...
    'matrix',kron([0,1;-1,0],eye(xdim)),'ref',true);
hvcond={hopf_scale,hopf_phase,ststcond,hopfcond};
hopf_scale_ini=@(p)dde_scale_cond(p,[],hopf_v_id,'res',0);
hopfargs={'norm',false,'phase_condition',false,'recompute',true,...
    'initcond',hopfcond,'usercond',hvcond,'v_scal',hopf_scale_ini};
%
[shfuncs,hopf,suc]=SetupHopf(funcs,eqs_wbifs,ihopf,'contpar',[ip.a,ip.b],'dir',ip.b,...
    bounds{:},'max_step',[0,2e-1],'outputfuncs',true,'extracolumns','auto',...
    'stops',{@(p)p(end).omega<0},hopfargs{:},'print_residual_info',1)
%% Continue Hopf bifurcations
figure(2);
hopf=br_contn(shfuncs,hopf,100,'plotaxis',ax2);
hopf=br_rvers(hopf);
hopf=br_contn(shfuncs,hopf,60,'plotaxis',ax2);
[hopf,nunst_hopf]=br_stabl(funcs,hopf);
%% Find special points for symmetric Hopf bifurcations
% Note that the special points found are not generic BT of zero-Hopf
% points. Their normal forms are the same as the 1:2 resonance
% interpolating flow cases for the 'BT', and Hopf-Hopf interaction for the
% 'zero-Hopf' (where one dimension is not averaged but genuinely
% one-dimensional). Thus, one expects a non-symmetric Hopf bifurcation, a
% symmetry-breaking bifurcation of periodic orbits and a torus bifurcation
% to emerge from the 'zero-Hopf' point.
[hopf_wbifs,hopffuncs,bif2ind,bif2type]=LocateSpecialPoints(...
    shfuncs,hopf,'printlevel',1,'debug',true);

%% plot 2d bif
figure(3);clf;ax3=gca;hold(ax3,'on');xlabel(ax3,'a');ylabel(ax3,'b');
Plot2dBranch(bpbr,'ax',ax3,'funcs',bpfuncs);
Plot2dBranch(hopf_wbifs,'ax',ax3);
legend(ax3,'location','eastoutside');
%% branch off symmetric periodic orbit
% enforce symmetry of periodic orbit at times 0:0.1:0.5 (6 times, resulting
% in 6x2 symmetry conditions). As part of SetupPsol, the Hopf point is
% recomputed. We can pass on arguments to SetupHopf by prefixing them with
% |SetupHopf|.
figure(1);clf;ax1=gca;xlabel(ax1,'a');ylabel(ax1,'period');
psolsym=@(p)dde_psol_lincond(p,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)'*[1,1]);
addprefix=@(p,args)reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
    'UniformOutput',false),args(2:2:end)),1,[]);
phopfargs=addprefix('SetupHopf',hopfargs);
[pfuncs,sympo,suc]=SetupPsol(funcs,eqs_wbifs,ihopf,...
    'degree',6,'intervals',50,'print_residual_info',1,...
    'point.matrix','sparse','point.extra_condition',1,...
    'point.extracolumns','auto',...
    'max_step',[0,1;ip.a,5e-3],...
    'continuation.plot_measure',{@(p)p.parameter(ip.a),@(p)p.period},...
    'continuation.stops',{@(p)p(end).period>50},'usercond',psolsym,'outputfuncs',true,...
    phopfargs{:})
sympo=br_contn(pfuncs,sympo,200,'plotaxis',ax1);
%% Determine special points along branch of symmetric periodic orbits
% Determining special points along branches of periodic orbits is
% restricted to bisecting for stability changes. Symmetry breaking and fold
% bifurcations are both found as having an exrta Floquet multiplier equal
% to 1.
[sympo_wbifs,sympo_nunst,sympobifs,sympobifind]=MonitorChange(pfuncs,sympo,...
    'range',2:length(sympo.point),'printlevel',1,'print_residual_info',0,...
    'min_iterations',5);
%% plot bifurcation diagram for sym. periodic orbits phi vs period
figure(4);clf;ax4=gca;
Plot2dBranch(sympo_wbifs,'y',@(p)p.period,'max_nunst',2,'ax',ax4,'funcs',pfuncs,...
        'color',[0.5,0,0.75],'lgname','symmetric psol');
hold(ax4,'on');
a_symev1=arrayfun(@(p)p.parameter(ip.a),sympobifs);
T_symev1=arrayfun(@(p)p.period,sympobifs);
plot(ax4,a_symev1,T_symev1,'ko','MarkerFaceColor','k','displayname','POEV1');
set(legend(ax4),'Location','best');
xlabel(ax4,'a');
ylabel(ax4,'period');
%% continue symmetry breaking
isb=sympobifind(1);
sbxsym=@(p)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)');
sbvsym=@(p)dde_psol_lincond(p,xdim,'v','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)');
sbvscale=@(p,pref)dde_scale_cond(p,pref,@(p,pref)sbvsym(p),'res',1);
sbvscale_ini=@(p,pref)dde_scale_cond(p,pref,@(p,pref)sbvsym(p),'res',0);
poev1args={'usercond',{sbxsym,sbvscale},'POEV1_norm',false,'v_scal',sbvscale_ini};
[sbfuncs,sbbranch,suc]=SetupPOEV1(funcs,sympo_wbifs,isb(1),...
    'contpar',[ip.a,ip.b],'dir',ip.b,'step',-0.1,'max_step',[],...
    'minimal_angle',0,'plot_measure',[],'print_residual_info',1,'use_tangent',true,...
    poev1args{:})
sbbranch.method.continuation.stops={@(p)max(max(p(end).profile,[],1)-min(p(end).profile,[],1))<1e-4};
figure(2);
sbbranch=br_contn(sbfuncs,sbbranch,200,'plotaxis',ax2);
sbbranch=br_rvers(sbbranch);
sbbranch=br_contn(sbfuncs,sbbranch,20,'plotaxis',ax2);
%% Transversal stability of points
[sbbranch,nunst_sb,dom_sb,triv_sb]=br_stabl(sbfuncs,sbbranch);
[sb_wbifs,sb_tests,sb_bifs,sb_bifind]=MonitorChange(sbfuncs,sbbranch,...
    'printlevel',2,'print_residual_info',0,'min_iterations',5,'output','branch');
sb_wbifs.point(sb_bifind)=arrayfun(@(p)setfield(p,'flag','POEV1'),sb_wbifs.point(sb_bifind));
%% Insert symmetry breaking of periodic orbits into 2d bifurcation diagram
figure(3);
Plot2dBranch(sb_wbifs,'ax',ax3,'funcs',sbfuncs,'lgname','PO pitchfork');
%% continue POfold
condproj=eye(xdim);
pfxsym=@(p)dde_psol_lincond(p,xdim,'x','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)');
pfvsym=@(p)dde_psol_lincond(p,xdim,'profile','trafo',Rsym,'shift',[1,2],...
    'condprojint',linspace(0,0.5,6)');
pofoldargs={'usercond',{pfxsym,pfvsym},'initcond',{pfxsym,pfvsym}};
[pffuncs,pfbranch,suc]=SetupPOfold(funcs,sympo_wbifs,sympobifind(3),...
    'contpar',[ip.a,ip.b],'dir',ip.b,'step',-0.1,'max_step',[],...
    'minimal_angle',0,'plot_measure',[],'print_residual_info',1,'use_tangent',true,...
    pofoldargs{:});
pfbranch.method.continuation.stops={@(p)p(end).period>50};
%%
figure(2);
pfbranch=br_contn(pffuncs,pfbranch,200,'plotaxis',ax2);
pfbranch=br_rvers(pfbranch);
pfbranch=br_contn(pffuncs,pfbranch,200,'plotaxis',ax2);
%% Transversal stability of points
[pfbranch,nunst_pf,dom_pf,triv_pf]=br_stabl(pffuncs,pfbranch);
[pf_wbifs,pftests,pf_bifs,pf_bifind]=MonitorChange(pffuncs,pfbranch,...
    'printlevel',2,'print_residual_info',0,'min_iterations',5);
[dum,dom]=GetStability(pf_wbifs,'exclude_trivial',true);
i_pfev1=find(real(dom')>0&abs([diff(pftests),0])==1);
pf_wbifs.point(i_pfev1)=arrayfun(@(p)setfield(p,'flag','POEV1'),pf_wbifs.point(i_pfev1));
%% include symmetric folds of periodic orbits into 2d bifurcation diagram
figure(3);
Plot2dBranch(pf_wbifs,'ax',ax3,'funcs',pffuncs,'lgname','PO fold','color',[0,0.5,0]);
%% Branch off toward non-symmetric orbit at first symmetry breaking PO
nspoev1args=addprefix('SetupPOEV1',poev1args);
[nsfuncs,nonsymper,suc]=SetupPsol(funcs,sympo_wbifs,sympobifind(1),'branch_off','POEV1',...
    'print_residual_info',1,'outputfuncs',true,nspoev1args{:},...
    'continuation.stops',{@(p)p(end).period>50},...
    'contpar',ip.a,...
    'continuation.plot_measure',{@(p)p.parameter(ip.a),@(p)p.period});
figure(1);
nonsymper=br_contn(nsfuncs,nonsymper,100,'plotaxis',ax1);
%% Check Stability for POs
[nonsymper,nunst_nsp,dom_nsp,triv_nsp]=br_stabl(nsfuncs,nonsymper,'recompute',1,...
    'max_number_of_eigenvalues',4);
[ns_wbifs,ns_tests,ns_bifs,ns_bifind]=MonitorChange(nsfuncs,nonsymper,'print_residual_info',0,...
    'range',2:length(nonsymper.point))
[dum,dom,triv]=GetStability(ns_wbifs,'exclude_trivial',true);
i_nspd=find(real(dom')<0&abs([diff(ns_tests),0])==1);
ns_wbifs.point(i_nspd)=arrayfun(@(p)setfield(p,'flag','PD'),ns_wbifs.point(i_nspd));
i_nsev1=find(real(dom')>00&abs([diff(ns_tests),0])==1);
ns_wbifs.point(i_nsev1)=arrayfun(@(p)setfield(p,'flag','POEV1'),ns_wbifs.point(i_nsev1));
%% Insert nonsymmetric orbit into 1d bifurcation diagram
figure(4)
Plot2dBranch(ns_wbifs,'y',@(p)p.period,'max_nunst',2,'ax',ax4,'funcs',nsfuncs,...
    'color',[0,0.5,0],'lgname','non-symmetric psol');
%% try to determine period doubling
pd_cond=@(p)dde_pd_cond(p,xdim,linspace(0,0.75,4));
[pdfuncs,pdbr,suc]=SetupPeriodDoubling(funcs,ns_wbifs,ns_bifind(1),...
    'contpar',[ip.a,ip.b],'step',0.01,'dir',ip.a,'print_residual_info',1,'TorusInit.initmethod','svd',...
    'dir',ip.b,'max_step',[],...
    'minimal_angle',0,'plot_measure',[],'print_residual_info',1,'use_tangent',true,...
    'usercond',pd_cond,'initcond',pd_cond,'extracolumns','auto')
pdbr.method.continuation.stops={@(p)p(end).period>50};
figure(2);
pdbr=br_contn(pdfuncs,pdbr,100,'plotaxis',ax2);
pdbr=br_rvers(pdbr);
pdbr=br_contn(pdfuncs,pdbr,100,'plotaxis',ax2);
pdbr=br_stabl(pdfuncs,pdbr,'recompute',1);
%%
[pd_wbifs,pd_tests,pd_bifs,pd_bifind]=MonitorChange(pdfuncs,pdbr,'print_residual_info',0);
[pd_wbifs,nunst_pd,dom,triv]=br_stabl(pdfuncs,pd_wbifs);
i_12=find(real(dom')<0&abs([diff(pd_tests),0])==1);
pd_wbifs.point(i_12)=arrayfun(@(p)setfield(p,'flag','R12'),pd_wbifs.point(i_12));
%% insert period doubling into 2d bifurcation diagram
figure(3);
Plot2dBranch(pd_wbifs,'ax',ax3,'funcs',pdfuncs,'lgname','PD','color',[0,0,0.5]);
%% Branch off at 1st 1:2 resonances
% generate torus bifurcation r.h.s trfuncs, create an initial point with
% fixed rotation number different from 1:2 (omega=1)
[trfuncs,tr0br,suc]=SetupTorusBifurcation(pdfuncs,pd_wbifs,pd_bifind(1),'correc',false)
tr0br.point=pd_wbifs.point(pd_bifind(1));
tr0br.point.parameter(trfuncs.ip.omega)=0.99;
[torusbif,suc]=ChangeBranchParameters(trfuncs,tr0br,1,'print_residual_info',1,...
    'plot_measure',[]);
figure(2)
torusbif=br_contn(trfuncs,torusbif,30,'plotaxis',ax2);
%%
figure(3);
Plot2dBranch(torusbif,'ax',ax3,'funcs',trfuncs,'lgname','TR','color',[0.5,0,0]);
%% trace symmetric periodic orbit with maximum equal to pi/6 in two parameters
% find extrema of $x$  along orbits on sympo branch and pick orbits that
% have two extrema
smaxval=pi/6;
sympo_x1roots=arrayfun(@(p)dde_coll_roots(p,[0,1])',sympo.point,'uniformoutput',false);
x_eval=@(p,t)[1,0]*dde_coll_eva(p.profile,p.mesh,t(:)',p.degree); % evaluate x1 at t in point p
sympomax=cellfun(@(p,t)max(x_eval(p,t)),num2cell(sympo.point(2:end)),sympo_x1roots(2:end));
icross=find(sympomax>smaxval,1,'first');
itcross=find(x_eval(sympo.point(icross+1),sympo_x1roots{icross+1}')>0);
ipe=ip;
ipe.t0=length(fieldnames(ip))+1;
ipe.val=ipe.t0+1;
sympo0=setfield(sympo,'point',sympo.point(icross+1));
sympo0.point.parameter([ipe.t0,ipe.val])=[sympo_x1roots{icross+1}(itcross),smaxval];
max_cond=@(p,pref)dde_extreme_cond(p,[1,0],ipe.val,ipe.t0);
[mfuncs,mbranch,suc]=ChangeBranchParameters(funcs,sympo0,1,...
    'contpar',[ipe.a,ipe.b,ipe.t0],...
    'usercond',{psolsym,max_cond},'outputfuncs',true,...
    'print_residual_info',1,'plot_measure',[]);
figure(2);
mbranch=br_contn(mfuncs,mbranch,200,'plotaxis',ax2);
mbranch=br_rvers(mbranch);
mbranch=br_contn(mfuncs,mbranch,50,'plotaxis',ax2);
[mbranch,nunst_m,dom_m,triv_m]=br_stabl(mfuncs,mbranch,'recompute',1,...
    'max_number_of_eigenvalues',4);
[mbr_wbifs,m_tests,m_bifs,m_bifind]=MonitorChange(mfuncs,mbranch,'print_residual_info',0)
[dum,m_dom,m_triv]=GetStability(mbr_wbifs,'exclude_trivial',true);
%% include family with fixed maximum int o2d bifurcation diagram
figure(3)
Plot2dBranch(mbr_wbifs,'ax',ax3,'funcs',mfuncs,'lgname','max=pi/6','color',[1,0,0],'linestyle','--');
%% Track approximate symmetric heteroclinic orbit by large-period orbit
fixperiod=@(p,pref)sys_cond_coll_fixperiod(p,sympo.point(end).period,false);
[hetfuncs,hetpo,suc]=ChangeBranchParameters(pfuncs,sympo,length(sympo.point),...
    'contpar',[ip.a,ip.b],...
    'max_step',[],...
    'continuation.plot_measure',[],...
    'continuation.stops',[],...
    'usercond',fixperiod,'outputfuncs',true)
hetpo=br_contn(hetfuncs,hetpo,100,'plotaxis',ax2);
hetpo=br_rvers(hetpo);
hetpo=br_contn(hetfuncs,hetpo,30,'plotaxis',ax2);
%% include symmetric heteroclinic into bifurcation diagram
figure(3)
Plot2dBranch(hetpo,'ax',ax3,'funcs',hetfuncs,'lgname','sym. heteroclinic','color',[1,1,1]*0.5,'stability',false);
