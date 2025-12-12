%% A web of bifurcation curves for a three-dimensional vector field
%
% We explore bifurcations of equilibria and periodic orbits in a vector
% field described in Freire, E., Rodriguez-Luis, A., Gamero, E. & Ponce, E.
% (1993), ''A case study for homoclinic chaos in an autonomous electronic
% circuit: A trip from Takens-Bogdanov to Hopf-Shilnikov,'' Physica D 62,
% 230?253, and also used in the AUTO manual to demonstrate functionality
% and syntax.

% Figure shows trival branch of equilibria at the origin under variations
% in one problem parameter, the detection and continuation under variations
% in two problem parameters of a family of Hopf bifurcations of the trivial
% equilibrium, a primary branch of periodic orbits emanating from the Hopf
% bifurcation under variations in one problem parameter, a secondary branch
% of periodic orbits emanating from a branch point along the primary
% branch, continuation of families of saddle-node, period-doubling and
% torus bifurcations under variations in two problem parameters from points
% along the secondary branf of periodic orbits, and continuation of a
% branch of period-doubled periodic orbits under variation in one problem
% parameter emanating from a point on the period-doubling bifurcation
% branch.

%% Continue branch of equilibria at origin

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation.
clear
format compact
% execute coco's startup for setting path
startup_coco(fullfile(pwd(),'..','coco_2025January28'));
addpath(fullfile(pwd(),'tools'))
parnames={'nu','be','ga','r','a3','b3'};
ip=structind_from_names(parnames);
%% Symbolic toolbox (see gentor.m) to provide derivative
%tor=sco_gen(@sym_tor);
%funcs  = {tor(''),tor('x'),tor('p'),tor({'x','x'}),tor({'x','p'}),tor({'p','p'})};
%% Do not provide derivative
f= @(x,nu,be,ga,r,a3,b3)[...
    ( -(be+nu).*x(1,:) + be.*x(2,:) - a3.*x(1,:).^3 + b3.*(x(2,:)-x(1,:)).^3 )./r;
        be.*x(1,:) - (be+ga).*x(2,:) - x(3,:) - b3.*(x(2,:)-x(1,:)).^3;
        x(2,:)];
tor=@(x,p)f(x,p(ip.nu,:),p(ip.be,:),p(ip.ga,:),p(ip.r,:),p(ip.a3,:),p(ip.b3,:));
funcs={tor};
%% Initial parameters and equilibrium
p0([ip.nu, ip.be,ip.ga,  ip.r,ip.a3,ip.b3]) =...
   [-0.65 ; 0.5 ; -0.6 ; 0.6 ; 0.3 ; 0.9];
p0=p0(:);
%% Track equilibria in nu
prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, [0;0;0],parnames, p0);
fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep');

bd_ep = coco(prob, 'ep', [], 1, 'nu', [-0.65, -0.55]);
%% Plot result
figure(1); clf; hold on; grid on
set(gca,'FontSize',18);
thm = struct('special', {{'EP', 'HB'}});
coco_plot_bd(thm, 'ep', 'nu', 'be', '||x||_2')
%% Continue branch of Hopf bifurcations of equilibrium at origin

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released and allowed to vary
% during continuation.
bd_ep=coco_bd_read('ep');
lab = coco_bd_labs(bd_ep, 'HB');
prob = coco_prob();
prob = ode_HB2HB(prob, '', 'ep', lab);

fprintf(...
  '\n Run=''%s'': Continue Hopf bifurcations from point %d in run ''%s''.\n', ...
  'hb', lab, 'ep');

bd_hb = coco(prob, 'hb', [], 1, {'nu', 'be'}, [-0.65, -0.55]);
%% Add result to plot
thm = struct('special', {{'EP', 'BP', 'FP'}});
coco_plot_bd(thm, 'hb', 'nu', 'be', '||x||_2')
%% Continue first branch of periodic orbits

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period, the
% ratio of the interpolation error estimate with the tolerance, and the
% stability indicator that counts the number of Floquet multipliers outside
% the unit circle.
bd_ep=coco_bd_read('ep');
lab = coco_bd_labs(bd_ep, 'HB');
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false);
prob = ode_HB2po(prob, '', 'ep', lab);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', [100 0]);
prob=coll_simple_measure(prob,'po.orb.coll','nrmx');
prob = coco_add_event(prob, 'MXNRM', 'MX', 'po.orb.coll.nrmx', '>', 1);
fprintf(...
  '\n Run=''%s'': Continue primary branch of periodic orbits from point %d in run ''%s''.\n', ...
  'po1', lab, 'ep');

bd1  = coco(prob, 'po1', [], 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF' 'po.orb.coll.nrmx' 'po.test.USTAB'}, ...
  {[-0.65, -0.55]});
%% Add result to plot
thm = struct('special', {{'EP', 'SN', 'BP', 'FP'}});
coco_plot_bd(thm, 'po1', 'nu', 'be', '||x||_{2,MPD}')
%% Branch off at BP from po1
% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period, the
% ratio of the interpolation error estimate with the tolerance, and the
% stability indicator that counts the number of Floquet multipliers outside
% the unit circle.
bd1=coco_bd_read('po1');
BPlabs = coco_bd_labs(bd1, 'BP');
prob = coco_prob();
prob=ode_po2po(prob,'','po1',BPlabs(1));
prob=coco_set(prob,'cont','branch','switch');
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', [100 0]);

fprintf(...
  '\n Run=''%s'': Continue secondary branch of periodic orbits from point %d in run ''%s''.\n', ...
  'po2', BPlabs(end), 'po1');

bd2 = coco(prob, 'po2', [], 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}, [-0.65, -0.55]);
%% Add result to plot
thm = struct('special', {{'EP', 'PD', 'TR'}});
coco_plot_bd(thm, 'po2', 'nu', 'be', '||x||_{2,MPD}')

%% Continuation of saddle-node bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.
bd1=coco_bd_read('po1');
SNlabs = coco_bd_labs(bd1, 'SN');
prob = coco_prob();
prob=ode_SN2SN(prob,'','po1',SNlabs(2));
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue saddle-node bifurcations from point %d in run ''%s''.\n', ...
  'lp1', SNlabs(2), 'po1');

bd3 = coco(prob, 'lp1',[],1, {'nu' 'po.period' 'be'}, [-0.65, -0.55]);
%% Add result to plot
thm = struct('special', {{'EP', 'FP'}});
coco_plot_bd(thm, 'lp1', 'nu', 'be', '||x||_{2,MPD}')

%% Continuation of period-doubling bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.

bd2=coco_bd_read('po2');
PDlabs = coco_bd_labs(bd2, 'PD');
prob = coco_prob();
prob=ode_PD2PD(prob,'','po2',PDlabs(1));
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue period-doubling bifurcations from point %d in run ''%s''.\n', ...
  'pd1', PDlabs(1), 'po2');

bd4 = coco(prob, 'pd1',[],1, {'nu' 'po.period' 'be'}, [-0.65, -0.55]);
%% Add result to plot
thm = struct('special', {{'EP', 'FP'}});
coco_plot_bd(thm, 'pd1', 'nu', 'be', '||x||_{2,MPD}')

%% Continuation of torus bifurcations

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released ('po.period' is
% already active) and allowed to vary during continuation.

bd2=coco_bd_read('po2');
TRlabs = coco_bd_labs(bd2, 'TR');
prob = coco_prob();
prob=ode_TR2TR(prob,'','po2',TRlabs(1));
prob = coco_set(prob, 'cont', 'NAdapt', 5);
%% add a variable and continuation parameter recording the rotation number
[tr.uidx,tr.u0,tr.data]=coco_get_func_data(prob,'po.TR','uidx','u0','data');
mon_alpha=@(a,b,alpha)b*cos(2*pi*alpha)-a*sin(2*pi*alpha);
Jmon_alpha=@(a,b,alpha)[-sin(2*pi*alpha),cos(2*pi*alpha),...
    2*pi*(-b*sin(2*pi*alpha)-a*cos(2*pi*alpha))];
ab=tr.u0([tr.data.po_tr.a_idx,tr.data.po_tr.b_idx]);
prob=coco_add_func(prob,'roteq',...
    f2coco(@(u)mon_alpha(u(1),u(2),u(3))),...    
    f2coco(@(u)Jmon_alpha(u(1),u(2),u(3))),...
    [],'zero',...
    'uidx',tr.uidx([tr.data.po_tr.a_idx,tr.data.po_tr.b_idx]),...
    'u0',atan2(ab(2),ab(1))/2/pi);
rot_uidx=coco_get_func_data(prob,'roteq','uidx');
prob=coco_add_pars(prob,'tr_rot',rot_uidx(3),{'rotation'});
%% add events for encounter of low-order resonances
farey_seq=generate_farey(13);
farey_signed=[diag([-1,1])*farey_seq(:,end:-1:2),farey_seq(:,2:end)];
prob=coco_add_event(prob,'RE','rotation',farey_signed(1,:)./farey_signed(2,:));
%%
fprintf(...
  '\n Run=''%s'': Continue torus bifurcations from point %d in run ''%s''.\n', ...
  'tr1', TRlabs(1), 'po2');
bd5 = coco(prob, 'tr1',[],1, {'nu' 'po.period' 'be' 'rotation' }, [-0.65, -0.55]);
%% Add result to plot
thm=coco_plot_theme();
thm.special={'EP', 'FP','RE'};
thm.RE={'o' 'Color' [0,0.5,0] 'MarkerFaceColor'  'y'  'MarkerSize'  [8]};
coco_plot_bd(thm, 'tr1', 'nu', 'be', '||x||_{2,MPD}')

%% Continuation along period-doubled branch

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'nu' is released and allowed to vary during
% continuation. The screen output also includes the orbital period and the
% ratio of the interpolation error estimate with the tolerance.

prob = coco_prob();
prob=ode_PD2po(prob,'','pd1',4); % use
prob = coco_set(prob, 'cont', 'NAdapt', 5);

fprintf(...
  '\n Run=''%s'': Continue period-doubled periodic orbits from point %d in run ''%s''.\n', ...
  'po_db', 4, 'pd1');

bd_db  = coco(prob, 'po_db',[], 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF'}, [-0.65, -0.55]);
%% Add result to plot
thm = struct('special', {{'EP', 'BP', 'FP', 'PD'}}, 'xlab', '\nu', 'ylab', '\beta');
coco_plot_bd(thm, 'po_db', 'nu', 'be', '||x||_{2,MPD}')
hold off; view(-192,44)
