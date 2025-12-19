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
% execute coco's startup for setting path
startup_coco(fullfile('..','coco_2025January28'));
addpath(fullfile('tools'));
parnames={'nu','be','ga','r','a3','b3'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
tor=sco_gen(@sym_tor);
funcs  = {tor(''),tor('x'),tor('p'),tor({'x','x'}),tor({'x','p'}),tor({'p','p'})};
p0([ip.nu, ip.be,ip.ga,  ip.r,ip.a3,ip.b3]) =...
   [-0.65 ; 0.5 ; -0.6 ; 0.6 ; 0.3 ; 0.9];

prob = coco_prob();
prob = ode_isol2ep(prob, '', funcs{:}, [0;0;0],parnames, p0);


fprintf('\n Run=''%s'': Continue family of equilibrium points.\n', ...
  'ep');

bd_ep = coco(prob, 'ep', [], 1, 'nu', [-0.65, -0.55]);

%% Continue branch of Hopf bifurcations of equilibrium at origin

% The continuation problem encoded below includes six monitor functions
% that evaluate to the problem parameters, and corresponding inactive
% continuation parameters 'nu', 'be', 'ga', 'r', 'a3', and 'b3'. Its
% dimensional deficit equals -1. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameters 'nu' and 'be' are released and allowed to vary
% during continuation.

lab = coco_bd_labs(bd_ep, 'HB');
prob = coco_prob();
prob = ode_HB2HB(prob, '', 'ep', lab);

fprintf(...
  '\n Run=''%s'': Continue Hopf bifurcations from point %d in run ''%s''.\n', ...
  'hb', lab, 'ep');

bd_hb = coco(prob, 'hb', [], 1, {'nu', 'be'}, [-0.65, -0.55]);

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

prob = coco_prob();
prob = ode_HB2po(prob, '', 'ep', lab);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'PtMX', [10 0]);

fprintf(...
  '\n Run=''%s'': Continue primary branch of periodic orbits from point %d in run ''%s''.\n', ...
  'po1', lab, 'ep');

bd1  = coco(prob, 'po1', [], 1, ...
  {'nu' 'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}, [-0.65, -0.55]);

%% try continuation in two parameters with small maximal stepsize
bd1=coco_bd_read('po1');
EPlabs = coco_bd_labs(bd1, 'EP');
coco_func_data.pointers('set',[]);
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false);
prob=ode_po2po(prob,'','po1',EPlabs(1));
prob = coco_set(prob, 'cont', 'NAdapt', 5,'NSV',1, 'PtMX', 3200,'R_max',0.025);
[prob,mod1ids]=coll_construct_modes(prob,'po.orb.coll',...
    'fid','modes','mid','mode','norm',true);
bd2  = coco(prob, 'po2d_dense', [], 2, ...
  [{'nu' 'be'}, mod1ids, {'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}], ...
  {[-0.65, -0.55],[0.45,0.65],[3e-2,1]});
%% plot atlas as trimesh
ltx={'Interpreter','LaTeX'};
[info,po_atlas,bd2]=info_from_run('po2d_dense');
xp=info.p;
cind=[info.pnames;num2cell(1:length(info.pnames))];
acp=@(s)cind{2,strcmp(s,cind(1,:))};
[tri,adj]=triangles_from_atlas('atlas',po_atlas);
figure(3);clf;hold on;ax3=gca;
colormap('parula')
trisurf(tri,...
    xp(:,acp('nu')),xp(:,acp('be')),xp(:,acp('po.orb.coll.mode_sum')),...
    xp(:,acp('po.period')),'linewidth',0.5,'FaceAlpha',0.5)
grid(ax3,'on');
xlabel(ax3,'$\nu$',ltx{:});
ylabel(ax3,'$\beta$',ltx{:});
zlabel(ax3,'$\sqrt{\sum \mathrm{modes}^2}$',ltx{:});
cb=colorbar(ax3);
cb.Label.String='period';
set(ax3,'Fontsize',18,'fontweight','bold','fontname','courier');
%% add stability of points
ustab=info.p(:,strcmp(info.pnames,'po.test.USTAB'));
ustab=ustab-min(ustab)+1;
clr=lines();
for i=1:length(ustab)
    plot3(ax3,xp(i,acp('nu')),xp(i,acp('be')),xp(i,acp('po.orb.coll.mode_sum')),...
    'o','color',clr(ustab(i),:),'MarkerFaceColor',clr(ustab(i),:));
end
%% add BP's and SN's (detected as SN's)
cpply=@(f,x)cell2mat(cellfun(f,x,'uniformoutput',false));
apply=@(f,x)cell2mat(arrayfun(f,x,'uniformoutput',false));
snlabs=coco_bd_lab2idx(bd2, coco_bd_labs(bd2,'SN'));
sncharts=apply(@(x)coco_read_solution('po2d_dense',x,'chart'),snlabs);
%%
lw={'linewidth',2};
bd2t=coco_bd_table('po2d_dense');
bpvals=bd2t{strcmp(bd2t.TYPE,'SN'),{'nu','be','po.orb.coll.mode_sum'}};
%bpvals=cpply(@(n)coco_bd_vals(bd2,coco_bd_labs(bd2,'SN'),n),...
%    {'nu','be','po.orb.coll.mode_sum'});
plot3(ax3,bpvals(:,1),bpvals(:,2),bpvals(:,3),...
    'ko','MarkerFaceColor','k',lw{:});

%% try continuation in two parameters with larger maximal stepsize
bd1=coco_bd_read('po1');
EPlabs = coco_bd_labs(bd1, 'EP');
coco_func_data.pointers('set',[]);
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL',false);
prob=ode_po2po(prob,'','po1',EPlabs(1));
prob = coco_set(prob, 'cont', 'NAdapt', 5,'NSV',1, 'PtMX', 3200,'R_max',0.1);
[prob,mod1ids]=coll_construct_modes(prob,'po.orb.coll',...
    'fid','modes','mid','mode','norm',true);
bd2a  = coco(prob, 'po2d', [], 2, ...
  [{'nu' 'be'}, mod1ids, {'po.period' 'po.orb.coll.err_TF' 'po.test.USTAB'}], ...
  {[-0.65, -0.55],[0.45,0.65],[3e-2,1]});
%% plot atlas as trimesh
ltx={'Interpreter','LaTeX'};
[info2,po_atlas2,bd2a]=info_from_run('po2d');
xp=info2.p;
cind2=[info2.pnames;num2cell(1:length(info2.pnames))];
acp2=@(s)cind2{2,strcmp(s,cind2(1,:))};
[tri2,adj2]=triangles_from_atlas('atlas',po_atlas2);
figure(4);clf;hold on;ax4=gca;
trisurf(tri2,...
    xp(:,acp2('nu')),xp(:,acp2('be')),xp(:,acp2('po.orb.coll.mode_sum')),...
    xp(:,acp2('po.period')),'linewidth',0.5,'FaceAlpha',0.5)
grid(ax4,'on');
xlabel(ax4,'$\nu$',ltx{:});
ylabel(ax4,'$\beta$',ltx{:});
zlabel(ax4,'$\sqrt{\sum \mathrm{modes}^2}$',ltx{:});
cb=colorbar(ax4);
cb.Label.String='period';
set(ax4,'Fontsize',18,'fontweight','bold','fontname','courier');
plot_atlas_kd(po_atlas2.charts,1,2,3,'basepoints')