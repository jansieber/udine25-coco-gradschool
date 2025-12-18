%% User-Initiated Enhancements
%
% Demos corresponding to Sections 7.1-7.4 in the Getting Started
% tutorial for COCO. For more information about the 'ep' and 'po'
% toolboxes, see EP-Tutorial.pdf and PO-Tutorial.pdf in the coco/help
% folder.

% Copyright (C) 2016-2023 Frank Schilder, Harry Dankowicz, Jan Sieber
clear
format compact
startup_coco(fullfile('..','coco_2025January28'))
%% Section 7.1: Monitoring properties of equilibria

% A model of oxidation of carbon monoxide on platinum analyzed in Tutorial
% IV: Two-parameter bifurcation analysis of equilibria and limit cycles
% with MATCONT, by Yu.A. Kuznetsov, September 20, 2011 (see
% http://www.staff.science.uu.nl/~kouzn101/NBA/LAB4.pdf).

% construct function generator
pnames={'Q1','Q2','Q3','Q4','Q5','Q6','K'};
unames={'x','y','s'};
[iu,ip]=deal(structind_from_names(unames),structind_from_names(pnames));
F=sco_fun({@(u,p)fun_bykov(0,iu,ip,u,p),...
    @(u,p,du,dp)fun_bykov(1,iu,ip,u,p,du,dp)},{'x','p'});
x0([iu.x, iu.y, iu.s])=...
   [0.24; 0.04; 0.51];
p0([ip.Q1,ip.Q2,ip.Q3,  ip.Q4,ip.Q5,ip.Q6,ip.K])=...
     [2.5;  0.5;   10; 0.0675;   1;   0.1; 0.4];
funcs  = {F(''), F('x'), F('p'), F({'x','x'}), F({'x','p'}), F({'p','p'})};

%% continue along a curve of equilibria
prob = coco_prob;
prob = ode_isol2ep(prob, '', funcs{:}, x0, pnames, p0);
eprunid = 'ep_run';
coco(prob, eprunid, [], 'Q2', [0.4 3]);

%% append computed data to the output cell array
prob = ep_add_bddat(prob, '', 'svds', ...
  @(d,x,p) min(svds(feval(F('x'),x,p))));
prob = ep_add_bddat(prob, '', 'sol', ...
  @(d,x,p) ep_create('x',x,'parameter',p));
coco(prob, eprunid, [], 'Q2', [0.4 3]);

%% continue along a curve of saddle-node bifurcations
labs = coco_bd_labs(eprunid, 'SN');
snrunid = 'SN-curve';
coco(snrunid, 'ode', 'SN', 'SN', ...
  eprunid, labs(2), {'Q2' 'K'}, {[0.4 3], [0 8]});

% continue along a curve of Hopf bifurcations
labs = coco_bd_labs(eprunid, 'HB');
prob = coco_prob;
prob = ode_HB2HB(prob, '', eprunid, labs(1));
prob = coco_set(prob, 'cont', 'PtMX', [90 0]);
hbrunid = 'HB-curve';
coco(prob, hbrunid, [], {'Q2' 'K'}, {[0.4 3], [0 8]});

% append computed data to the output cell array
pos_or_nan=@(k)subsref([NaN,k],substruct('()',{1,(k>=0)+1}));
prob = ep_HB_add_bddat(prob, '', 'freq', ...
  @(d,x,p,v,k) sqrt(pos_or_nan(k)));
coco(prob, hbrunid, [], {'Q2' 'K'}, {[0.4 3], [0 8]});

% append additional computed data to the output cell array
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
  'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_bddat(prob, '', 'L1', @lyapunov, 'data', data);
prob = coco_set(prob, 'cont', 'PtMX', [90 0]);
coco(prob, hbrunid, [], {'Q2' 'K'}, {[0.4 3], [0 8]});

% visualize the results
figure(1)
clf
hold on
thm = struct();
thm.special = {'HB', 'SN'};
coco_plot_bd(thm, eprunid, 'Q2', 'x')
coco_plot_bd(snrunid, 'Q2', 'x')
thm = struct();
thm.special = {'BTP'};
thm.ustab = 'L1';
thm.ustabfun = @(x) 1+(~isnan(x) & x>0)+2*isnan(x);
thm.usept = {};
thm.lspec = {{'r-', 'LineWidth', 1.5}, {'r-.', 'LineWidth', 1.5}};
thm.xlab  = 'Q_2';
coco_plot_bd(thm, hbrunid, 'Q2', 'x')
hold off
grid on
axis([0.5 2 0 0.16])

%% Section 7.2: Continuing generalized Hopf bifurcations

% continue along the family of Hopf bifurcations
labs = coco_bd_labs(eprunid, 'HB');
prob = coco_prob;
prob = ode_HB2HB(prob, '', eprunid, labs(2));
% append regular embedded monitor function and detect a zero crossing
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
  'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, ...
  'regular', 'L1');
prob = coco_add_event(prob, 'GH', 'L1', 0);
prob = coco_set(prob, 'cont', 'PtMX', [85 0]);
coco(prob, 'HB-curve', [], {'Q2' 'K', 'L1'}, {[0.4 3], [0 8]});

figure(2)
clf
hold on
thm = struct();
thm.special = {'HB', 'SN'};
coco_plot_bd(thm, eprunid, 'Q2', 'x')
coco_plot_bd('SN-curve', 'Q2', 'x')
thm = struct();
thm.special = {'BTP', 'GH'};
thm.GH = {'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 10};
thm.ustab = 'L1';
thm.ustabfun = @(x) 1+(~isnan(x) & x>0)+2*isnan(x);
thm.usept = {'BTP', 'GH'};
thm.lspec = {{'r-', 'LineWidth', 1.5}, {'r-.', 'LineWidth', 1.5}};
thm.xlab  = 'Q_2';
coco_plot_bd(thm, 'HB-curve', 'Q2', 'x')
hold off
grid on
axis([0.5 2 0 0.16])

% continue along a curve of generalized Hopf bifurcations using an inactive monitor function
labs = coco_bd_labs('HB-curve', 'GH');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'HB-curve', labs(1));
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
  'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, ...
  'inactive', 'L1');
prob = coco_set_parival(prob, 'L1', 0);
prob = coco_set(prob, 'cont', 'PtMX', 50);
coco(prob, 'GH-curve1', [], {'Q2' 'K' 'Q1'}, {[0.4 3], [0 8]});

% continue along a curve of generalized Hopf bifurcations using a zero function
labs = coco_bd_labs('HB-curve', 'GH');
prob = coco_prob;
prob = ode_HB2HB(prob, '', 'HB-curve', labs(2));
data = struct('dfdx', F('x'), 'Dfdxdx', F({'x*v','x*v'}), ...
  'Dfdxdxdx', F({'x*v','x*v','x*v'}), 'nanflag', 1);
prob = ep_HB_add_func(prob, '', 'lyap', @lyapunov, data, 'zero');
prob = coco_set(prob, 'cont', 'PtMX', 50);
coco(prob, 'GH-curve3', [], {'Q2' 'K' 'Q1'}, {[0.4 3], [0 8]});


%% Section 7.3: Monitoring properties of periodic orbits

% construct function (generalized Hopf normal form)
r2=@(x)x(1,:).^2+x(2,:).^2;
f=@(x,p)[...
    x(1,:).*(p(1,:)+p(2,:).*r2(x)-r2(x).^2)-x(2,:);...
    x(2,:).*(p(1,:)+p(2,:).*r2(x)-r2(x).^2)+x(1,:)];
% initial run for finding branch of equilibria
prob = coco_prob;
prob = ode_isol2ep(prob, '', f, [0; 0], {'p1', 'p2'}, [-1; 1]);
coco(prob, 'ep_run', [], 'p1', [-1 1]);

% continue along a family of periodic orbits from a Hopf bifurcation
HB = coco_bd_labs('ep_run', 'HB');
prob = coco_prob;
prob = ode_HB2po(prob, '', 'ep_run', HB);
prob = coco_set(prob, 'cont', 'PtMX', [50 0],'NAdapt',1);
prob = po_add_bddat(prob, '', 'cx1', @fourier, 'data', struct('n', 1));
coco(prob, 'po_run', [], {'p1' 'p2'}, [-1 1]);

% continue along a family of saddle-node bifurcations of periodic orbits
SN = coco_bd_labs('po_run', 'SN');
prob = coco_prob;
prob = ode_po2SN(prob, '', 'po_run', SN(1));
prob = po_add_bddat(prob, '', 'cx1', @fourier, 'data', struct('n', 1));
coco(prob, 'po_SN_run', [], {'p1', 'p2'}, {[-1 1] [0.001 3]});

figure(1)
clf
hold on
coco_plot_bd('po_run', 'p1', 'p2', 'cx1', @(x) 2*real(x(1,:)))
thm = struct();
thm.lspec = {{'ro', 'MarkerFaceColor', 'r'}, {'ro'}};
coco_plot_bd(thm, 'po_run', ...
  {'cx1', 'p2'}, @(x,y) (2*real(x(1,:))).^4-y.*(2*real(x(1,:))).^2, ...
  'p2', 'cx1', @(x) 2*real(x(1,:)));
coco_plot_bd('po_SN_run', 'p1', 'p2', 'cx1', @(x) 2*real(x(1,:)))
thm.lspec = {'ko', 'MarkerFaceColor', 'k'};
thm.xlab = 'p1';
thm.ylab = 'p2';
thm.zlab = 'r';
coco_plot_bd(thm, 'po_SN_run', ...
  'cx1', @(x) -(2*real(x(1,:))).^4, ...
  'cx1', @(x) 2*(2*real(x(1,:))).^2, ...
  'cx1', @(x) 2*real(x(1,:)));
hold off
grid on
view(3)

%% Section 7.4: Tracking and constraining orbit maxima

% construct function for Marsden problem
f=@(x,p)[...
    p(1,:).*x(1,:)+x(2,:)+p(2,:).*x(1,:).^2;...
    -x(1,:)+p(1,:).*x(2,:)+x(2,:).*x(3,:);...
    (p(1,:).^2-1).*x(2,:)-x(1,:)-x(3,:)+x(1,:).^2];
prob = coco_prob;
prob = ode_isol2ep(prob, '', f, [0; 0; 0], {'p1', 'p2'}, [-1; 6]);
coco(prob, 'ep_run', [], 'p1', [-1 1]);

% continue along a family of periodic orbits from a Hopf bifurcation
HB = coco_bd_labs('ep_run', 'HB');
prob = coco_prob;
prob = ode_HB2po(prob, '', 'ep_run', HB);
prob = coco_set(prob, 'cont', 'PtMX', [50 0], 'NAdapt', 1);
coco(prob, 'po_run', [], {'p1' 'po.period'}, [-1 1]);

figure(1)
clf
coco_plot_sol('po_run', '', 't', 'x')
grid on
axis([0 inf -0.3 0.25])

% continue along family of periodic orbits while tracking x1 and x1' for fixed t
prob = coco_prob;
prob = ode_po2po(prob, '', 'po_run', 5);
prob = coco_set(prob, 'cont', 'PtMX', [0 50], 'NAdapt', 1);
prob = po_add_func(prob, '', 'cmax', @slope,@dslope, struct('idx', 1), ...
  {'trs', 'xrs', 'xtrs'}, 'inactive', 'u0', 5.17);
coco(prob, 'new_po_run1', [], {'p1' 'po.period','xrs', 'xtrs'}, [-1 1]);
figure(2)
clf
hold on
coco_plot_sol('new_po_run1', '', 't', 'x')
tx = coco_bd_vals('new_po_run1', 'all', {'trs', 'xrs'});
plot(tx(1,:), tx(2,:), 'ro', 'MarkerFaceColor', 'r');
hold off
grid on
box on
axis([0 inf -0.3 0.25])

% continue along family of periodic orbits while tracking maximum in x1
prob = coco_set_parival(prob, 'xtrs', 0);
coco(prob, 'new_po_run2', [], {'p1' 'po.period', 'trs', 'xrs'}, [-1 1]);
figure(3)
clf
hold on
coco_plot_sol('new_po_run2', '', 't', 'x')
tx = coco_bd_vals('new_po_run2', 'all', {'trs', 'xrs'});
plot(tx(1,:), tx(2,:), 'ro', 'MarkerFaceColor', 'r');
hold off
grid on
box on
axis([0 inf -0.3 0.25])

% continue along family of periodic orbits with fixed maximum in x1
coco(prob, 'new_po_run3', [], {'p1' 'po.period', 'trs', 'p2'}, [-1 1]);
figure(4)
clf
hold on
coco_plot_sol('new_po_run3', '', 't', 'x')
tx = coco_bd_vals('new_po_run3', 'all', {'trs', 'xrs'});
plot(tx(1,:), tx(2,:), 'ro', 'MarkerFaceColor', 'r');
hold off
grid on
box on
axis([0 inf -0.3 0.25])
