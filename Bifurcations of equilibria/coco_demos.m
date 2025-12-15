%% Bifurcation Diagrams - Equilibria
%
% Demos corresponding to Sections 3.1-3.5 in the Getting Started tutorial
% for COCO. For more information about the 'ep' toolbox, see
% EP-Tutorial.pdf in the coco/help folder. All the demos below use the
% inline construction of a continuation problem in the call to the coco
% entry-point function.

% Copyright (C) 2016 Frank Schilder, Harry Dankowicz, Jan Sieber
%% setup path
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))

%% Section 3.1: Detecting saddle-node points

f = @(x,p) p-x.^2; % vectorized form of the fold normal form

% continue along a single curve of equilibria
% initial solution for equilibrium point
% provide equation name, vector field, initial state value, parameter name, initial parameter value
prob=coco_prob();
prob=ode_isol2ep(prob,'',f,0.5,'p',1); 
coco(prob,'run1',[],'p', [-1 1]); % branch label, [], branch dimension, active continuation parameter, computational domain

% visualize the result
figure(1)
clf
theme = struct('special', {{'SN', 'FP'}}); % plotting theme (check with ep_plot_theme())
coco_plot_bd(theme, 'run1', 'p', 'x') % 'x' denotes the state vector, which is scalar here
axis tight
grid on
% print out bifurcation diagram
bd_run1=coco_bd_table('run1','numlab',true) %#ok<*NOPTS>
%% Section 3.2: Detecting branch points

f = @(x,p) p.*x-x.^2; % vectorized form of the transcritical normal form

% continue along a first curve of equilibria
prob=coco_prob();
prob=ode_isol2ep(prob,'', f, 0, 'p', -1);
coco(prob,'run2a', [],'p', [-1 1]); 

% continue along a second curve of equilibria through a branch point
BP = coco_bd_labs('run2a', 'BP'); % labels for BP points in run1
prob=coco_prob();
prob=ode_ep2ep(prob,'','run2a',BP);    % reload solution from label
prob = coco_set(prob, 'cont', 'branch','switch'); % set parameters:enforce switching
coco(prob, 'run2b',[],'p', [-1 1]); % branch label, toolbox family, initial point, branch type
  

% visualize the results
figure(2)
clf
theme1 = struct('special', {{'BP', 'EP'}});
theme2 = struct('special', {{'EP'}});
hold on
coco_plot_bd(theme1, 'run2a', 'p', 'x')
coco_plot_bd(theme2, 'run2b', 'p', 'x')
hold off
axis tight
grid on
% look up files coco_screen.txt, coco_log.txt in folders data/run2a, data/run2b
%% Section 3.3: Continuing saddle-node points

f = @(x,p) p(2,:)+x.*(p(1,:)-x.^2); % vectorized form of the pitchfork normal form

% continue along a first curve of equilibria
prob=coco_prob();
prob=ode_isol2ep(prob,'', f, 0, {'la' 'ka'}, [-2; 0]);
coco(prob,'run3a',[],'la', [-2 2]);

% continue along a second curve of equilibria through a branch point
BP = coco_bd_labs('run3a', 'BP');
prob=coco_prob();
prob=ode_ep2ep(prob,'','run3a', BP);
prob = coco_set(prob, 'cont', 'branch', 'switch');
coco(prob, 'run3b',[], 'la', [-2 2]);

% start from other point, continue along a third curve of equilibria
prob=coco_prob();
prob=ode_isol2ep(prob,'', f, 0,{'la' 'ka'}, [0.5; 0]);
coco(prob,'run3c',[], 'ka', [-2 2]);

% continue along a curve of saddle-node bifurcations of equilibria in two
% parameters
SN = coco_bd_labs('run3c', 'SN');
prob=coco_prob();
prob=ode_ep2SN(prob,'','run3c', SN(1));
coco(prob,'run3d',[], {'la' 'ka'}, {[-2 2] [-2 2]});
bd_run3d=coco_bd_table('run3d','numlab',true) %#ok<*NOPTS>

% visualize the results
figure(3)
clf
theme1 = struct('special', {{'BP', 'EP'}});
theme2 = struct('special', {{'EP'}});
theme3 = struct('special', {{'SN' 'EP'}});
hold on
coco_plot_bd(theme1, 'run3a', 'la', 'ka', 'x')
coco_plot_bd(theme2, 'run3b', 'la', 'ka', 'x')
coco_plot_bd(theme3, 'run3c', 'la', 'ka', 'x')
coco_plot_bd(theme2, 'run3d', 'la', 'ka', 'x')
hold off
axis tight
grid on
view(-15,25)
%% Section 3.4: Continuing Hopf bifurcations in the Brusselator model
% vectorized encoding of vector field in brus.m

% continue along a curve of equilibria
prob=coco_prob();
prob=ode_isol2ep(prob,'', @brus, [1; 0], {'A' 'B'}, [1; 0]);
coco(prob,'run4a',[],'B', [0 3]); 

% continue along a curve of Hopf bifurcations of equilibria
HB = coco_bd_labs('run4a', 'HB');
prob=coco_prob();
prob=ode_ep2HB(prob,'','run4a',HB(1));
coco(prob,'run4b',[], {'A' 'B'}, {[] [0 3]});

% visualize the results
figure(4)
clf
hold on
theme1 = struct('special', {{'HB', 'EP'}});
coco_plot_bd(theme1, 'run4a', 'B', 'A', '||x||_2') % ||x||_2 denotes the Euclidean norm of the state vector
theme2 = struct('special', {{'EP'}});
coco_plot_bd(theme2, 'run4b', 'B', 'A', '||x||_2')
hold off
axis tight
grid on
view(-65,40)

%% Section 3.5: Performing parameter sweeps in the ABC reaction
% vectorized encoding of vector field in abc.m

% continue along multiple families of equilibria and visualize the results
figure(5)
clf
hold on
axis([0.12 0.22 1 7])
grid on
Nrun = 0;
theme = struct('special', {{'SN' 'HB'}});
for beta = linspace(1.20, 1.42, 23)
  Nrun = Nrun+1;
  runid = sprintf('run5_beta_%d', Nrun);
  prob=coco_prob();
  prob=ode_isol2ep(prob,'',@abc,[0;0;0],...
      {'al', 'si', 'D', 'B', 'be'},...
        [1; 0.04;  0.0;  8;  beta]);
  coco(prob,runid,[], {'D' 'be'}, [0.0 0.25]);
  coco_plot_bd(theme, runid, 'D','||x||_2','be');
  drawnow
end

% continue along a curve of saddle-node bifurcations and visualize the result
runid = 'run5_beta_1';
SN = coco_bd_labs(runid, 'SN');
prob=coco_prob();
prob=ode_ep2SN(prob,'',runid,SN(1));
coco(prob,'run5_SN',[], {'D' 'be'}, [0.1 0.25]);
theme_SN = struct('special', {{'FP'}});
coco_plot_bd(theme_SN, 'run5_SN', 'D','||x||_2','be');

% continue along a curve of Hopf bifurcations and neutral saddles and visualize the result
HB = coco_bd_labs(runid, 'HB');
prob=coco_prob();
prob=ode_ep2HB(prob,'',runid, HB(1));
% set continuation parameters look up parameter meaning in
% atlas_1d_settings(), ep_settings(), coll_settings(), ode_settings(),
% po_settings(),...
prob = coco_set(prob, 'cont', 'ItMX', 150);
coco(prob, 'run_HB',[],{'D' 'be'}, [0.1 0.8]);
theme_HB = struct('special', {{'FP' 'BTP'}});
coco_plot_bd(theme_HB, 'run_HB', 'D','||x||_2','be');
