%% Neural network model with delay
% 
% Demo illustrating how to branch off a genric Bogdanov-Takens bifurcation point
%

%% Differential equations neural network model with delay
% 
% From: Giannakopoulos, Fotios and Zapp, Andreas
% Bifurcations in a planar system of differential delay equations modeling 
% neural activity
% Physica D: Nonlinear Phenomena, 159 (2001), no 3, 215-232.
%
% $$\mu\dot{u}_1(t) = -u_1(t) + q_{11}\alpha(u_1(t\text{-}T))-q_{12}u_2(t\text{-}T) + e_1,$$
%
% $$u_2'=\mu\dot{u}_2(t) = -u_2(t) + q_{21}\alpha(u_1(t\text{-}T))-q_{22}u_2(t\text{-}T) + e_2.$$
%

%% Clean workspace and add DDE-BifTool scripts to the MATLAB search path
clear;      % clear variables
close all;	% close figures
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_psol'),...
        strcat(ddebiftoolpath,'ddebiftool_extra_nmfm'),...
        strcat(ddebiftoolpath,'ddebiftool_utilities'),...
        '../');

%% Set parameter names
% This could also be loaded from mat file |'FHN_parnames.mat'|.
parnames={'Q','E','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});

%% Set the funcs structure
% We load the precalculated multilinear forms. These have been generated
% with the file gen_sym_neural_network.m.
funcs=set_symfuncs(@sym_neural_network_mf,'sys_tau',@()ind.tau);
% Alternatively, uncomment the line below to use directional derivaties.
% funcs=set_symfuncs(@sym_FHN,'sys_tau',@()ind.tau);

%% Set bifurcation parameter range and step size bounds
brpars={'max_bound', [ind.Q 1.7],...
        'min_bound', [ind.Q 1.2; ind.E 0],...
        'max_step', [ind.Q 0.005; ind.E 0.005]};

%% Define analytically derived Bogdanov-Takens point
% construct steady-state point
Q0 = 13/10; E0 = (sqrt(39) - 10*atanh(sqrt(3/13)))/20; tau0 = 1;
stst=dde_stst_create('x',[1/4*log((8-sqrt(39))/5); -1/2*sqrt(3/13)]);
stst.parameter(ind.Q)  = Q0;
stst.parameter(ind.E) = E0;
stst.parameter(ind.tau)   = tau0;
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1

%% Calculate coefficients
free_pars=[ind.Q, ind.E];
bt = p_tobt(funcs,stst);
bt = nmfm_bt_orbital(funcs, bt,'free_pars', free_pars);

%% Plot comparing profiles computed and predicted homoclinic orbits
figure(1); clf; hold on;
cm=colormap('lines');
tiledlayout(2,2)
ylabels = {'$x$','$y$'};
step = 0.2;
for order=[1,3]
    [~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
        'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0);
    [~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
        'codim2','BT','codim1','hcli','step',step,'debug',0);
    for j=1:2
        ax = nexttile; hold on; title(sprintf('order %d',order))
        plot(ax,hcli_br_approx.point(1).mesh,hcli_br_approx.point(1).profile(j,:),'*')
        plot(ax,hcli_br_correc.point(1).mesh,hcli_br_correc.point(1).profile(j,:))
        xlabel('$t$','Interpreter','LaTeX'); 
        ylabel(string(ylabels(j)),'Interpreter','LaTeX')
    end
end

%% Continue homoclinic solutions emanating from the BT point
[~,hcli_br,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli',brpars{:},'plot',0,'step',0.1);
assert(all(suc(:)>0))
nop=300; [hcli_br,suc]=br_contn(funcs,hcli_br,nop);assert(suc>0)

%% Continue Hopf curve
[~,hbr,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf',brpars{:},'step',1e-05,'plot',0);
assert(all(suc(:)>0))
nop=300; [hbr,suc]=br_contn(funcs,hbr,nop);assert(suc>0)

%% Continue fold curve emanating from Bogdanov-Takens point
[~,fold_br,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold',brpars{:},'step',1e-03,'plot',0);
assert(all(suc(:)>0))
nop=300; [fold_br,suc]=br_contn(funcs,fold_br(1),nop);assert(suc>0)
fold_br = br_rvers(fold_br);
nop=300; [fold_br,suc]=br_contn(funcs,fold_br(1),nop);assert(suc>0)

%% Plot the bifurcation curves
figure(2);clf;hold on
title('Codimension 1 curves emanating from the transcritical Bogdanov-Takens point');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hcli_br_pm = [getpars(hcli_br,ind.Q), getpars(hcli_br,ind.E)]';
hbr_pm = [getpars(hbr,ind.Q); getpars(hbr,ind.E)];
fold_br_pm = [getpars(fold_br,ind.Q); getpars(fold_br,ind.E)];
plot(hcli_br_pm(1,:),hcli_br_pm(2,:),'.','Color',cm(1,:),'DisplayName','Homoclinic branches')
plot(hbr_pm(1,:),hbr_pm(2,:),'Color',cm(2,:),'DisplayName','Hopf branches')
plot(fold_br_pm(1,:),fold_br_pm(2,:),'Color',cm(3,:),'DisplayName','Transcritical branch')
h(1) = plot(bt.parameter(ind.Q),bt.parameter(ind.E),'.k', 'MarkerSize', 18, ...
        'DisplayName','generic Bogdanov-Takens point');
xlabel('$Q$','Interpreter','LaTeX')
ylabel('$E$','Interpreter','LaTeX')
% axis([1.2928    1.4323    0.0209    0.0522])
legend

%% Predictors for homoclinic, Hopf, and transcritical curves
step = linspace(1e-4,0.6,50);
[~,hcli_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,'debug',0,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true);
step = linspace(0.0001,0.2,50);
[~,hbr_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf','step',step,'predictor',true);
step = linspace(-0.1,0.3,40);
[~,fold_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold','step',step,'predictor',true);
hcli_br_pred_pm  = [getpars(hcli_br_pred, ind.Q), getpars(hcli_br_pred, ind.E)];
hbr_pm_pred = [getpars(hbr_pred,ind.Q); getpars(hbr_pred,ind.E)];
fold_br_pm_pred = [getpars(fold_br_pred,ind.Q); getpars(fold_br_pred,ind.E)];

%% Plot the bifurcation curves with predictors
figure(3);clf;hold on
title(['Codimension 1 curves emanating from the generic ', ... 
       'Bogdanov-Takens point with predictors']);
plot(hcli_br_pm(1,:),hcli_br_pm(2,:),'.','Color',cm(1,:),'DisplayName','Homoclinic branches')
plot(hcli_br_pred_pm(:,1),hcli_br_pred_pm(:,2),'Color', cm(7,:),'DisplayName','Third order homoclinic asymptotics')
plot(hbr_pm(1,:),hbr_pm(2,:),'o','Color',cm(1,:),'DisplayName','Hopf branch')
plot(hbr_pm_pred(1,:),hbr_pm_pred(2,:),'--','Color',cm(7,:),'DisplayName','Hopf asymptotics')
plot(fold_br_pm(1,:),fold_br_pm(2,:),'d','Color',cm(1,:),'DisplayName','Fold branch')
plot(fold_br_pm_pred(1,:),fold_br_pm_pred(2,:),':','Color',cm(7,:),'DisplayName','Fold asymptotics')
plot(bt.parameter(ind.Q),bt.parameter(ind.E),'.k', 'MarkerSize', 18, ...
        'DisplayName','generic Bogdanov-Takens point');
xlabel('$Q$','Interpreter','LaTeX')
ylabel('$E$','Interpreter','LaTeX')
axis([1.2928    1.4323    0.0209    0.0522])
legend

%% Compare predicted and corrected orbits in phase-space
step = linspace(0.1,0.3,10);
[~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0);
[~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
    'codim2','BT','codim1','hcli','step',step,'debug',0);
figure(4); clf; hold on
title('Continued homoclinic orbits in phase-space')
for i=1:length(hcli_br_approx.point)
    profile = hcli_br_approx.point(i).profile;
    plot(profile(1,:),profile(2,:), 'Color', cm(2,:))
    profile = hcli_br_correc.point(i).profile;
    plot(profile(1,:),profile(2,:), '.', 'Color', cm(1,:))
end
legend({'Predicted homoclinic order', 'Corrected homoclinic orbir'})
xlabel('$u_1$','Interpreter','LaTeX');
ylabel('$u_2$','Interpreter','LaTeX');

%% Convergence plot
amplitudes = logspace(-2.6, -1, 20);
orders = [1:3];
relativeerrors = convergence_plot(funcs, bt, orders, amplitudes, 'free_pars', free_pars,...
    'orders',orders,'TTolerance',1e-4,'ntst',82,'debug',false,'dir',1);
figure(5); clf; hold on
title('Convergence plot comparing first and third order homoclinic asymptotics')
plot(log10(amplitudes), log10(relativeerrors{1}(:)), '*-')
plot(log10(amplitudes), log10(relativeerrors{3}(:)), '*-')
legend({'first order', 'thrid order'})
xlabel('amplitude')
ylabel('relative error')
