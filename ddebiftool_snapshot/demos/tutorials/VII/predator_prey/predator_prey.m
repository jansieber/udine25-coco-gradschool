%% Predator-prey system with double Allee effect 
% 
% Demo illustrating how to branch off a genric Bogdanov-Takens bifurcation point
%

%% Differential equations of a delayed predator-prey system with double Allee effect
% 
% From: Jianfeng Jiao and Can Chen
% Bogdanov-Takens bifurcation analysis of a delayed predator-prey system with 
% double Allee effect
% Nonlinear Dynamics, 104 (2021), pp. 1697–1707.
%
% $$\dot x(t) = \dfrac{rx}{x+n_0}\left(1-\dfrac1 K\right)\left(x - m_0\right) - \dfrac{cxy}{x+\varrho y},$$
%
% $$y'=\\dot y(t) = -dy + \dfrac{c_1 x(t-\tau)y}{x(t-\tau) + \varrho y(t-\tau)}.$$

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
parnames={'theta','delta','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});

%% Set the funcs structure
% We load the precalculated multilinear forms. These have been generated
% with the file gen_sym_predator_prey.m.
funcs=set_symfuncs(@sym_predator_prey_mf,'sys_tau',@()ind.tau);

%% Set bifurcation parameter range and step size bounds
brpars={'max_bound', [ind.delta 0.7],...
        'min_bound', [ind.delta 0.4],...
        'max_step', [ind.theta 0.002; ind.delta 0.002]};

%% Define constants
m = 1.502983803; alpha = 0.9; gamma = 0.15;

%% Define analytically derived Bogdanov-Takens point
% construct steady-state point
x0 = (m*gamma-m*alpha+m+alpha)/(2*m);
stst=dde_stst_create('x',[x0;  (m-1)*x0]);
stst.parameter(ind.theta)  = (m^2*gamma^2 - 2*m*(m*alpha + m - alpha)*gamma + (m*alpha - m - alpha)^2)/(4*m*alpha*(m-1));
stst.parameter(ind.delta) = alpha/m;
stst.parameter(ind.tau)   = 0.01;
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1

%% Calculate normal form coefficients
free_pars=[ind.theta, ind.delta];
bt = p_tobt(funcs,stst);
bt = nmfm_bt_orbital(funcs, bt,'free_pars', free_pars);

%% Plot comparing profiles computed and predicted homoclinic orbits
figure(1); clf; hold on;
cm=colormap('lines');
tiledlayout(2,2)
ylabels = {'$x$','$y$'};
step = 0.3;
for order=[1,3]
    [~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
        'codim2','BT','codim1','hcli','step',step,'predictor',true);
    [~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,'order',order,...
        'codim2','BT','codim1','hcli','step',step);
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
    'codim2','BT','codim1','hcli',brpars{:},'step',[0.1;0.11],'plot',0);
assert(all(suc(:)>0))
nop=300; hcli_br=br_contn(funcs,hcli_br,nop);assert(suc>0)

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
title('Codimension 1 curves emanating from the generic Bogdanov-Takens point');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hcli_br_pm = [getpars(hcli_br,ind.theta), getpars(hcli_br,ind.delta)]';
hbr_pm = [getpars(hbr,ind.theta); getpars(hbr,ind.delta)];
fold_br_pm = [getpars(fold_br,ind.theta); getpars(fold_br,ind.delta)];
plot(hcli_br_pm(1,:),hcli_br_pm(2,:),'Color',cm(1,:),'DisplayName','Homoclinic branches')
plot(hbr_pm(1,:),hbr_pm(2,:),'Color',cm(2,:),'DisplayName','Hopf branches')
plot(fold_br_pm(1,:),fold_br_pm(2,:),'Color',cm(3,:),'DisplayName','Fold branch')
h(1) = plot(bt.parameter(ind.theta),bt.parameter(ind.delta),'.k', 'MarkerSize', 18, ...
        'DisplayName','generic Bogdanov-Takens point');
xlabel('$\theta$','Interpreter','LaTeX');
ylabel('$\delta$','Interpreter','LaTeX');
legend

%% Predictors for homoclinic, Hopf, and fold curve
step = linspace(1e-4,0.6,50);
[~,hcli_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,'debug',0,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true);
step = linspace(0.0001,0.1,50);
[~,hbr_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf','step',step,'predictor',true);
step = linspace(-0.05,0.05,40);
[~,fold_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold','step',step,'predictor',true);
hcli_br_pred_pm  = [getpars(hcli_br_pred, ind.theta), getpars(hcli_br_pred, ind.delta)];
hbr_pm_pred = [getpars(hbr_pred,ind.theta); getpars(hbr_pred,ind.delta)];
fold_br_pm_pred = [getpars(fold_br_pred,ind.theta); getpars(fold_br_pred,ind.delta)];

%% Predictors from [Jiao2021]
mu1 = -linspace(1e-5,0.003,250);
mu2 = -1.386173725.*mu1 - 0.9901240893*sqrt(mu1.*(mu1 - 3.652165730));
hbr_pm_pred_2021 = bt.parameter(free_pars)' + [mu1; mu2];
mu2 = -1.386173725.*mu1 - 1.386173725*sqrt(mu1.*(mu1 - 3.652165730));
hcli_pm_pred_2021 = bt.parameter(free_pars)' + [mu1; mu2];

%% Plot the bifurcation curves with predictors
figure(3);clf;hold on
title('Codimension 1 curves emanating from the generic Bogdanov-Takens point with predictors');
plot(hcli_br_pm(1,:),hcli_br_pm(2,:),'.','Color',cm(1,:),'DisplayName','Homoclinic branches')
plot(hcli_br_pred_pm(:,1),hcli_br_pred_pm(:,2),'Color', cm(7,:),'DisplayName','Third order homoclinic asymptotics')
plot(hcli_pm_pred_2021(1,:), hcli_pm_pred_2021(2,:), 'Color',cm(3,:),'DisplayName','Asymptotics from Jiao2021');
plot(hbr_pm(1,:),hbr_pm(2,:),'o','Color',cm(1,:),'DisplayName','Hopf branch')
plot(hbr_pm_pred(1,:),hbr_pm_pred(2,:),'--','Color',cm(7,:),'DisplayName','Hopf asymptotics')
plot(hbr_pm_pred_2021(1,:),hbr_pm_pred_2021(2,:),'--','Color',cm(3,:),'DisplayName','Hopf asymptotics from Jiao2021')
plot(fold_br_pm(1,:),fold_br_pm(2,:),'d','Color',cm(1,:),'DisplayName','Fold branch')
plot(fold_br_pm_pred(1,:),fold_br_pm_pred(2,:),':','Color',cm(7,:),'DisplayName','Fold asymptotics')
plot(bt.parameter(ind.theta),bt.parameter(ind.delta),'.k', 'MarkerSize', 18, ...
        'DisplayName','generic Bogdanov-Takens point');
legend('Location','NorthWest')
xlabel('$\theta$','Interpreter','LaTeX');
ylabel('$\delta$','Interpreter','LaTeX');
axis([0.0978    0.1002    0.5086    0.6042])

%% Compare homoclinic solutions in phase-space
step = linspace(0.1,0.3,10);
[~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0);
[~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'debug',0);
figure(4); clf; hold on
title('Compare homoclinic orbits in phase-space')
for i=1:length(hcli_br_approx.point)
    profile = hcli_br_approx.point(i).profile;
    plot(profile(1,:),profile(2,:), 'Color', cm(2,:))
    profile_correc = hcli_br_correc.point(i).profile;
    plot(profile_correc(1,:),profile_correc(2,:), '.', 'Color', cm(1,:))
end
legend({'Predicted homoclinic order', 'Corrected homoclinic orbir'})
xlabel('$x$','Interpreter','LaTeX');
ylabel('$y$','Interpreter','LaTeX');


%% Convergence plot
amplitudes = logspace(-4, -1, 20);
orders = [1:3];
relativeerrors = convergence_plot(funcs, bt, orders, amplitudes, 'free_pars', free_pars,...
    'orders',orders,'ttolerance',1e-6,'ntst',82,'debug',false,'dir',1);
figure(5); clf; hold on
title('Convergence plot comparing first and third order homoclinic asymptotics')
plot(log10(amplitudes), log10(relativeerrors{1}(:)), '*-')
plot(log10(amplitudes), log10(relativeerrors{3}(:)), '*-')
legend({'first order', 'thrid order'})
xlabel('amplitude')
ylabel('relative error')
