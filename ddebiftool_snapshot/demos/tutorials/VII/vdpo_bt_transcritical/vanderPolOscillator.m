%% Van der Pol Oscillator with delay
% 
% Demo illustrating how to branch off a transcritical Bogdanov-Takens bifurcation point
%

%% Differential equations
% 
% From: W. Jiang and Y. Yuan
% Bogdanov–Takens singularity in van der Pol’s oscillator with delayed feedback
% neural system with delay.
% Physica D: Nonlinear Phenomena, 227 (2007), pp. 149–161.
%
% $$\dot ẋ_1 = x_2,$$
%
% $$\dot x_2 = \epsilon g(x_1(t-\tau))-\epsilon(x_1^2-1)x_2-x_1,$$
%
% where
%
% g(x) = \frac{e^x-1}{c_1e^x + c_2}.

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
parnames={'epsilon','tau','vartau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});

%% Set the funcs structure
funcs=set_symfuncs(@sym_vdpo_mf,'sys_tau',@()ind.vartau);

%% Set bifurcatiun parameter range and step size bounds
brpars={'min_bound', [ind.epsilon 0.7122; ind.tau 0.7154],...
        'max_bound', [ind.epsilon 0.7998; ind.tau 0.7905],...
        'max_step',  [ind.epsilon 0.0001; ind.tau 0.0001]};

%% Define analytically derived Bogdanov-Takens point
% manually construct steady-state point
epsilon = 0.75; tau = 0.75;
stst=dde_stst_create('x',[0;0]);
stst.parameter(ind.epsilon)  = epsilon;
stst.parameter(ind.tau) = tau;
stst.parameter(ind.vartau) = 1;

%% Compute stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1

%% Calculate parameter-dependent normal form coefficients
free_pars=[ind.epsilon, ind.tau];
bt = p_tobt(funcs,stst);
bt = nmfm_bt_orbital(funcs, bt,'free_pars', free_pars, 'generic_unfolding', false);
bt.nmfm

%% Plot comparing profiles computed and predicted homoclinic orbits
figure(1); clf; hold on;
cm=colormap('lines');
tiledlayout(2,4)
ylabels = {'$x$','$y$'};
step = 0.1;
for order=[1,3]
    [~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
        'debug',0,'codim2','BT','codim1','hcli','step',step,'plot',0, ...
        'generic_unfolding',false,'predictor',true,'order',order);
    [~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
        'order',order,'codim2','BT','codim1','hcli','step',step, ...
        'plot',0,'generic_unfolding',false,'debug',0);
    for i=1:2
        for j=1:2
            ax = nexttile; hold on; title(sprintf('order %d',order))
            plot(ax,hcli_br_approx(i).point(1).mesh, ...
                hcli_br_approx(i).point(1).profile(j,:),'*')
            plot(ax,hcli_br_correc(i).point(1).mesh, ...
                hcli_br_correc(i).point(1).profile(j,:))
            xlabel('$t$','Interpreter','LaTeX'); 
            ylabel(string(ylabels(j)),'Interpreter','LaTeX')
        end
    end
end

%% Continue homoclinic curve emanating from Bogdanov-Takens point
[~,hcli_br,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
    'codim2','BT','codim1','hcli',brpars{:},'step',0.02, ...
    'plot',0,'generic_unfolding',false);
assert(all(suc(:)>0))
nop=300;[hcli_br(1),suc]=br_contn(funcs,hcli_br(1),nop);assert(suc>0)
nop=300;[hcli_br(2),suc]=br_contn(funcs,hcli_br(2),nop);assert(suc>0)

%% Continue Hopf curve emanating from Bogdanov-Takens point
[~,hbr,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
    'codim2','BT','codim1','hopf',brpars{:},'step',1e-05, ...
    'plot',0,'generic_unfolding',false);
assert(all(suc(:)>0))
nop=300; [hbr(1),suc]=br_contn(funcs,hbr(1),nop); assert(suc>0)
nop=300; [hbr(2),suc]=br_contn(funcs,hbr(2),nop); assert(suc>0)
[hbr(2),suc]=br_contn(funcs,hbr(2),nop); assert(suc>0)

%% Continue transcritical curve emanating from Bogdanov-Takens point
[~,tc_br,suc]=C1branch_from_C2point(funcs,bt,free_pars,'dir',2, ...
    'codim2','BT','codim1','fold',brpars{:},'step',-1e-4, ...
    'plot',0,'generic_unfolding',false);
assert(all(suc(:)>0))
nop=300; [tc_br,suc]=br_contn(funcs,tc_br(1),nop);assert(suc>0)
tc_br = br_rvers(tc_br);
nop=300; [tc_br,suc]=br_contn(funcs,tc_br(1),nop);assert(suc>0)

%% Plot the bifurcation curves
figure(2);clf;hold on
title('Codimension 1 curves emanating from the transcritical Bogdanov-Takens point');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hcli_br1_pm = [getpars(hcli_br(1),ind.epsilon), getpars(hcli_br(1),ind.tau)]';
hcli_br2_pm = [getpars(hcli_br(2),ind.epsilon), getpars(hcli_br(2),ind.tau)]';
hbr1_pm = [getpars(hbr(1),ind.epsilon); getpars(hbr(1),ind.tau)];
hbr2_pm = [getpars(hbr(2),ind.epsilon); getpars(hbr(2),ind.tau)];
tc_br_pm = [getpars(tc_br(1),ind.epsilon); getpars(tc_br(1),ind.tau)];
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'Color',cm(1,:),'DisplayName','Homoclinic branches')
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'Color',cm(1,:),'HandleVisibility','off')
plot(hbr1_pm(1,:),hbr1_pm(2,:),'Color',cm(2,:),'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'Color',cm(3,:),'DisplayName','Transcritical branch')
h(1) = plot(bt.parameter(ind.epsilon),bt.parameter(ind.tau),'.k', 'MarkerSize', 18, ...
        'DisplayName','transcritical Bogdanov-Takens point');
xlabel('$\epsilon$','Interpreter','LaTeX');
ylabel('$\tau$','Interpreter','LaTeX');
legend

%% Predictors for homoclinic, Hopf, and transcritical curves
step = [linspace(1e-4,0.1,20);linspace(1e-4,0.1,20)];
[~,hcli_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,'debug',0, ...
    'codim2','BT','codim1','hcli','step',step,'generic_unfolding',false, ...
    'predictor',true);
step = [linspace(0.0001,0.1,20);linspace(0.0001,0.1,20)];
[~,hbr_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf','step',step,'generic_unfolding',false, ...
    'predictor',true);
step = linspace(-0.04,0.04,40);
[~,tc_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold','step',step,'generic_unfolding',false, ...
    'predictor',true);
hcli_br1_pm_pred  = [getpars(hcli_br_pred(1), ind.epsilon), ...
    getpars(hcli_br_pred(1), ind.tau)];
hcli_br2_pm_pred  = [getpars(hcli_br_pred(2), ind.epsilon), ...
    getpars(hcli_br_pred(2), ind.tau)];
hbr1_pm_pred = [getpars(hbr_pred(1),ind.epsilon); getpars(hbr_pred(1),ind.tau)];
hbr2_pm_pred = [getpars(hbr_pred(2),ind.epsilon); getpars(hbr_pred(2),ind.tau)];
tc_br_pm_pred = [getpars(tc_br_pred,ind.epsilon); getpars(tc_br_pred,ind.tau)];

%% Plot the bifurcation curves with predictors
figure(3);clf;hold on
title(['Codimension 1 curves emanating from the transcritical', ...
    'Bogdanov-Takens point with predictors']);
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'.','Color',cm(1,:), ...
    'DisplayName','Homoclinic branches')
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'.','Color',cm(1,:), ...
    'HandleVisibility','off')
plot(hbr1_pm(1,:),hbr1_pm(2,:),'o','Color',cm(2,:),'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'o','Color',cm(2,:),'HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'d','Color',cm(3,:), ...
    'DisplayName','Transcritical branch')
plot(hcli_br1_pm_pred(:,1),hcli_br1_pm_pred(:,2),'Color', cm(7,:), ...
    'DisplayName','Third order homoclinic asymptotics')
plot(hcli_br2_pm_pred(:,1),hcli_br2_pm_pred(:,2),'Color', cm(7,:), ...
    'HandleVisibility','off')
plot(hbr1_pm_pred(1,:),hbr1_pm_pred(2,:),'--','Color',cm(7,:), ...
    'DisplayName','Hopf asymptotics')
plot(hbr2_pm_pred(1,:),hbr2_pm_pred(2,:),'--','Color',cm(7,:), ...
    'HandleVisibility','off')
plot(tc_br_pm_pred(1,:),tc_br_pm_pred(2,:),':','Color',cm(7,:), ...
    'DisplayName','Transcritical asymptotics')
plot(bt.parameter(ind.epsilon),bt.parameter(ind.tau),'.k', 'MarkerSize', 18, ...
    'DisplayName','generic Bogdanov-Takens point');
legend
xlabel('$\epsilon$','Interpreter','LaTeX');
ylabel('$\tau$','Interpreter','LaTeX');

%% Plot continued homoclinic orbits
figure(4); clf; hold on;
title('Continued homoclinic orbits in (epsilon,x1,x2)-space');
for i=1:length(hcli_br(1).point)
       profile = hcli_br(1).point(i).profile;
       plot3(hcli_br(1).point(i).parameter(ind.epsilon)*ones(size(profile(1,:))), ...
           profile(1,:),profile(2,:), 'Color',cm(1,:))
end
for i=1:length(hcli_br(2).point)
       profile = hcli_br(2).point(i).profile;
       plot3(hcli_br(2).point(i).parameter(ind.epsilon)*ones(size(profile(1,:))), ...
           profile(1,:),profile(2,:), 'Color',cm(1,:))
end
xlabel('$\epsilon$','Interpreter','LaTeX')
ylabel('$x_1$','Interpreter','LaTeX')
zlabel('$x_2$','Interpreter','LaTeX')
grid on
view(-46,7)

%% Compare homoclinic solutions in phase-space
step = linspace(0.01,0.03,10);
[~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0, ...
    'generic_unfolding', false);
[~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'debug',0, ...
    'generic_unfolding', false);
figure(5); clf; hold on
title('Compare homoclinic orbits in phase-space')
for i=1:length(hcli_br_approx(1).point)
    for j=1:2
        point = hcli_br_approx(j).point(i);
        profile = point.profile;
        plot3(point.parameter(ind.epsilon)*ones(size(profile(1,:))), ...
            profile(1,:),profile(2,:), 'Color',cm(2,:))
        point_correc = hcli_br_correc(j).point(i);
        profile_correc = point.profile;
        plot3(point_correc.parameter(ind.epsilon)*ones(size(profile(1,:))), ...
            profile_correc(1,:),profile_correc(2,:), '.', 'Color',cm(1,:))
    end
end
xlabel('$\epsilon$','Interpreter','LaTeX')
ylabel('$x_1$','Interpreter','LaTeX')
zlabel('$x_2$','Interpreter','LaTeX')
grid on
view(-46,7)
legend({'Predicted homoclinic order', 'Corrected homoclinic orbir'})

%% Convergence plot
amplitudes = logspace(-3.4, -1, 20);
orders = [1:3];
relativeErrors = convergence_plot(funcs, bt, orders, amplitudes, 'free_pars', free_pars,...
    'orders',orders,'TTolerance',9e-6,'generic_unfolding',0,'debug',false);
figure(6); clf; hold on
title('Convergence plot comparing first and third order homoclinic asymptotics')
plot(log10(amplitudes), log10(relativeErrors{1}(:)), '*-')
plot(log10(amplitudes), log10(relativeErrors{3}(:)), '*-')
legend({'first order', 'thrid order'})
xlabel('amplitude')
ylabel('relative error')
