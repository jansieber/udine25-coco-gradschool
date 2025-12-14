%% Tri-neuron BAM neural network model with delay
% 
% Demo illustrating how to branch off a transcritical Bogdanov-Takens bifurcation point
%

%% Differential equations
% 
% From: Dong, Tao and Liao, Xiaofeng
% Bogdanov-Takens bifurcation in a tri-neuron BAM neural network model with multiple delays
% Nonlinear Dynamics.(2013), no. 3, 583-595.
%
% $$\dot{x}_{1}(t) = -\mu_{1}x_{1}(t)+c_{21}f_{1}(x_{2}(t-\tau_{2}))+c_{31}f_{1}(x_{3}(t-\tau_{2})), $$
% $$\dot{x}_{2}(t) = -\mu_{2}x_{2}(t)+c_{12}f_{2}(x_{1}(t-\tau_{1})),$$
% $$\dot{x}_{3}(t) = -\mu_{3}x_{3}(t)+c_{13}f_{3}(x_{1}(t-\tau_{1})).$$
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
parnames={'alpha1','alpha2','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});

%% Set the funcs structure
funcs=set_symfuncs(@sym_BAMnn_mf,'sys_tau',@()ind.tau);

%% Set bifurcatiun parameter range and step size bounds
brpars={'max_bound', [ind.alpha1 0.5],...
        'min_bound', [ind.alpha1 -0.35],...
        'max_step',  [ind.alpha1 0.01; ind.alpha2 0.01]};

%% Define analytically derived Bogdanov-Takens point
% construct steady-state point
stst=dde_stst_create('x',zeros(3,1));
stst.parameter(ind.alpha1) = 0;
stst.parameter(ind.alpha2) = 0;
stst.parameter(ind.tau) = 5;
% Calculate stability
method=df_mthod(funcs,'stst');
stst.stability=p_stabil(funcs,stst,method.stability);
stst.stability.l1

%% Calculate parameter-dependent normal form coefficients
free_pars=[ind.alpha1, ind.alpha2];
bt = p_tobt(funcs,stst);
bt = nmfm_bt_orbital(funcs, bt,'free_pars', free_pars, 'generic_unfolding', false); 
bt.nmfm

%% Plot comparing profiles computed and predicted homoclinic orbits
figure(1); clf; hold on;
cm=colormap('lines');
tiledlayout(2,6)
ylabels = {'u1','u2','u3'};
step = 0.02;
for order=[1,3]
    [~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
        'debug',0,'codim2','BT','codim1','hcli','step',step,'plot',0, ...
        'generic_unfolding',false,'predictor',true,'order',order);
    [~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars, ...
        'order',order,'codim2','BT','codim1','hcli','step',step, ...
        'plot',0,'generic_unfolding',false,'debug',0);
    for i=1:2
        for j=1:3
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

%% Continue homoclinic curves emanating from Bogdanov-Takens point
[~,hcli_br,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli',brpars{:},'step',[0.003;0.012], ...
    'plot',0,'generic_unfolding',false); assert(all(suc(:)>0))
hcli_br(1)=br_contn(funcs,hcli_br(1),80);
hcli_br(2)=br_contn(funcs,hcli_br(2),80);

%% Continue Hopf curves emanating from Bogdanov-Takens point
[~,hbr,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf',brpars{:},'step',1e-04, ...
    'plot',0,'generic_unfolding',false); assert(all(suc(:)>0))
nop=300; [hbr(1),suc]=br_contn(funcs,hbr(1),nop);assert(suc>0)
nop=80; [hbr(2),suc]=br_contn(funcs,hbr(2),nop);assert(suc>0)

%% Continue transcritical curve emanating from Bogdanov-Takens point
[~,tc_br,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold',brpars{:},'step',1e-03, ...
    'plot',0,'generic_unfolding',false); assert(all(suc(:)>0))
nop=50; [tc_br,suc]=br_contn(funcs,tc_br(1),nop);assert(suc>0)
tc_br = br_rvers(tc_br);
nop=50; [tc_br,suc]=br_contn(funcs,tc_br(1),nop);assert(suc>0)

%% Plot the bifurcation curves
figure(2);clf;hold on
title('Codimension 1 curves emanating from the transcritical Bogdanov-Takens point');
getpars=@(br,ind) arrayfun(@(p)p.parameter(ind),br.point);
hcli_br1_pm = [getpars(hcli_br(1),ind.alpha1), getpars(hcli_br(1),ind.alpha2)]';
hcli_br2_pm = [getpars(hcli_br(2),ind.alpha1), getpars(hcli_br(2),ind.alpha2)]';
hbr1_pm = [getpars(hbr(1),ind.alpha1); getpars(hbr(1),ind.alpha2)];
hbr2_pm = [getpars(hbr(2),ind.alpha1); getpars(hbr(2),ind.alpha2)];
tc_br_pm = [getpars(tc_br(1),ind.alpha1); getpars(tc_br(1),ind.alpha2)];
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'Color',cm(1,:), ...
    'DisplayName','Homoclinic branches')
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'Color',cm(1,:),'HandleVisibility','off')
plot(hbr1_pm(1,:),hbr1_pm(2,:),'Color',cm(2,:),'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'Color',cm(3,:),'DisplayName', ...
    'Transcritical branch')
plot(bt.parameter(ind.alpha1),bt.parameter(ind.alpha2),'.k', ...
    'MarkerSize', 18, 'DisplayName','transcritical Bogdanov-Takens point');
xlabel('$\alpha_1$','Interpreter','LaTeX');
ylabel('$\alpha_2$','Interpreter','LaTeX');
legend

%% Predictors for homoclinic, Hopf, and transcritical curves
step = [linspace(1e-4,0.02,20);linspace(1e-4,0.07,20)]
[~,hcli_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,'debug',0,...
    'codim2','BT','codim1','hcli','step',step,'generic_unfolding',false,'predictor',true);
step = [linspace(0.0001,0.014,20);linspace(0.0001,0.006,20)];
[~,hbr_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf','step',step,'generic_unfolding',false,'predictor',true);
step = linspace(-0.001,0.001,40);
[~,tc_br_pred,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','fold','step',step,'generic_unfolding',false,'predictor',true);
hcli_br1_pm_pred  = [getpars(hcli_br_pred(1), ind.alpha1), ...
    getpars(hcli_br_pred(1), ind.alpha2)];
hcli_br2_pm_pred  = [getpars(hcli_br_pred(2), ind.alpha1), ...
    getpars(hcli_br_pred(2), ind.alpha2)];
hbr1_pm_pred = [getpars(hbr_pred(1),ind.alpha1); getpars(hbr_pred(1),ind.alpha2)];
hbr2_pm_pred = [getpars(hbr_pred(2),ind.alpha1); getpars(hbr_pred(2),ind.alpha2)];
tc_br_pm_pred = [getpars(tc_br_pred,ind.alpha1); getpars(tc_br_pred,ind.alpha2)];

%% Plot the bifurcation curves with predictors
figure(3);clf;hold on
title(['Codimension 1 curves emanating from the transcritical ', ... 
       'Bogdanov-Takens point with predictors']);
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'.','Color',cm(1,:), ...
    'DisplayName','Homoclinic branches')
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'.','Color',cm(1,:), ...
    'HandleVisibility','off')
plot(hbr1_pm(1,:),hbr1_pm(2,:),'o','Color',cm(1,:),'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'o','Color',cm(1,:),'HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'d','Color',cm(1,:), ...
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
plot(bt.parameter(ind.alpha1),bt.parameter(ind.alpha2),'.k', ...
    'MarkerSize', 18, 'DisplayName','generic Bogdanov-Takens point');
legend
xlabel('$\alpha_1$','Interpreter','LaTeX');
ylabel('$\alpha_2$','Interpreter','LaTeX');
axis(1.0e-03*[-0.3696    0.0714   -0.2838    0.2934])

%% Detect bifurcations on Hopf branch II
[hbr2_wbifs,hopftests,bifindx,~]=LocateSpecialPoints(funcs,hbr(2));
bt2 = hbr2_wbifs.point(bifindx(2));
bt2 = p_tobt(funcs,bt2);
bt2 = nmfm_bt_orbital(funcs, bt2,'free_pars',free_pars);
plot(bt2.parameter(1),bt2.parameter(2), '.k', 'MarkerSize', 18, ...
        'DisplayName','generic Bogdanov-Takens point');

%% Plot Bogdanov-Takens test function along Hopf curve
figure(4);clf;hold on
title('Bogdanov--Takens testfunction along Hopf bifurcation curve')
plot3(getpars(hbr2_wbifs,ind.alpha1),getpars(hbr2_wbifs,ind.alpha2),hopftests.bt)
plot3(bt.parameter(ind.alpha1), bt.parameter(ind.alpha2), ...
    hopftests.bt(bifindx(4)), '*k')
plot3(bt2.parameter(ind.alpha1), bbt2.parameter(ind.alpha2)t2.parameter(ind.alpha2), ...
    hopftests.bt(bifindx(2)), '*k')
xlabel('$\alpha_1$','Interpreter','LaTeX');
ylabel('$\alpha_2$','Interpreter','LaTeX');
zlabel('\omega')
grid on
view(-47, 29)

%% Continue homoclinic solutions emanating from the detected BT point
[~,hcli_br(3),suc]=C1branch_from_C2point(funcs,bt2,free_pars,...
    'codim2','BT','codim1','hcli',brpars{:},'step',0.002,'plot',0);
assert(all(suc(:)>0))
hcli_br(3) = br_contn(funcs,hcli_br(3),300);
hcli_br3_pm = [getpars(hcli_br(3),ind.alpha1), getpars(hcli_br(3),ind.alpha2)]';

%% Predictor for homoclinic curve
step = linspace(1e-4,0.015,20);
[~,hcli_br_pred(3),suc]=C1branch_from_C2point(funcs,bt2,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0);
hcli_br3_pm_pred  = [getpars(hcli_br_pred(3), ind.alpha1), ...
    getpars(hcli_br_pred(3), ind.alpha2)];

%% Plot the bifurcation curves with predictors
figure(5);clf;hold on
title(['Codimension 1 curves emanating from the generic ', ... 
       'and transcritical Bogdanov-Takens point with predictors']);
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'.','Color',cm(1,:), ...
    'DisplayName','Homoclinic branches', 'MarkerSize', 20)
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'.','Color',cm(1,:), ...
    'HandleVisibility','off', 'MarkerSize', 20)
plot(hcli_br3_pm(1,:),hcli_br3_pm(2,:),'.','Color',cm(1,:), ...
    'HandleVisibility','off', 'MarkerSize', 20)
plot(hbr1_pm(1,:),hbr1_pm(2,:),'Color','#ccc', ...
    'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'Color','#ccc','HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'--','Color','#ccc', ...
    'DisplayName','Transcritical branch')
plot(hcli_br1_pm_pred(:,1),hcli_br1_pm_pred(:,2),'Color', cm(7,:), ...
    'DisplayName','Third order homoclinic asymptotics')
plot(hcli_br2_pm_pred(:,1),hcli_br2_pm_pred(:,2),'Color', cm(7,:), ...
    'HandleVisibility','off')
plot(hcli_br3_pm_pred(:,1),hcli_br3_pm_pred(:,2),'Color', cm(7,:), ...
    'HandleVisibility','off')
plot(bt.parameter(ind.alpha1),bt.parameter(ind.alpha2),'ok', ...
    'DisplayName','Transcritical Bogdanov-Takens point', ...
    'MarkerFaceColor',[0 0 0], 'MarkerSize', 09);
plot(bt2.parameter(ind.alpha1),bt2.parameter(ind.alpha2),'sk', ...
    'DisplayName','generic Bogdanov-Takens point','MarkerFaceColor', ...
    [0 0 0], 'MarkerSize', 10);
legend
xlabel('$\alpha_1$','Interpreter','LaTeX');
ylabel('$\alpha_2$','Interpreter','LaTeX');
axis([-0.0030, 0.0005, -0.0010, 0.0021])

%% Plot larger bifurcation bifurcation diagram without predictors
figure(6);clf;hold on
tiledlayout(1,2)
hold on;
title(['Codimension 1 curves emanating from the generic ', ... 
       'and transcritical Bogdanov-Takens point with predictors']);
plot(hcli_br1_pm(1,:),hcli_br1_pm(2,:),'Color',cm(1,:), ...
    'DisplayName','Homoclinic branches')
plot(hcli_br2_pm(1,:),hcli_br2_pm(2,:),'Color',cm(1,:),'HandleVisibility','off')
plot(hcli_br3_pm(1,:),hcli_br3_pm(2,:),'Color',cm(1,:),'HandleVisibility','off')
plot(hbr1_pm(1,:),hbr1_pm(2,:),'Color',cm(2,:),'DisplayName','Hopf branches')
plot(hbr2_pm(1,:),hbr2_pm(2,:),'Color',cm(2,:),'HandleVisibility','off')
plot(tc_br_pm(1,:),tc_br_pm(2,:),'Color',cm(3,:), ...
    'DisplayName','Transcritical branch')
plot(bt.parameter(ind.alpha1),bt.parameter(ind.alpha2),'ok', ...
    'DisplayName','Transcritical Bogdanov-Takens point', ...
    'MarkerFaceColor',[0 0 0], 'MarkerSize', 09);
plot(bt2.parameter(ind.alpha1),bt2.parameter(ind.alpha2),'sk', ...
    'DisplayName','generic Bogdanov-Takens point', ...
    'MarkerFaceColor',[0 0 0], 'MarkerSize', 10);
legend
xlabel('$\alpha_1$','Interpreter','LaTeX');
ylabel('$\alpha_2$','Interpreter','LaTeX');

%% Plot continued homoclinic orbits I
figure(7); clf
title('Continued homoclinic orbits in phase-space');
tiledlayout(1,5)
nexttile; hold on;
for i=1:length(hcli_br(1).point)
    profile = hcli_br(1).point(i).profile;
    plot3(profile(1,:),profile(2,:),profile(3,:), 'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot continued homoclinic orbits I (rotated)
nexttile; hold on;
R = @(alpha) [cos(alpha) -sin(alpha); sin(alpha) cos(alpha)]
for i=1:length(hcli_br(1).point)
    profile = hcli_br(1).point(i).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(profileRotated(1,:),profileRotated(2,:),profile(3,:), 'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot continued homoclinic orbits II
nexttile; hold on;
for i=1:length(hcli_br(2).point)
    profile = hcli_br(2).point(i).profile;
    plot3(profile(1,:),profile(2,:),profile(3,:), 'Color',cm(1,:))
end
profile = hcli_br(2).point(end).profile;
plot3(profile(1,:),profile(2,:),profile(3,:), 'Color',cm(2,:))
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot continued homoclinic orbits III
nexttile; hold on;
for i=1:length(hcli_br(3).point)
    profile = hcli_br(3).point(i).profile;
    plot3(profile(1,:),profile(2,:),profile(3,:), 'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot continued homoclinic orbits III (rotated)
nexttile; hold on;
for i=1:length(hcli_br(3).point)
    profile = hcli_br(3).point(i).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(profileRotated(1,:),profileRotated(2,:),profile(3,:), 'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot continued homoclinic orbits in (alpha_1,t,u_1) space
figure(8); clf; hold on;
title(''LaTeX' continued homoclinic orbits in (alpha_1,t,u_1) space')
for i=1:length(hcli_br(1).point)
    profile = hcli_br(1).point(i).profile;
    plot3(hcli_br(1).point(i).parameter(ind.alpha1)*ones(size(profile(1,:))), ....
        hcli_br(1).point(i).mesh,profile(1,:), 'Color',cm(1,:))
end
for i=1:length(hcli_br_pred(1).point)
    profile = hcli_br_pred(1).point(i).profile;
    plot3(hcli_br_pred(1).point(i).parameter(ind.alpha1)*ones(size(profile(1,:))), ...
        hcli_br_pred(1).point(i).mesh,profile(1,:), '.', 'Color',cm(2,:))
end
grid on
view(-12,12)
xlabel('$\alpha_1$','Interpreter','LaTeX')
ylabel('$t$','Interpreter','LaTeX')
zlabel('$u_1$','Interpreter','LaTeX')

%% Compare homoclinic solutions in phase-space
step = linspace(0.003,0.009,10);
[~,hcli_br_approx,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'predictor',true,'debug',0, ...
    'generic_unfolding', false);
[~,hcli_br_correc,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hcli','step',step,'debug',0, ...
    'generic_unfolding', false);
figure(9); clf; hold on
title('Compare homoclinic orbits in phase-space')
for i=1:length(hcli_br_approx(1).point)
    for j=1:2
        point = hcli_br_approx(j).point(i);
        profile = point.profile;
        mesh = point.mesh;
        plot3(point.parameter(ind.alpha2)*ones(size(profile(1,:))), ...
            mesh,profile(2,:), 'Color',cm(2,:))
        point_correc = hcli_br_correc(j).point(i);
        profile_correc = point.profile;
        mesh_correc = point_correc.mesh;
        plot3(point_correc.parameter(ind.alpha2)*ones(size(profile(1,:))), ...
            mesh,profile_correc(2,:), '.', 'Color',cm(1,:))
    end
end
xlabel('$\epsilon$','Interpreter','LaTeX')
ylabel('$\tilde t$','Interpreter','LaTeX')
zlabel('$u_1$','Interpreter','LaTeX')
grid on
view(-46,7)
legend({'Predicted homoclinic order', 'Corrected homoclinic orbir'})

%% Continue periodic orbit from Hopf point on the first Hopf branch
psol_br1=SetupPsol(funcs,hbr(1),length(hbr(1).point),'contpar',ind.alpha2, ... 
    'degree',3,'intervals',50,brpars{1:4},'max_step',[0,inf],'plot',0);
psol_br1=br_contn(funcs,psol_br1,40);
psol_br1_pm = [getpars(psol_br1,ind.alpha1); getpars(psol_br1,ind.alpha2)]';

%% Add parameters of periodic branch to bifurcation diagram
figure(6);
plot(psol_br1_pm(:,1),psol_br1_pm(:,2),'.','Color',cm(3,:), ...
    'HandleVisibility','off')

%% Plot periodic branch
% we see that the perioidic orbit converges to a homoclinic orbit
figure(10); clf; hold on;
title('Continued periodic orbits in phase-space');
for i=1:length(psol_br1.point)
    profile = psol_br1.point(i).profile;
    plot3(profile(1,:),profile(2,:),profile(3,:), 'Color',cm(1,:))
end
profile = psol_br1.point(end).profile;
plot3(profile(1,:),profile(2,:),profile(3,:),'.', 'Color',cm(2,:))
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Set bifurcatiun parameter range and step size bounds
brpars={'max_bound', [ind.alpha1 0.5],...
        'min_bound', [ind.alpha1 -0.35],...
        'max_step',  [ind.alpha1 1e-5; ind.alpha2 1e-5]};

%% Continue Hopf curve emanating from Bogdanov-Takens point again
% but with smaller stepsize
[~,hbr_small_stepsize,suc]=C1branch_from_C2point(funcs,bt,free_pars,...
    'codim2','BT','codim1','hopf',brpars{:},'step',1e-04, ...
    'plot',0,'generic_unfolding',false); assert(all(suc(:)>0))
nop=285; [hbr_small_stepsize(2),suc]=br_contn(funcs,hbr_small_stepsize(2),nop);
assert(suc>0)
hbr_small_stepsize_pm = [getpars(hbr_small_stepsize(2),ind.alpha1); 
    getpars(hbr_small_stepsize(2),ind.alpha2)];

%% Add parameters of new Hopf curve to bifurcation diagram
figure(6);
plot(hbr_small_stepsize_pm(1,:),hbr_small_stepsize_pm(2,:),'*','Color',cm(4,:), ...
    'HandleVisibility','off')

%% Continue periodic orbits from every Hopf point on the second Hopf branch
% (takes a few minutes)
start_ind = 19;
end_ind = length(hbr_small_stepsize_pm)-5;
for i=start_ind:end_ind
    psol_br(i-start_ind+1)=SetupPsol(funcs,hbr_small_stepsize(2),i,'contpar',...
        ind.alpha2, ... 
        'degree',3,'intervals',50,brpars{1:4},'max_step',[0,inf],'plot',0, ...
        'radius', 0.001);
    psol_br(i-start_ind+1)=br_contn(funcs,psol_br(i-start_ind+1),50);
end

%% Plot last points of every psol_br
figure(11); clf; hold on
for i=1:end_ind-start_ind+1
    profile = psol_br(i).point(end).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(profileRotated(1,:),profileRotated(2,:),profile(3,:),'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$u_1$','Interpreter','LaTeX')
ylabel('$u_2$','Interpreter','LaTeX')
zlabel('$u_3$','Interpreter','LaTeX')

%% Plot homoclinic solutions in (\alpha_2,\tilde u_1, \tilde u_2) space
figure(12); clf; hold on
for i=1:end_ind-start_ind+1
    profile = psol_br(i).point(end).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(psol_br(i).point(end).parameter(ind.alpha1)*ones(size(profile(1,:))), ...
        profileRotated(1,:),profileRotated(2,:),'Color',cm(2,:))
end
for i=1:length(hcli_br(1).point)
    profile = hcli_br(1).point(i).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(hcli_br(1).point(i).parameter(ind.alpha1)*ones(size(profile(1,:))), ...
        profileRotated(1,:),profileRotated(2,:), 'Color',cm(1,:))
end
for i=1:length(hcli_br(3).point)
    profile = hcli_br(3).point(i).profile;
    profileRotated = R(0.29)*profile([1,2],:)
    plot3(hcli_br(3).point(i).parameter(ind.alpha1)*ones(size(profile(1,:))), ...
        profileRotated(1,:),profileRotated(2,:), 'Color',cm(1,:))
end
grid on
view(-12,12)
xlabel('$\alpha_1$','Interpreter','LaTeX')
ylabel('$\tilde u_1$','Interpreter','LaTeX')
zlabel('$\tilde u_2$','Interpreter','LaTeX')
plot3(bt.parameter(ind.alpha1),bt.x(1),bt.x(2),'.k', ...
    'MarkerSize', 18, 'DisplayName','transcritical Bogdanov-Takens point');
bt2rotated = R(0.29)*bt2.x([1,2])
plot3(bt2.parameter(ind.alpha1),bt2rotated(1),bt2rotated(2),'.k', ...
    'MarkerSize', 18, 'DisplayName','Bogdanov-Takens point');

%% Extract homoclinic parameters from the psol branches
hcli_bifurcation_br_pm = []
for i=1:end_ind-start_ind+1
    parameter = psol_br(i).point(end).parameter;
    hcli_bifurcation_br_pm = [hcli_bifurcation_br_pm; parameter(ind.alpha1:ind.alpha2)];
end
figure(6)
plot(hcli_bifurcation_br_pm(:,1), hcli_bifurcation_br_pm(:,2), 'Color', cm(1,:))

%% Convergence plot
amplitudes = logspace(-3, -1.4, 20);
orders = [1:3];
relativeErrors = convergence_plot(funcs, bt, orders, amplitudes, 'free_pars', free_pars,...
    'orders',orders,'TTolerance',1e-5,'ntst',82,'generic_unfolding',0,'debug',false);
figure(13); clf; hold on
title('Convergence plot comparing first and third order homoclinic asymptotics')
plot(log10(amplitudes), log10(relativeErrors{1}(:)), '*-')
plot(log10(amplitudes), log10(relativeErrors{3}(:)), '*-')
legend({'first order', 'thrid order'})
xlabel('amplitude')
ylabel('relative error')
