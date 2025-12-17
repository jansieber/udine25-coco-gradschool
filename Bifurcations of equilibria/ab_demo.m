%% Demo AB reaction
% x(1)=concentration fraction
% x(2)=temperature
% parameters [alpha,beta,gamma]=reaction rate, exothermicity, cooling
%% set path
clear;
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
%% Define parameters and variables
pnames={'alpha','beta','gamma'};
xnames={'c','T'}; % concentration and temperature
[iv,ip]=deal(structind_from_names(xnames),structind_from_names(pnames));
%% define right-hand side
ab=@(c,T,alpha,beta,gamma)[...     % ODE r.h.s.
    -c+alpha.*(1-c).*exp(T);...
    -T+beta.*alpha.*(1-c).*exp(T)-gamma.*T];
rhs=@(x,p)ab(x(iv.c,:),x(iv.T,:),p(ip.alpha,:),p(ip.beta,:),p(ip.gamma,:));
%% initial guess for equilibrium
par0([ip.alpha,ip.beta,ip.gamma])=...
     [     0;       14;     2];
par0=par0(:);
x0=[0;0];
%% primary bifurcation alpha, suggested range [0,0.15]
%% secondary continuation parameter beta [0,30]
prob=coco_prob();
prob=ode_isol2ep(prob,'',rhs,x0,pnames,par0); 
coco(prob,'ab_ep',[],'alpha', [0,0.15]); % branch label, [], branch dimension, active continuation parameter, computational domain
% visualize the result
figure(1)
clf
theme = struct('special', {{'SN', 'HB'}}); % plotting theme (check with ep_plot_theme())
coco_plot_bd(theme, 'ab_ep', 'alpha', 'x') % 'x' denotes the state vector, which is scalar here
axis tight
grid on
%%
SN = coco_bd_labs('ab_ep', 'SN');
prob=coco_prob();
prob=ode_ep2SN(prob,'','ab_ep', SN(1));
coco(prob,'ab_SN',[], {'alpha' 'beta'}, {[0,0.15] [0,30]});
bd_abSN=coco_bd_table('ab_SN','numlab',true) %#ok<*NOPTS>
%%
figure(2);clf;
plot(bd_abSN.alpha,bd_abSN.beta,'-','LineWidth',2)
%%
%%
HB = coco_bd_labs('ab_ep', 'HB');
prob=coco_prob();
prob=ode_ep2HB(prob,'','ab_ep', HB);
prob=coco_set(prob,'cont','PtMX',[-200,200])
coco(prob,'ab_HB',[], {'alpha' 'beta'}, {[0,0.15] [0,30]});
bd_abHB=coco_bd_table('ab_HB','numlab',true) %#ok<*NOPTS>
%%
figure(2);hold on;
plot(bd_abHB.alpha,bd_abHB.beta,'-','LineWidth',2)
