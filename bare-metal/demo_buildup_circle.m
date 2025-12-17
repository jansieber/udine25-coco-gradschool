%% Illustration of use for coco_add_func for gradually building up continuation problems
%
% For more extensively commented demos covering the material below, see 
% 
% * |recipes/cmds_demo/demo_basic.m|
% * |recipes/cmds_demo/demo_stages.m| 
% * |recipes/cmds_demo/demo_events.m| 
% 
%% Perform startup of coco somewhere to set paths
% coco has startup.m file in its installation folder. Executing it adds the
% relevant folders to the path.
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'));
addpath(fullfile(pwd(),'tools'))
screen=@(varargin)cellfun(@(x)disp(x),varargin); % hack for disp with description
%% Implicit definitions of a circle - Zero problem Phi1 and initial guess
Phi1=@circle;
u0=[1; 1.1];
%% Create empty coco problem and add Phi zero problem and initial guess for u0
% The guess also determines input dimension for Phi1.
% Inputs for coco_add_func: problem, fcn-name, function handle, (initial)
% data, flag for type of function {'zero','regular','singular','a}
prob=coco_prob();
prob=coco_add_func(prob,'Phi1',Phi1,[],'zero','u0',u0); % add 'f+df' if Phi1 provides Jacobian 
% set parameters for run (atlas, corrector, toolbox options)
prob=coco_set(prob,'cont','PtMX',100);
%% continuation run
% arguments: problem, run name, toolbox selector (usually empty), dimension
% of manifold ...
bd=coco(prob, 'run1', [], 1);
%% We could add another zero function
% This zero function depends o,n the second variable of the already defined
% u. This is indicated by the uidx input pair. Now we have 2 eqs, 2 vars,
% so solution manifold has dim 0.
Phi2=@(u)u-1; % Phi_2=u_2-1
JPhi2=@(u)1;
prob2=coco_add_func(prob,'Phi2',f2coco(Phi2),f2coco(JPhi2),[],'zero','uidx',2);
bd2=coco(prob2, 'run2', [],0);
%% We may have dependence on previous variables and on new variables in new zero functions
% overall system for run2a: 0=u1^2+(u2-1)^2-1, u2+u3=0
Phi2a=@(u)u(1)+u(2);
prob2a=coco_add_func(prob,'Phi2a',f2coco(Phi2a),[],'zero','uidx',2,'u0',-1);
bd2a=coco(prob2a, 'run2a', [], 1);
% 
%% Add monitor function Psi1
% (the squared norm of u) and call it continuation parameter p
% 'active' means 'can vary during run'
%%
Psi1=@(u)u(1)^2+u(2)^2;
prob_p=coco_add_func(prob2a,'Psi1',f2coco(Psi1),[],'active','p','uidx',[1,2]);
bd_p=coco(prob_p,'p_demo',[],1,{'p'},{[1,5]});
%% The above is bad
% One should not rely on the ordering of of variables, but extract them
% from the problem:
[uidx_Phi1,data_Phi1]=coco_get_func_data(prob2a,'Phi1','uidx','data');
prob_p=coco_add_func(prob2a,'Psi1',f2coco(Psi1),[],'active','p','uidx',uidx_Phi1);
bd_p=coco(prob_p,'p_demo',[],1,{'p'},{[1,5]});
%% Check dependence of Phi2a
uidx_Phi2a=coco_get_func_data(prob2a,'Phi2a','uidx')
%% Equate some variables with continuation parameters
% or simply introduce names for variables
% continuation parameters have names and (optional) boundaries
prob_pxy=coco_add_pars(prob_p,'xy',uidx_Phi1,{'x','y'});
bd_pxy=coco(prob_pxy,'pxy_demo',[],1,{'x','y','p'},{[0,2]}); % swap x and y position to see effect on 'FP' outputs
%% Check functions present, their dependence on variables and data
pxy_info=prob_fcn_info(prob_pxy)
%% Read solution from disk
% file is, eg, data/p_demo/sol1.mat
%%
chart_p = coco_read_solution('p_demo', 1,'chart'); %Extract solution from disk
screen('solution in LAB1 of p_demo: chart_p.x',chart_p.x); % solution of dimension-0 manifold
chartphi2a = coco_read_solution('Phi2a','run2a', 1,'chart'); %Extract solution from disk
screen('u for Phi1a in LAB1 in 1dim mf: chartphi2a.x',chartphi2a.x); % solution of dimension-1 manifold
chartrun2a = coco_read_solution('run2a', 1,'chart');
screen('complete solution in LAB1 in 1dim mf: chartrun2a.x',chartrun2a.x); % solution of dimension-1 manifold
screen('complete tangent in LAB1 in 1dim mf: chartrun2a.t',chartrun2a.t); % tangent at x for dimension-1 manifold
%% Other useful things
%% Add events: in the simplest case, when a continuation parameter reaches a value
% also set output properties: switch off printing of every 10th point
prob_pxy=coco_add_pars(prob_p,'xy',[1,2,3],{'x','y','u3'});
prob_pxy=coco_add_event(prob_pxy,'UZ','x',-2:0.1:2);
prob_pxy=coco_set(prob_pxy,'cont','NPR',inf,'ItMX',40);
bd_w_ev=coco(prob_pxy,'pxy_w_event',[],1,{'x','y','u3','p'}); % swap x and y position to see effect on 'FP' outputs
%% Process bifurcation diagram with coco_bd_ functions
% (there are many more coco_bd_ functions)
bda=coco_bd_table('pxy_w_event');   % read from data/pxy_w_event/bd.mat
xy=bda{:,{'x','y'}};     % get columns of cell array headed by x and y
getidx=@(bd,type)find(strcmp(bd.TYPE,type)); % extract rows with special point types
uzidx=getidx(bda,'UZ'); % extract UZ rows
fpidx=getidx(bda,'FP'); % extract FP rows
figure(1);clf;
lw={'Linewidth',2};
plot(xy(:,1),xy(:,2),'.-',...
    xy(uzidx,1),xy(uzidx,2),'o',...
    xy(fpidx,1),xy(fpidx,2),'ks',lw{:});grid on
legend({'branch','UZ','FP'});
set(gca,'FontSize',18,'DataAspectRatio',[1,1,1]);
xlabel('x');ylabel('y')
%%
%%
%% functions for right-hand side for circle
function [data,y]=circle(prob,data,u) %#ok<INUSD>
y=u(1)^2+(u(2)-1)^2-1;
end
%%  functions for right-hand side for circle with Jacobian provided
function [data,y,J]=circle_J(~,data,u)
y=u(1)^2+(u(2)-1)^2-1;
J=2*[u(1),u(2)-1];
end
