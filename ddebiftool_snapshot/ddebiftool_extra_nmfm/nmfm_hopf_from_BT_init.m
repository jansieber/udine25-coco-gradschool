function [hopf_br,augmented] = nmfm_hopf_from_BT_init(funcs,bt,eps,free_pars,varargin)
%% Initialize branch for continuing the homoclinic curve
% emanating from the generic/transcritical Bogdanov-Takens point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: nmfm_hopf_from_BT_init 314 2019-01-24 14:28:23Z mmbosschaert $

%% Process options
default = {'debug',1,'perturbationparameter',0.01,...
    'generic_unfolding',true,'step',1.1,'correc',false,'dir',0};
[options,pass_on] = dde_set_options(default,varargin,'pass_on');

%% critical normal form coefficients
a  = bt.nmfm.a;
b  = bt.nmfm.b;

%% centermanifold transformation
H = bt.nmfm.H;
K = bt.nmfm.K;

hopf_br1 = df_brnch(funcs,free_pars,'hopf');
hopf_br1 = replace_branch_pars(hopf_br1,free_pars,pass_on);
if ~options.generic_unfolding
    hopf_br2 = hopf_br1;
end
hopf = repmat(p_tohopf(funcs,bt),1,length(eps));
mhopf=df_mthod(funcs,'hopf');

for i=1:size(eps,2)
    if options.generic_unfolding
        hopf(i).parameter(free_pars) = hopf(i).parameter(free_pars) + K(-eps(1,i)^4/(4*a),b/(2*a)*eps(1,i)^2)';
        hopf(i).x = hopf(i).x + H(-eps(1,i)^2/(2*a),0,-eps(1,i)^4/(4*a),b/(2*a)*eps(1,i)^2);
    else
        hopf(i).parameter(free_pars) = hopf(i).parameter(free_pars) + K(-eps(1,i)^2,0)';
    end
    hopf(i).omega = eps(1,i);
    [p,q,~]=nmfm_nullvector(funcs,hopf(i),eps(1,i));
    hopf(i).v = q;
end
hopf_br1.point=hopf;
if options.generic_unfolding
    hopf_br = hopf_br1;
    augmented=true;
    return
end

% initialize Hopf branch to the non-trivial equilibrium for the transcritical Bogdanov-Takens bifurcation
hopf = repmat(p_tohopf(funcs,bt),1,length(eps));
if size(eps,1) == 1
    j=1;
else
    j=2;
end
for i=1:size(eps,2)
    hopf(i).parameter(free_pars) = hopf(i).parameter(free_pars) + K(eps(j,i).^2,b/a*eps(j,i).^2)';
    alpha =  hopf(i).parameter(free_pars);
    hopf(i).x = hopf(i).x + H(-eps(j,i)^2/a,0,eps(j,i)^2,b/a*eps(j,i)^2);
    hopf(i).omega = eps(j,i);
    [p,q,~]=nmfm_nullvector(funcs,hopf(i),eps(j,i));
    hopf(i).v = q;
end
hopf_br2.point=hopf;

hopf_br = [hopf_br1,hopf_br2];
augmented=true;
end
