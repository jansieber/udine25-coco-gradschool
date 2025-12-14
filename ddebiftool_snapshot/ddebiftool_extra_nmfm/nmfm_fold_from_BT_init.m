function [fold_br,augmented] = nmfm_fold_from_BT_init(funcs,bt,eps,free_pars,varargin)
%% Initialize branch for continuing the fold curve
% emanating from the generic/transcritical Bogdanov-Takens point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_BT_H.m 314 2019-01-24 14:28:23Z mmbosschaert $

%% Process options
default = {'debug',1,'perturbationparameter',0.01,...
    'generic_unfolding',true,'step',1.1,'correc',true,'dir',1};
[options,pass_on] = dde_set_options(default,varargin,'pass_on');

%% critical normal form coefficients
a  = bt.nmfm.a;
b  = bt.nmfm.b;

%% centermanifold transformation
H = bt.nmfm.H;
K = bt.nmfm.K;

fold_br = df_brnch(funcs,free_pars,'fold');
fold = repmat(p_tofold(funcs,bt),1,length(eps));

for i=1:length(eps)
    if options.dir == 1
        fold(i).parameter(free_pars) = fold(i).parameter(free_pars) + K(0,eps(i))';
        [p,q,~]=nmfm_nullvector(funcs,fold(i),eps(i));
        fold(i).v = q;
    else
        % for some reason 1e-08 is needed for some system to converge
        fold(i).parameter(free_pars) = fold(i).parameter(free_pars) + K(1e-08,eps(i))';
    end
end
fold_br.point=fold;
augmented = true;
end
