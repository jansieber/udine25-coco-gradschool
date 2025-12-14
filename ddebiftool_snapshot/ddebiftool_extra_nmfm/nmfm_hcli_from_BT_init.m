function [hcli_br,augmented] = nmfm_hcli_from_BT_init(funcs,bt,eps,free_pars,varargin)
%% Initialize branch for continuing the homoclinic curve(s)
% emanating from the generic/transcritical Bogdanov-Takens point.
%
% @author Maikel Bosschaert, maikel.bosschaert -at- uhasselt.be
% @Id $Id: init_BT_H.m 314 2019-01-24 14:28:23Z mmbosschaert $

%% Process options
default = {'debug',1,'generic_unfolding',true, 'correc',false,'order',3};
[options,pass_on] = dde_set_options(default,varargin,'pass_on');

if length(eps)==1
    eps(2)=eps(1)*1.1;
end
hcli_br1 = df_brnch(funcs,free_pars,'hcli');
hcli_br1 = replace_branch_pars(hcli_br1,free_pars,pass_on);
if options.generic_unfolding
    number_of_branches = 1;
else
    number_of_branches = 2;
end
hcli_br = repmat(hcli_br1, 1, number_of_branches);
mhcli=df_mthod(funcs,'hcli');

options.free_pars = free_pars;
for i=1:size(eps,2)
    options.perturbationparameter = eps(:,i);
    hcli = init_BT_Hom_orbital(funcs,bt,options);
    for j=1:number_of_branches
        if options.generic_unfolding
            hcli_pred = hcli;
            hcli_pred = hcli_pred{1};
        else
            hcli_pred = hcli(j);
            hcli_pred = hcli_pred{1};
        end
        if options.correc
            secant = p_tangent(funcs,mhcli.point,hcli_pred,free_pars);
            [hcli_correc,s] = p_correc(funcs,hcli_pred,free_pars,secant,mhcli.point);
            if ~s
                return
            end
            hcli_br(j).point = [hcli_br(j).point; hcli_correc];
        else
            hcli_br(j).point = [hcli_br(j).point; hcli_pred];
        end
    end
end

augmented = true;
end
