function relativeErrors = convergence_plot(funcs, bt, orders, amplitudes, varargin)

    default = {'dir',0,'free_pars',[],'orders',3,'TTolerance',1e-5,'ntst',40,'generic_unfolding',1,'debug',false};
    [options,pass_on] = dde_set_options(default,varargin,'pass_on');

    mhcli=df_mthod(funcs,'hcli');
    relativeErrors = {};
    for i=orders
        relativeErrors{i} = zeros(size(amplitudes));
        options.order = i;
        for j=1:length(amplitudes)
            if options.generic_unfolding
                options.amplitude = amplitudes(j);
            else
                options.amplitude = amplitudes(j)*[1,1];
                options.TTolerance = options.TTolerance.*[1,1];
            end
            hcli = init_BT_Hom_orbital(funcs,bt,options);
            if length(hcli) == 1
                hcli_pred = hcli;
            else
                hcli_pred = hcli(2);
            end
            hcli_pred = hcli_pred{1};
            if options.dir > 0
                [hcli_corrected,s]=p_correc(funcs,hcli_pred,options.free_pars(options.dir),[],mhcli.point);
            else
                secant=p_tangent(funcs,mhcli.point,hcli_pred,options.free_pars);
                [hcli_corrected,s]=p_correc(funcs,hcli_pred,options.free_pars,secant,mhcli.point);
            end
            if s
                relativeErrors{i}(j) = norm(hcli_corrected.profile-hcli_pred.profile,Inf)/norm(hcli_corrected.profile,Inf);
            end
      end
    end

end
