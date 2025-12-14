%% determine tangent to rhs in a point
%
% inputs:
%
% * |funcs|: right-hand side structure
% * |mth|: point method structure of branch (for Jacobian)
% * |pt|: point in which tangent is taken
% * |free_par|: indices of free parameters
%
% output:
% |pnull| array of point structures of same type as pt spanning the
% nullspace of the Jacobian of the zero problem defined by funcs, pt and
% mth.
%
% This routine assumes that the |nr x nc| Jacobian J has full rank such
% that the array |pnull| has length |nc-nr|.
%
% $Id: p_tangent.m 369 2019-08-27 00:07:02Z jansieber $
%%
function pnull=p_tangent(funcs,mth,pt,free_par,varargin)
default={'tangent_dim',1,'border',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if strcmp(mth.preprocess,'dde_jac2square_preprocess')
    mth.preprocess='';
    mth.postprocess='';
end
border=options.border;
if isstruct(border)
    border=dde_x_from_point(border,free_par);
    for i=1:size(border,2)
        border(:,i)=border(:,i)/norm(border(:,i),'inf');
    end
end
[fJ,x0,data]=p_correc_setup(funcs,pt,free_par,mth,'previous',pt,'output','J',...
    'row_minus_column',options.tangent_dim,pass_on{:});
J=fJ(x0);
nulldim=size(J,2)-size(J,1); % assume full rank
if all(size(border)==[size(J,2),nulldim])
    rhs=cat(1,zeros(size(J,1),nulldim),eye(nulldim));
    jnull=cat(1,J,border.')\rhs;
    [jnull,dum]=qr(full(jnull),0); %#ok<ASGLU>
else
    jnull=dde_nullspaces_lr(J,'nulldim',nulldim);
end
for i=1:size(jnull,2)
    jnull(:,i)=jnull(:,i)/norm(jnull(:,i),'inf');
end
pnull=dde_point_from_x(jnull,p_axpy(0,data.point,[]),data.free_par);
pnull=dde_postprocess(pnull,data);
end