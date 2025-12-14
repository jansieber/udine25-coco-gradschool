function [J,res,tT,extmesh]=dde_coll_jac(funcs,pt,free_par,varargin)
%% return Jacobian and residual of collocation problem
% function [J,res,tT,extmesh]=coll_dde_jac(funcs,pt,free_par,...)
%% INPUT:
% *  funcs problem functions
% * pt point structure (psol or hcli)
% * free_par indices of free parameters numbers in N^np
% 
%% OUTPUT:
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+np)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% modified from psol_jac to permit arbitrary nesting of state-dependent
% delays, vectorisation and re-use for computation of Floquet multipliers,
% extra conditions for state-dependent delay (p_tau, etc), connecting
% orbits
%
%% Optional inputs:
%
% * c collocation parameters in [0,1]^m
% * If 'wrapJ' (default true) the returned jacobian is augmented with
% derivative wrt period and free_par and wrapped around periodically inside
% [0,1] if ~wrapJ the returned jacobian puts its entries into the interval
%  [-min(delay,0),max(1-[0,delays])], no augmentation is done. This
%  Jacobian can be used for computation of Floquet multipliers and modes
%
% If 'c_is_tvals' (default false) then the values in c (must be not empty) are
% taken as those values of time in [0,1] of the entire period, where
% residual and Jacobian are calculated.
%
% The argument 'Dtmat' (default=Id) is multiplied with the time
% derivative. Dtmat==zeros() permits algebraic equations.
%
% The argument
%%
%% optional (args are also passed on to coll_dde_map)
default={'Dtmat',funcs.lhs_matrixfun(size(pt.profile,1)),'output','J','matrix','sparse',...
    'wrapJ',strcmp(pt.kind,'psol'),'max_xderiv',1,'maxsize',2^23};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% are delayed arguments derivatives?
funcs.lhs_matrixfun=@(nx)options.Dtmat;
[ntaum1,max_rhs_xderiv]=dde_num_delays(funcs); %#ok<ASGLU>
max_xderiv=max_rhs_xderiv+max(1,options.max_xderiv);
%% get values at collocation points
[tc,submesh_limit]=coll_map_times(funcs,pt,pass_on{:});
[tau,xxd]=dde_coll_delays(funcs,pt,'t',tc,...
    'submesh_limit',submesh_limit,pass_on{:},'max_xderiv',max_xderiv);
fmt=[3,2,1];
baseargs={xxd,pt.parameter,pt.period};
arg_exp=@(varargin)arg_array_expand(fmt,baseargs{:},varargin{:});
%% residual
res=dde_coll_rhs(funcs,arg_exp());
dres=dde_coll_bd_diff(pt,pass_on{:});
nforce_smooth=numel(dres);
res=[res(:);dres(:)];
J.nforce_smooth=nforce_smooth;
J.tc=tc;
if strcmp(options.output,'res')
    [J,res]=deal(res,J);
    return
end
[nx,nderivp1]=size(xxd,[1,3]);% dim x and number of derivatives provided by coll_delays
tT=tau/pt.period;
[W,ptext]=coll_map(pt,tc,tT,nderivp1,options.wrapJ,submesh_limit,pass_on{:});
extmesh=ptext.mesh;
%% find derivatives wrt to profile points
xpattern=add_Mxp_pattern(funcs,nx);
[argstruc,fillfun]=dde_xjac_expand(xpattern,fmt,true,options.maxsize,baseargs{:});
for i=length(argstruc.chunks):-1:1
    argx=argstruc.fun(argstruc.chunks{i});
    [argx,dimx]=arg_flatten(fmt,argx);
    Jxc{i}=dde_coll_rhs(funcs,argx{:});
    Jxc{i}=reshape(Jxc{i},[size(Jxc{i},1),dimx]);
end
Jxall=fillfun(cat(2,Jxc{:}));
dtJprofile=dde_coll_bd_diff(ptext,'output','matrix',pass_on{:});
J.profile=cat(1,sparse_blkdiag(reshape(Jxall,size(Jxall,1),[],size(Jxall,5)))*W,...
    dtJprofile);
assert(size(J.profile,1)==length(res));
if ~strcmp(options.matrix,'sparse')
    J.profile=full(J.profile);
end
if strcmp(options.output,'profile')
    J=J.profile;
    return
end
%% derivative wrt parameter p and period T
[argp,dimp,nvecp]=arg_flatten(fmt,arg_exp(0,{1,{free_par}},0));
[argT,dimT,nvecT]=arg_flatten(fmt,arg_exp(0,0,1));
all_c=reshape(permute(cat(3,argp,argT),[3,2,1]),2,[]);
all_args=arrayfun(@(i){cat(ndims(all_c{1,i}),all_c{:,i})},1:size(all_c,2));
all_args=permute(reshape(all_args,2,[]),[2,1]);
JpT=dde_coll_rhs(funcs,all_args{:});
JpTc=cellfun(@(a,f){reshape(a,[size(a,1),f])},...
    mat2cell(JpT,size(JpT,1),[nvecp,nvecT]),{dimp,dimT});
%% derivative wrt parameter
J.parameter=reshape(permute(JpTc{1},[1,3,2]),[],length(free_par));
J.parameter=cat(1,J.parameter,zeros(nforce_smooth,length(free_par)));
%% derivative wrt period
J.period=reshape(permute(JpTc{2},[1,3,2]),[],1);
J.period=cat(1,J.period,zeros(nforce_smooth,1));
assert(...
    (isempty(J.period)   ||   size(J.period,1)==length(res)) && ...
    (isempty(J.parameter)||size(J.parameter,1)==length(res)))
end
%%
function [W,ptext]=coll_map(pt,tc,tT,nderivp1,wrapJ,submesh_limit,varargin)
[d,neqs]=size(tT);
nx=size(pt.profile,1);
t_eval=tc(ones(d,1),:)-tT;
Wc=cell(nderivp1,1);
for k=0:nderivp1-1
    [Wc{k+1},dum,ptext]=dde_coll_wrap_eva(pt,t_eval(:).','kron',true,varargin{:},...
        'diff',k,'wrapJ',wrapJ,'submesh_limit',submesh_limit,'output','matrix'); %#ok<ASGLU>
    Wc{k+1}=reshape(Wc{k+1},nx*d,[]);
end
W=reshape(cat(1,Wc{:}),nx*d*nderivp1*neqs,[]);
end
%%
function [tc,submesh_limit] = coll_map_times(funcs,pt,varargin)
%% Determine where residuals are computed: mesh tc 
% only compute res & jac at specified points?
default={'c_is_tvals',false,{'c','collocation_parameters'},[],'submesh_limit',0,'include0',false,...
    'mesh_range',1:pt.degree:length(pt.mesh)};
options=dde_set_options(default,varargin,'pass_on');
submesh_limit=options.submesh_limit;
c=options.c;
if ~options.c_is_tvals
    if isempty(c) || ischar(c)
        c=dde_coll_set_grid('collocation',pt.degree,'type',options.c,...
            'lhs_num',funcs.lhs_matrixfun(size(pt.profile,1)));
    end
    if 1==max(c)
        submesh_limit=1;
    end
    tc=dde_coll_meshfill(pt.mesh(options.mesh_range),pt.degree,...
        'purpose','collocation','grid',c);
else
    tc=c(:)';
end
if options.include0 && strcmp(options.c,'force_smooth')
    tc=sort([0,tc]);
end
end
%% add Mx' to funcs.xpattern
function xpat=add_Mxp_pattern(funcs,nx)
orig_xpattern=dde_funcs_xpattern(funcs,nx);
M=funcs.lhs_matrixfun(nx);
ixd=find(any(M~=0,1));
Mxp_pattern=[ixd;[1;2]*ones(1,size(ixd,2))];
xpat=dde_join_xpattern(Mxp_pattern,orig_xpattern);
end