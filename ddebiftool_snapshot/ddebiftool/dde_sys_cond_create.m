function [cond,name]=dde_sys_cond_create(varargin)
%% create sys_cond structure
% named inputs: 'name' used for debugging, 'fun' function handle to
% [r,J]=sys_cond(p[,pref]) or arbitrary function
% y=f(x1,x2,..,xref1,xref2,...), where format of x1,x2,... is controlled by
% 'args'.
% 
% 'args': cell array, sequence of field names in ind structure as obtained
% from ind_from_point(point,1:length(point.parameter)) and indices. For
% example, if point has kind 'hopf', args={'parameter',3:4,'v','re',1}
% implies that f has the form f(p,v) where p is point.parameter(3:4) and v
% is real(point.v(1)). If point has kind 'psol',
% args={'profile',{1:3,1},'parameter',1:3,'period',1}, then f has form
% f(x,p,T), where x=profile(1:3,1), p=parameter(1:3) and T=period.
%
% 'deriv': number for numerical differentiation (stepsize) or function for
% returning directional derivative of f.
%
% 'isvec' flag indicating if f is vectorized (only relevant if deriv is
% numeric for finite differences)
%% call with name-value pairs for construction of single sys_cond
if isempty(varargin)||ischar(varargin{1})
    [cond,name]=loc_sys_cond_create(varargin{:});
    return
end
%% input is already array of sys_cond's
if isstruct(varargin{1})
    cond=varargin{1};
    name={cond.name};
    if isscalar(cond)
        name=name{1};
    end
    return
end
%% input is cell array or tuple of functions that need conversion to sys_cond
if isscalar(varargin) && iscell(varargin{1})
    args=varargin{1};
else
    args=varargin;
end
cond=repmat(dde_sys_cond_create(),length(args),1);
name=cell(length(args),1);
for i=1:length(args)
    name{i}=sprintf('sys_cond_%d',i);
    cond(i)=loc_sys_cond_create('name',name{i},'fun',args{i},'reference',nargin(args{i})==2);
end
end
function [cond,name]=loc_sys_cond_create(varargin)
default={'name','dde_dummy_cond','fun',@dde_dummy_cond,'reference',false,...
    'args',{},'deriv',1e-6,'isvec',true};
[options,dum,used]=dde_set_options(default,varargin,'pass_on'); %#ok<ASGLU>
cond=struct('name',options.name,'reference',options.reference,...
    'fun',options.fun);
name=cond.name;
if ~isempty(options.args)
    if ~used.name && isa(options.fun,'function_handle')
        cond.name=func2str(options.fun);
    end
    cond.fun=@(p,pref)gen_cond(options.fun,p,pref,options.reference,...
        options.args,options.deriv,options.isvec);
    cond.reference=true;
end
end
%% dummy sys_cond condition
function [resi,condi]=dde_dummy_cond(point,pref) %#ok<INUSD>
resi=zeros(0,1);
condi=repmat(point,0,1);
end
%% general simplified interface for sys_cond
function [r,J]=gen_cond(fun,p,pref,hasref,args,deriv,isvec)
freepar=1:length(p.parameter);
[ip,xlen]=dde_ind_from_point(p,freepar);
x=dde_x_from_point(p,freepar);
if hasref
    xref=dde_x_from_point(pref,freepar);
end
indxfield=ip;
nargs=sum(cellfun(@(a)iscell(a)|isnumeric(a),args));
[xarg,iparg,ipall]=deal(cell(1,nargs));
fmt=ones(1,nargs);
if hasref
    xargref=xarg;
else
    xargref=cell(1,0);
end
ia=0;
for i=1:length(args)
    if ischar(args{i})
        indxfield=indxfield.(args{i});
    elseif iscell(args{i})||isnumeric(args{i})
        ia=ia+1;
        if isempty(args{i})
            sz=size(indxfield);
            args{i}=arrayfun(@(a){1:a},sz);
        elseif isnumeric(args{i})
            args{i}=args(i);
        end
        iparg{ia}=args{i};
        fmt(ia)=numel(iparg{ia});
        ipall{ia}=indxfield(iparg{ia}{1});
        [xarg{ia},ipall{ia}]=extract_ind(x,indxfield,args{i});
        if hasref
            xargref{ia}=extract_ind(xref,indxfield,args{i});
        end
        indxfield=ip;
    end
end
r=fun(xarg{:},xargref{:});
r=r(:);
ncond=numel(r);
%% if derivative not provided, use finite differences
if isnumeric(deriv)
    fun1=@(varargin)fun(varargin{:},xargref{:});
    vec_fun = @(varargin)fun_vec(fun1,fmt,isvec,ncond,varargin{:});
    df=@(varargin)num_dirderiv(vec_fun,fmt,varargin{:},'hjac',deriv);
else
    df=@(varargin)deriv(varargin{1:end/2},xargref{:},varargin{end/2+1:end});
end
%% find all partial derivatives
dev0=num2cell(zeros(1,ia));
Jnum=cell(1,ia);
for i=1:ia
    dev=dev0;
    dev{i}={1,'I'};
    Jnum{i}=mult_deriv(df,fmt,xarg{:},dev{:});
end
%% assign partial derivatives to point structures
Jx=zeros(ncond,xlen);%repmat(p_axpy(0,p,[]),ncond,1);
for i=1:ia
    Jx(:,ipall{i})=Jx(:,ipall{i})+reshape(Jnum{i},[],numel(ipall{i}));
end
J=dde_point_from_x(Jx',p_axpy(0,p,[]),freepar);
end
function [xarg,indarg]=extract_ind(x,indxfield,inds)
xargfield=reshape(x(indxfield(:)),size(indxfield));
if ~isempty(inds)
    xarg=subsref(xargfield,struct('type','()','subs',{inds}));
    indarg=subsref(indxfield,struct('type','()','subs',{inds}));
else
    xarg=xargfield;
end
end