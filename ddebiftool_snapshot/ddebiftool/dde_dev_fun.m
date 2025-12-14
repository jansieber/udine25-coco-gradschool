function out=dde_dev_fun(cmd,varargin)
%% History function of form sum_k vk theta^(ak) exp(lambda_k theta)
%
% x=sum_2 fun.v(:,i,:) * theta.^fun.t(i,:) * exp(fun.lambda(i,:)*theta)
switch cmd
    case {'fun','define','set'}
        out=loc_dev_fun(varargin{:});
    case 'ax'
        out=loc_dev_ax(varargin{:});
    case {'matmul','*'}
        out=loc_dev_matmul(varargin{:});
    case {'call','','eval'}
        out=loc_dev_call(varargin{:});
    case {'conj'}
        out=loc_dev_conj(varargin{:});
    case {'dtheta','d_theta','dt'}
        out=loc_dev_dtheta(varargin{:});
    case {'deriv1', 'd1dt'}
        out=loc_dev_deriv1(varargin{:});
    case 'dvlam'
        out=loc_dev_dvlam(varargin{:});
    case {'==','equal','is_equal'}
        out=loc_dev_equal(varargin{:});
    case {'is0','==0'}
        out=loc_dev_equal(varargin{1},setfield(varargin{1},'v',varargin{1}.v*0)); %#ok<SFLD>
    case 'group'
        out=loc_dev_group(varargin{:});
    otherwise
        error('dde_dev_fun:cmd','dde_dev_fun: unknown command');
end
end
%%
function fun=loc_dev_fun(in1,varargin)
%% Create history function of form sum vk theta^(ak) exp(lambda_k theta)
%
% x=sum_2 fun.v(:,i,:) * theta.^fun.t(1,i,:) * exp(fun.lambda(1,i,:)*theta)
%
if isstruct(in1)
    [v,t,lambda]=deal(in1.v,in1.t,in1.lambda);
else
    v=in1;
    %% Optional arguments
    default={'lambda',0,'t',0};
    options=dde_set_options(default,varargin,'pass_on');
    lambda=rs_rep(v, options.lambda);
    t=rs_rep(v,options.t);
end
%% combine duplicates (if lambda_k and t_k are exactly equal)
[nv,nvec]=size(v,2:3);
tlam=reshape(permute(cat(1,t,lambda),[2,1,3]),nv,[]);
[tl,dum,ic]=unique(tlam,'rows'); %#ok<ASGLU>
tl=reshape(tl,[],2,nvec);
ic=reshape(ic,nv,2,nvec);
tnew=reshape(tl(:,1,:),1,[],nvec);
lambdanew=reshape(tl(:,2,:),1,[],nvec);
vnew=zeros(size(v,1),size(lambdanew,2));
for i=1:nv
    vnew(:,ic(i),:)=vnew(:,ic(i),:)+v(:,i,:);
end
vnz=any(any(vnew~=0,1),3);
if any(vnz) %|| isempty(vnz)
    fun=struct('v',vnew(:,vnz,:),'lambda',lambdanew(1,vnz,:),'t',tnew(1,vnz,:));
else
    fun=struct('v',vnew(:,1,:),'lambda',zeros(1,1,nvec),'t',zeros(1,1,nvec));
end
end
%%
function out=loc_dev_matmul(mat,funs)
for i=size(mat,1):-1:1
    out(i)=loc_dev_ax(mat(i,:),funs);
end
end
%%    
function fun=loc_dev_ax(avec,funs)
%% linear combination of history functions
%
%%
nf=length(avec);
for i=1:nf
    funs(i).v=avec(i)*funs(i).v;
end
vvec=cat(2,funs.v);
lvec=cat(2,funs.lambda);
tvec=cat(2,funs.t);
fun=loc_dev_fun(struct('v',vvec,'lambda',lvec,'t',tvec));
end
%%
function x=loc_dev_call(fun,theta,deriv,rg)
%% Evaluate history function of form sum vk t^(ak) exp(lambda_k theta)
%
% x=sum fun.v(:,i) * theta.^fun.t(i) * exp(fun.lambda(i)*theta)
%
% if present, 3rd arg deriv (0 default): differentiate so many times before
% evaluating
%
%%
if nargin>=3
    for i=1:deriv
        fun=loc_dev_deriv1(fun);
    end
end
[n,nv,nvec]=size(fun.v);
if nargin<4
    rg=1:n;
end
n=length(rg);
nt=length(theta);
ont=ones(nt,1);
on=ones(n,1);
theta=reshape(theta,1,1,1,nt);
theta=theta(1,ones(nv,1),ones(nvec,1),:);
v=fun.v(rg,:,:,ont);
lam=fun.lambda(1,:,:,ont);
explamtheta=ones(1,nv,nvec,nt);
explamtheta(lam~=0)=exp(lam(lam~=0).*theta(lam~=0));
explamtheta=explamtheta(on,:,:,:);
theta=theta(on,:,:,:);
tp=fun.t(on,:,:,ont);
tp_gt_0=tp>0;
ttp=ones(n,nv,nvec,nt);
ttp(tp_gt_0)=theta(tp_gt_0).^tp(tp_gt_0);
s=v.*explamtheta.*ttp;
x=permute(reshape(sum(s,2),[n,nvec,nt]),[1,3,2]);
end
%%
function fout=loc_dev_conj(fin)
%% create complex conjugate of history function
%
%% 
fout=loc_dev_fun(...
    struct('v',conj(fin.v),'lambda',conj(fin.lambda),'t',fin.t));
end
%%
function dfun=loc_dev_dtheta(fun,order)
dfun=fun;
if nargin==1
    order=1;
end
for i=1:order
    dfun=loc_dev_deriv1(dfun);
end
end
%%
function dfun=loc_dev_deriv1(fun)
%% Differentiate history function (once) wrt theta
%
%%
n=size(fun.v,1);
lambda=fun.lambda;
v=fun.v.*lambda(ones(n,1),:,:);
t=fun.t;
t_gt_0=t>0;
if ~any(t_gt_0)
    dfun=loc_dev_fun(struct('v',v,'lambda',lambda,'t',t));
    return
end
tex=zeros(size(t));
tex(1,t_gt_0)=t(1,t_gt_0);
vex=zeros(size(fun.v));
vex(:,t_gt_0)=fun.v(:,t_gt_0).*tex(ones(n,1),t_gt_0);
lex=zeros(size(lambda));
lex(1,t_gt_0)=lambda(1,t_gt_0);
tex=max(tex-1,0);
v=cat(2,v,vex);
lambda=cat(2,lambda,lex);
t=cat(2,t,tex);
dfun=loc_dev_fun(struct('v',v,'lambda',lambda,'t',t));
end
%%
function dfun=loc_dev_dvlam(fun,ord,dv,dlam)
%% Directional deriv of history function wrt v and lam of order ord
% in direction dv (size(v)), dlam (size(lambda))
%%
if ord==0
    dfun=fun;
    return
end
[n,nvec]=size(fun.v,[1,3]);
dv=repmat(dv,1,nv/size(dv,2),nvec/size(dv,3));
dlam=rs_rep(fun.v,dlam);
dlamkm1=dlam.^(ord-1);
dlamk=dlam.^ord;
v1=ord*dv.*dlamkm1(ones(n,1),:,:);
t1=fun.t+ord-1;
t2=fun.t+ord;
v2=fun.v.*dlamk(ones(n,1),:,:);
dfun=loc_dev_fun(struct('v',cat(2,v1,v2),...
    'lambda',cat(2,fun.lambda,fun.lambda),'t',cat(2,t1,t2)));
end
%%
function iseq=loc_dev_equal(f1,f2)
%% Check if two history functions are (exactly) equal
%
%%
iseq=false;
f1=loc_dev_fun(f1);% compress if duplicates
f2=loc_dev_fun(f2);% compress if duplicates
[n1,nv1]=size(f1.v);
[n2,nv2]=size(f2.v);
if n2~=n1 || nv1~=nv2
    return
end
if any(f1.v(:)~=f2.v(:)) || any(f1.lambda~=f2.lambda) || any(f1.t~=f2.t)
    return
end
iseq=true;
end
%%
function groups=loc_dev_group(funs)
%% Group history functions into groups of equals
% output cell array of index sets where each cell contains indices with
% equal funs.
%%
nf=length(funs);
remaining=1:nf;
groups=cell(1,nf);
count=0;
%% for test purposes (such that funs can be complex numbers or dev functions)
if isnumeric(funs)
    iseq=@(x,y)all(x(:)==y(:));
else
    iseq=@nmfm_dev_equal;
end
%% collect into groups
while ~isempty(remaining)
    count=count+1;
    groups{count}=remaining(1);
    remaining(1)=[];
    for i=length(remaining):-1:1
        if iseq(funs(groups{count}(1)),funs(remaining(i)))
            groups{count}(end+1)=remaining(i);
            remaining(i)=[];
        end
    end
end
groups=groups(1:count);
end

function out = rs_rep(v, arg)
%% reshape and repmat argument to shape 1 x size(v,2) x size(v,3)
out=arg;
[nv,nvec]=size(v,2:3);
if isempty(out)
    out=0;
end
if size(out,1)>1
    out=reshape(out,1,size(out,1),[]);
end
out=repmat(out,1,nv/size(out,2),nvec/size(out,3));
end