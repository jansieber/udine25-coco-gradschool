function varargout=dde_psol_lincond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_psol_lincond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_psol_lincond(p,varargin{2:end}));
end
end
function [r,J]=loc_psol_lincond(point,varargin)
if isscalar(varargin)
    S=varargin{1};
elseif isnumeric(varargin{1})
    S=dde_lincond_struct(varargin{:});
else
    S=dde_lincond_struct(size(point.profile,1),varargin{:});    
end
if ~strcmp(S.fieldname,'profile')
    S=dde_psol_lincond_correct(point,S);
end
%% linear condition for kind psol
% impose linear condition on solution
ns=size(S.condprojmat,1);
fn=S.fieldname;
%% check if condition is empty
if ~isfield(point,fn) || ns==0
    r=zeros(0,1);
    J=repmat(point,0,1);
    return
end
%% extract solution profile & determine matrix dimensions
x=point.(fn);
[xdim,nt]=size(x);
[np,nxs]=size(S.stateproj);
%% append zero matrix if stateproj does not fit profile dimension
if nxs<xdim
    S.stateproj=[S.stateproj,zeros(size(S.stateproj,1),xdim-nxs)];
end
%% create rotation matrix (inserted after trafo)
op=ones(np,1);
frac=S.rotation(1)/S.rotation(2);
rot=exp(2*pi*1i*frac);
if S.rotation(1)==S.rotation(2)||S.rotation(1)==0
    matrot=eye(np);
else
    op=op(1:end/2);
    rrot=diag(real(rot(op,1)));
    irot=diag(imag(rot(op,1)));
    matrot=[rrot,-irot; irot,rrot];
end
%% check which conditions are at discrete times and which are integrals
ms=S.condprojint;
nc=size(ms,1);
isint=find(ms(:,1)<ms(:,2));
nint=length(isint);
ispt=find(ms(:,1)==ms(:,2));
npt=length(ispt);
tshift=S.shift(1)/S.shift(2);
Jprof=NaN(ns,nc,xdim*nt);
%% discrete part of conditions
tpt=[ms(ispt,1);ms(ispt,1)+tshift];
Jpt0=dde_coll_eva(point,tpt.','kron',true,'output','matrix');
Jprof(:,ispt,:)=mshape(Jpt0,'rs',{xdim,npt,2,xdim*nt},'p',[1,3,2,4],'rs',{xdim*2,npt*xdim*nt},...
   'mat', 'apply',@(M)S.condprojmat*[S.trafo*matrot*S.stateproj,-S.stateproj]*M,...
    'struct','rs',{ns,npt,xdim*nt},'state','full','mat');
%% integral part of conditions
for i=1:nint
    J0=dde_coll_int(point,ms(isint(i,:)),       'output','matrix','kron',true);
    J1=dde_coll_int(point,ms(isint(i,:))+tshift,'output','matrix','kron',true);
    Jloc=S.condprojmat*(S.trafo*matrot*S.stateproj*J0-S.stateproj*J1);
    Jprof(:,isint(i),:)=reshape(full(Jloc),ns,1,xdim*nt);
end
%% residual
r=reshape(Jprof,ns*nc,xdim*nt)*x(:);
if ~isempty(S.res_parameters)
    par_sel=S.res_parameters~=0;
    par=zeros(size(r));
    par(par_sel)=point.parameter(S.res_parameters(par_sel));
    r=r-par;
end
%% Jacobian
Jprof=permute(reshape(Jprof,ns,nc,xdim,nt),[3,4,1,2]);
J=p_axpy(0,point,[]);
J=repmat(J,ns,nc);
for k=1:nc
    for i=1:ns
        J(i,k).(fn)=Jprof(:,:,i,k);
        if ~isempty(S.res_parameters)&& par_sel(i,k)
            J(i,k).parameter(S.res_parameters(i,k))=-1;
        end
    end
end
J=J(:);
end
