function varargout=dde_stst_lincond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_stst_lincond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_stst_lincond(p,varargin{2:end}));
end
end
%%
function [r,J]=loc_stst_lincond(point,varargin)
if length(varargin)==1
    S=varargin{1};
elseif isnumeric(varargin{1})
    S=dde_lincond_struct(varargin{:});
else
    S=dde_lincond_struct(size(point.x,1),varargin{:});    
end
nc=size(S.condprojmat,1);
fn=S.fieldname;
if ~isfield(point,fn) || nc==0
    r=zeros(0,1);
    J=repmat(point,0,1);
    return
end
x=point.(fn);
[np,nxs]=size(S.stateproj);
if nxs<size(x,1)
    S.stateproj=[S.stateproj,zeros(size(S.stateproj,1),size(x,1)-nxs)];
end
frac=S.rotation(1)/S.rotation(2);
rot=exp(2*pi*1i*frac);
if 2*frac==round(2*frac)
    rot=sign(real(rot));
end
matrot=diag(rot(ones(np,1),1));
matsym=S.condprojmat*(S.trafo-matrot)*S.stateproj;
r=matsym*x;
if ~isempty(S.res_parameters)
    r=r-reshape(point.parameter(S.res_parameters),[],1);
end
J0=p_axpy(0,point,[]);
J=repmat(J0,nc,1);
for i=1:nc
    J(i).(fn)=reshape(matsym(i,:)',size(x,1),1);
    if ~isempty(S.res_parameters)
        J(i).parameter(S.res_parameters(i))=-1;
    end
end
ind=dde_ind_from_point(point,1:length(point.parameter));
ix=ind.(fn);
if isstruct(ix)
    r=[real(r);imag(r)];
    J(nc+(1:nc))=J0;
    for i=nc:-1:1
        J(i+nc).(fn)=reshape(1i*J(i).(fn),size(x,1),1);
    end
end
end
