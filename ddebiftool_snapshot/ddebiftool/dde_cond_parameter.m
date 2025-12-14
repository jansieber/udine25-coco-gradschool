function varargout=dde_cond_parameter(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_cond_parameter(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_cond_parameter(p,varargin{2:end}));
end
end
function [r,J]=loc_cond_parameter(p,minbounds,maxbounds)
rmin=minbounds(:,2)-reshape(p.parameter(minbounds(:,1)),[],1);
if nargin>2
    rmax=reshape(p.parameter(maxbounds(:,1)),[],1)-maxbounds(:,2);
else
    rmax=NaN(0,1);
end
r=[rmin;rmax];
if nargout==1
    return
end
nmin=length(rmin);
nmax=length(rmax);
J=repmat(p_axpy(0,p,[]),nmin+nmax,1);
for i=1:nmin
    J(i).parameter(minbounds(i,1))=-1;
end
for i=1:nmax
    J(nmin+i).parameter(maxbounds(i,1))=1;
end
end
