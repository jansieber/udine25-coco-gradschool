function varargout=dde_pd_cond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_pd_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_pd_cond(p,varargin{2:end}));
end
end
function [r,J]=loc_pd_cond(p,xdim,t0)
Id=eye(xdim);
Z=zeros(xdim);
for i=length(t0):-1:1
    trafo=cat(1,[(1+cos(pi*t0(i)))*Id,(1-sin(pi*t0(i)))*Id],...
        [Z,Z]);
    [rc{i},Jc{i}]=dde_psol_lincond(p,xdim,'profile','trafo',trafo,'shift',[0,1],...
        'condprojint',t0(i)*[1,1],'condprojmat',[Id,Id],'stateproj',[Z,Id,Z;Z,Z,Id]);
end
r=cat(1,rc{:});
J=cat(1,Jc{:});
end
