function [ip,freepar]=dde_ip_addnullpar(ip,point,nullparind,extend,varargin)
default={'lastname','nuserpar','contpar',[],'periodnullpar',zeros(1,0)};
options=dde_set_options(default,varargin,'pass_on');
if ~isfield(ip,'nuserpar')
    ip.nuserpar=length(point.parameter); % number of original system parameters
end
nnull=length(nullparind)*extend;
ip.nullparind=NaN(nnull,2); % location of parameters for nullspace vector
ip.nullparind(:,1)=nullparind;
ip.nullparind(:,1+(1:extend))=ip.(options.lastname)+reshape(1:nnull,[],extend);
ip.periodnullpar=options.periodnullpar;
freepar=[options.contpar,reshape(ip.nullparind(:,2:end),1,[])];
end
