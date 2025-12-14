function varargout=sys_cond_coll_fixperiod(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_cond_coll_fixperiod(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_cond_coll_fixperiod(p,varargin{2:end}));
end
end
%%
function [r,J]=loc_cond_coll_fixperiod(point,period,ispar)
%% fix period to given parameter (visible to user fcns)
% implemented as user condition
%% 
if nargin<3
    ispar=true;
end
J=p_axpy(0,point,[]);
J.period=-1;
if ispar
    ind_period=period;
    period=point.parameter(ind_period);
    J.parameter(ind_period)=1;
end
r=period-point.period;
end
