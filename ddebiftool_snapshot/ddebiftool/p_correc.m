function [newpoint,success]=p_correc(funcs,point,free_par,step_cnd,method,remesh_flag,...
    previous,varargin)
%% correct point using Newton iteration
% function [point,success]=p_correc(point0,free_par,step_cnd,method,adapt,
%                                 previous)
% INPUT:
%   funcs problem functions
%   point0 initial point guess
%   free_par free parameter numbers in N^d
%       step_cnd steplength condition(s) as point(s)
%   method method parameters
%   remeshflag if zero or absent, do not adapt mesh; if one, always adapt
%       previous (optional) previously computed branch point (used, in case
%            of periodic solutions or connecting orbits, to
%            minimize phase shift)
% OUTPUT:
%       point solution
%       success nonzero for successfull computation
%
% (c) DDE-BIFTOOL v. 2.02, 16/6/2002
%
% $Id: p_correc.m 348 2019-06-19 13:09:01Z jansieber $
%
%%
%% optional parameters:
if nargin<=5
    remesh_flag=0;
end
if nargin<=6
    previous=[];
end
[f,x0,data]=p_correc_setup(funcs,point,free_par,method,...
    'remesh_flag',remesh_flag,'previous',previous,'step_cnd',step_cnd,varargin{:});
%% perform Newton-Raphson iterations:
[x,success]=dde_nsolve(f,x0,data.method,'rJini',data.rJini);
data.point=dde_point_from_x(x,data.point,data.free_par);
newpoint=data.point;
%% store extracolumns and artificial parameters if present
if isfield(data,'extracolumns') && size(data.extracolumns,2)>0
    newpoint.nvec.extracolumns=data.extracolumns;
    newpoint.nvec.extrapars=x(end-size(data.extracolumns,2)+1:end);
end
newpoint=dde_postprocess(newpoint,data);
%% recorrect with adapted mesh if necessary
if ~success || ~isfield(method,'remesh') || ~method.remesh || ...
        ~isfield(newpoint,'mesh') || isempty(newpoint.mesh)
    return
end
ma=method.adapt_mesh_after_correct;
% do not adapt mesh when p_nr=0
% adapt mesh when p_nr=1
% adapt mesh when p_nr>1 & mod(p_nr,ma)=0
if ~ (remesh_flag==1 || (remesh_flag>1 && mod(remesh_flag,ma)==0 ))
    return
end
method2=method;
method2.adapt_mesh_before_correct=1;
method2.adapt_mesh_after_correct=0;
[newpoint_remeshed,success]=p_correc(funcs,newpoint,free_par,step_cnd,method2,2,previous,...
    varargin{:});
if success
    newpoint=newpoint_remeshed;
end
end
