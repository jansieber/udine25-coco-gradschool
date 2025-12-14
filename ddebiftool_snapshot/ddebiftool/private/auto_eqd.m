function eqf=auto_eqd(point,varargin)
%% find cumulative discretization error on given mesh, formula from AUTO
% function [eqf]=auto_eqd(point,...)
%% INPUT:
%   * point with full mesh, degree, profile
%   * extension: function ext_err=fcn(intlen,err) where intlen is a ntst vector
%   of interval lengths and err is a ntst-1 vector of error estimates. the
%   function should return an array ext_err of length ntst+1 that estimates
%   the error also to the left and right boundary. e.g
%
% @(ms,err)[err(:,end),err,err(:,1)] for periodic functions
% @(ms,err)[err(:,1),err,err(:,end)] for connecting orbits
%  
%% OUTPUT:
%       eqf values of monotonically increasing function eqdf
%% COMMENT: 
%       this function is a matlab translation of the AUTO
%       fortran function EQDF, except
%	- it leaves the extrapolation to the caller (two examples are prepared
%	for periodic and connecting orbits) with the function extension
%       - it uses different ordering of the data in ups
%	- the integral is approximated using the trapezium rule

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000
%
% $Id: auto_eqd.m 308 2018-10-28 15:08:12Z jansieber $
%
%% find out dimensions of solution array
%	ntst number of intervals
%   ndim system dimension
%	ncol number of collocation points (per interval)
default={'new_degree',point.degree,'min_err',1e-7};
options=dde_set_options(default,varargin,'pass_on');
tcoarse=point.mesh(1:point.degree:end);
%% if extension beyondf boundary is defined (eg for periodic orbits)
% extend
eqf=0:numel(point.mesh(1:point.degree:end))-1;
if exist(['dde_',point.kind,'_extend_profile'],'file')
    t=[tcoarse(1)-tcoarse(end)+tcoarse(end-1),tcoarse,tcoarse(end)+tcoarse(2)-tcoarse(1)];
    [point.mesh,point.profile]=dde_apply({'dde_',point.kind,'_extend_profile'},point,t);
    tcoarse=point.mesh(1:point.degree:end);
    extended=true;
else
    extended=false;
end
dtm=diff(tcoarse);
tmid=(tcoarse(2:end)+tcoarse(1:end-1))/2;
%% calculate highest derivative of collocation polynomials
deg=min(point.degree,options.new_degree);
hd=dde_coll_eva(point,tmid,'diff',deg);
%% Take care of "small derivative" case:
if max(abs(hd(:)))<options.min_err
    return
end
%% extrapolate to get one more interval
if ~extended
    [dtm,hd]=dde_coll_extend_profile(dtm,hd);
end
%% compute approximation to (ncol+1)-st derivative, 
% by taking divided differences of neighboring hd values
scav=2./(dtm(1:end-1)+dtm(2:end));
scav=repmat(scav,size(point.profile,1),1);
hdp1=diff(hd,[],2).*scav;

%% define the cumulative error err_dt, to be equidistributed
err=sum(abs(hdp1).^(1/(point.degree+1)),1);
err_dt=0.5*(err(1:end-1)+err(2:end)).*dtm(2:end-1);
%% avoiding zero error change (for interpolation)
err_dt=max(err_dt,options.min_err*max(err_dt));
eqf=[0,cumsum(err_dt)];
end
function [dtext,err_ext]=dde_coll_extend_profile(dt,err_est)
%% function extending highest derivative by linear extrapolation  
% (used in auto_eqd)
%
% Input: dt: time increments, hd: error estimate
% Outputs: dtext,err_ext: extended time increments and error estimates
%
% $Id: dde_coll_extend_profile.m 308 2018-10-28 15:08:12Z jansieber $
%%
dtext=[2*dt(1)-dt(2),dt,2*dt(end)-dt(end-1)];
err_ext=[2*err_est(:,1)-err_est(:,2),err_est,2*err_est(:,end)-err_est(:,end-1)];
end
