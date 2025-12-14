function varargout=dde_psol_delay_zero_cond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_psol_delay_zero_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_psol_delay_zero_cond(p,varargin{2:end}));
end
end
function [r,Jpt]=loc_psol_delay_zero_cond(point,funcs,free_par_ind,d_nr,itz)
%% residual & derivative wrt x, T and parameter for tau_d_nr(tz)=tau'_dnr(tz)=0
% formulated in sys_cond format
%
% inputs
% funcs: system functions
% point: point (type psol) for which delay=0 is to be
% determined
% tz: time for psol at which delay is supposed to be zero (for psol)
% d_nr: delay number
% free_par_ind: indices of free parameters
%
% outputs
% r residual: tau for stst,hopf,fold (1d), tau(tz) and tau'(tz) for psol
% J Jacobian (struct of type point): single row for stst,hopf,fold, two
% rows for psol
%
cond=@(xx,p)funcs.wrap_tau(d_nr,xx,p);
get_nr=@(y)y(d_nr+1,:);
dcond=@(xx,p,dxx,dp)get_nr(funcs.dtau_dir(1,xx,p,dxx,dp));
[r,Jpt]=dde_nlin_extreme_cond(point,funcs,cond,'dcond',dcond,...
    'free_par',free_par_ind,'itz',itz,'ntau',d_nr);
end
