%% residual & derivative wrt x, T and parameter for tau_d_nr=0 for equilibria
% formulated in sys_cond format
%
% inputs
% funcs: system functions
% point: point (type psol,stst,hopf,fold) for which delay=0 is to be
% determined
% tz: 
% d_nr: delay number
% free_par_ind: indices of free parameters
%
% outputs
% r residual: tau for stst,hopf,fold (1d), tau(tz) and tau'(tz) for psol
% J Jacobian (struct of type point): single row for stst,hopf,fold, two
% rows for psol
%
%%
function varargout=dde_stst_delay_zero_cond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_stst_delay_zero_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_stst_delay_zero_cond(p,varargin{2:end}));
end
end
function [r,Jp]=loc_stst_delay_zero_cond(point,funcs,free_par_ind,d_nr,varargin)
ntau=dde_num_delays(funcs);
xx=point.x(:,ones(ntau+1,1));
r=dde_stst_delays(funcs,point,'d_nr',d_nr,'repeat',false);
Jp=p_axpy(0,point,[]);
Jxall=funcs.dtau_mf(xx,point.parameter,{1,'I'},0);
Jp.x=reshape(sum(Jxall(d_nr,:,1:d_nr), 3),size(Jp.x));
Jpall=funcs.dtau_mf(xx,point.parameter,0,{1,{free_par_ind}});
Jp.parameter(free_par_ind)=Jpall(d_nr,:);
end
