function funcs=set_symfuncs(fcnrhs,varargin)
%% convert output of symbolic toolbox into funcs structure recognized by DDE-Biftool
%
%% Inputs
%
% |fcnrhs|: function handle of function created by dde_sym2funcs(...)
%
%% Optional inputs (name value pairs)
%
% * |'sys_tau'|: some inputs are passed on to |set_funcs|. In particular,
% |sys_tau| may need to be specified, if the delays are parameters.
% * |'x_vectorized'| (|true|): if resulting functions can be called with 3d
% array of |xx| arguments [dimension x (number of delays +1) x nvec].
% Typically, the symbolic toolbox generated functionsa permit this.
% * |'p_vectorized'| (|true|): can all the nvec different |xx| inputs also
% have different parameter values? If yes, then the function gets called
% with a parameter array of size 1 x np x nvec. 
% * |'lhs_matrix'|: specify left-hand side matrix for DDAEs
%
%% Output
%
% * |funcs|: structure containing system r.h.s, delays, derivatives and
% other information.
%%
default={'x_vectorized',true,'p_vectorized',true};%,'extrafuncs',[]}; % change default on vectorization
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ischar(fcnrhs)
    f=str2func(fcnrhs);
else
    f=fcnrhs;
end
directional_derivative=f('directional_derivative');
sys_rhs=dde_sym2fun(f,'rhs');
sys_dirderi=dde_sym2fun(f,'dirderi');
sys_mfderi={};
if ~directional_derivative
    sys_mfderi=dde_sym2fun(f,'mfderi');
end
xpattern=reshape(f('xpattern'),2,[]);
fcn_args={'wrap_rhs',sys_rhs,... %'sys_deri',sys_deri,
    'wrap_dirderi',sys_dirderi,'sys_mfderi',sys_mfderi,'xpattern',xpattern};
%% construct functions specifying variable delay if requested
sd_delay=f('tp_del');
ntau=f('ntau');
if sd_delay
    assert(directional_derivative);
    sys_ntau=@()ntau;
    sys_tau=dde_sym2fun(f,'tau');
    sys_dirdtau=dde_sym2fun(f,'dirdtau');
    tau_args={'wrap_tau',sys_tau,...'%sys_dtau',sys_dtau,
        'sys_ntau',sys_ntau,...
        'wrap_dirdtau',sys_dirdtau,'sys_tau_seq',f('sys_tau_seq')};
    fcn_args={fcn_args{:},tau_args{:}}; %#ok<CCAT>
end
funcs=set_funcs(fcn_args{:},...
    'x_vectorized',options.x_vectorized,'p_vectorized',options.p_vectorized,... %'format',format,...
    pass_on{:});
end
%%
function y=wrap_mfderi(f,xx,p,devs)
order=length(devs);
np=length(p);
[n,ntau]=size(xx);
dx=zeros(n,ntau,order);
dp=zeros(np,order);
for i=1:order
   devi=reshape(devs{i},n+np,ntau);
   dx(:,:,i)=devi(1:n,:);
   dp(:,i)=devi(n+1:end,1);
end
xx=xx(1:n,:);
y=dde_sym_rhs_wrap(f,order,xx,p,dx,dp);
end
