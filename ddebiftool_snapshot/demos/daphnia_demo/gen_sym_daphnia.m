%% Daphnia model equations
%%
clear
addpath('../../ddebiftool',...
    '../../ddebiftool_extra_symbolic');
%% Set number of delays and parameter names
parnames={'szB','szA','szmax','growth','sigma','feed','rep','mu','cap',...
    'flow','amax','amaxrel'};
cs=[parnames;num2cell(1:length(parnames))];
ip=struct(cs{:});
%% Create symbols for parameters, states and delayed states
% The array |par| is the array of symbols in the same order as parnames.
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now beta is par(1) etc
dr0=@(r)flow*r*(1-r/cap); % resource rate of change
fr=@(r)sigma*r/(1+sigma*r); % Holling type II resource dependence
cgrow=@(r,s)growth*(szmax*fr(r)-s); % consumer growth rate
ceff=@(r,s)fr(r)*s^2; % consumer effect per individual
%% consumer size for consumers of age a at time t
% is given by s(a,t)=exp(-growth a)*szB+sd(a,t) where sd(a,t) is the purely
% cumulative part of the size
%
% sd(t,a)=exp(-growth a)*szB+\int_0^a growth szmax exp(-growth al)f_r(r(t-al) dal
%% reproduction rate per individual at age a is
% rep*s(t,a)^2 f_r(r(t))
% this gets multiplied by the number of individuals of age a
% exp(-mu a)b(t-a), which is the integrand:
% rep*s(a,t)^2*exp(-mu*a)*fr(r(t))*b(t-a)
%
% so
%
% c_k(t)=\int_0^a_k rep*s(a,t)^2*exp(-mu*a)*fr(r(t))*b(t-a) da
%
% for various a_k (results c_k(t)). The age-independent term fr(r(t)) can
% be factored out and included into the c_k(t):
%
% Defining the discounted birth rate bd(t)=b(t)/f_r(r(t)), we get that
%
% bd(t)=\int_aA^amax rep*s(a,t)^2*exp(-mu*a)*fr(r(t-a))*bd(t-a) da
%
% and consumption is
%
% c(t)=\int_0^amax feed*s(a,t)^2*exp(-mu*a)*fr(r(t))*b(t-a) da
%
% Again, extracting factor f_r(r(t)) and incorporating it into the discounted
% consumption cd0(t)=c(t)/f_r(r(t))
%
% cd0(t)=\int_0^amax feed*s(a,t)^2*exp(-mu*a)*fr(r(t-a))*bd(t-a) da
%
% where the resource equation is  
%
% r'(t)=fres(r(t))-f_r(r(t)) cd0(t). Exploiting that consumption per age
% cohort is a multiple of the reproduction, we compute only one integral,
% the consumer effect
%
% ef(t,al)=\int_0^{al*amax} rep*s(a,t)^2*exp(-mu*a)*fr(r(t-a))*bd(t-a) da
%
% and then set bd(t)=ef(t,1)-ef(t,raA),
% cd(t)=ef(t,1)
%
% We also have discounted size at maturation age, given by
% sdA(t)=sd(t,raA(t)*amax).
%
% which satisfies the algebraic equation
%
% 0=sdA(t)+exp(growth*aA(t))-szA
%
%% Variables and delays
% we have two delays, aA and amax
%
% 5 dynamic variables: resource r(t), discounted birth rate bd(t),
% discounted consumption rate cd(t), maturation age aA(t), discounted
% maturation size sdA(t).
% They satisfy 1 DE, 1 AE and 3 equations defined by the integral:
ntau=0;
ntp1=ntau+1;
xnames={'r','raA','bd','cd','sdA'};
csx=[xnames(:)';num2cell(1:length(xnames))];
ix=struct(csx{:});
r=  sym('r',  [1,ntp1]); % resource
raA=sym('raA',[1,ntp1]); % maturation age relative to amax
bd= sym('bd', [1,ntp1]); % discounted birth rate
cd= sym('cd', [1,ntp1]); % discounted consumption
sdA=sym('sdA',[1,ntp1]); % discounted maturation size
%% Equations
dr=   dr0(r(1))-(feed/rep)*fr(r(1))*cd(1);    % r.h.s. DE for resource
ageth=sdA(1)+exp(-growth*amax*raA(1))*szB-szA; % implicit definition of maturation age
%% DDE
dde_sym2funcs(...
    [dr;ageth],...       % nf x 1 array of derivative symbolic expressions
    [r;raA;bd;cd;sdA],... % nx x (ntau+1) array of symbols for states (current & delayed)
    par,...              % 1 x np (or np x 1) array of symbols used for parameters
    'sys_tau',sym([]),...
    'filename','sym_daphnia_dde',...
    'directional_derivative',true,...
    'maxorder',3,'print_progress',true);
%% Integrands: 1st integrand for growth growth of individual
syms ra a sd bda
sd_int=exp(-growth*a)*cgrow(ra,0);
[dum1,ips]=intersect(par,symvar(sd_int));
ip_sd=sort(ips);
x_sd=[a;ra];
x_sd_names=dde_names_from_sym(x_sd);
cs_sd=[x_sd_names(:)';num2cell(1:length(x_sd_names))];
ix_sd=struct(cs_sd{:});
dde_sym2funcs(...
    sd_int,... % n x 1 array of symbolic expressions
    x_sd,...   % n or n+1 array of symbols for arguments
    par(ip_sd),...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_sd_int',...
    'directional_derivative',true,...
    'maxorder',3);
%% Integrands: 2nd integrand for effect off population
sa=sd+exp(-growth*a)*szB;
ceff_int=simplify(rep*sa^2*exp(-mu*a)*fr(ra)*bda);
[dum2,ipc]=intersect(par,symvar(ceff_int));
ip_ceff=sort(ipc);
x_ceff=[a;ra;bda;sd];
x_ceff_names=dde_names_from_sym(x_ceff);
cs_ceff=[x_ceff_names(:)';num2cell(1:length(x_ceff_names))];
ix_ceff=struct(cs_ceff{:});
dde_sym2funcs(...
    ceff_int,...      % n x 1 array of symbolic expressions
    x_ceff,...        % n or n+1 array of symbols for arguments
    par(ip_ceff),...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_ceff_int',...
    'directional_derivative',true,...
    'maxorder',3);
%% formulas for equilibrium
syms r0 raA0 bd0 cd0 sdA0
xeq=[r0,raA0,bd0,cd0,sdA0];
eqsubs={[r(1),raA(1),bd(1),cd(1),sdA(1)],xeq};
dreq=subs(dr,eqsubs{:});
cdeq=solve(dreq,cd0);
sdAeq1=solve(subs(ageth,eqsubs{:}),sdA0);
syms al
sdeqint=int(subs(sd_int,[a,ra],[al,r0]),al,0,a)+szB;
sdAeq2=subs(sdeqint,a,amax*raA0);
raAeq=simplify(solve(sdAeq1-sdAeq2,raA0));
ceff_eqintg=simplify(subs(ceff_int,[bda,ra,sd],[bd0,r0,sdeqint]));
ceff_eqint=simplify(int(subs(ceff_eqintg/bd0,a,al),al,0,a));
ceff_eqintmax=simplify(subs(ceff_eqint,a,amax));
bdeq=simplify(solve(ceff_eqintmax*bd0-cdeq,bd0));
sdAeq=subs(sdAeq1,raA0,raAeq);
stst_eq=simplify(subs(ceff_eqint,a,amax)-subs(ceff_eqint,a,raAeq*amax));
dde_sym2funcs(...
    stst_eq,...       % n x 1 array of symbolic expressions
    r0,...        % n or n+1 array of symbols for arguments
    par,...          % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_ststeq',...
    'directional_derivative',true,...
    'maxorder',1);
dde_sym2funcs(...
    [raAeq;bdeq;cdeq;sdAeq],...       % n x 1 array of symbolic expressions
    r0,...        % n or n+1 array of symbols for arguments
    par,...          % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_stst',...
    'directional_derivative',true,...
    'maxorder',1);


%%
save('index_daphnia','ip','ix','ip_sd','ix_sd','ip_ceff','ix_ceff');