%% test bc's generated
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
[ip,npars]=structind_from_names(pnames);
% initial parameters
um_pars([ip.r, ip.cm, ip.beta, ip.omega, ip.cp, ip.a, ip.phi])=...
    [      1,    0,     0.5,       2*pi,     4,     1,    0];
um_pars=um_pars(:);
% initialize (x,p,q) such that the satisfy all boundary conditions
[                   a0,            phi0,          cm0]=...
    deal(um_pars(ip.a),um_pars(ip.phi),um_pars(ip.cm));
um_x=[-sqrt(a0)-cm0;cos(phi0)/sqrt(2);sin(phi0)/sqrt(2)];
bdprep=@(fname)cellfun(@(f){@(varargin)f(varargin{2:end})},...
    arrayfun(@(i){sco_gen(fname,i)},0:2));
bc_m=bdprep(@r_tipping_bc_minus);
[t0,T]=deal(0,1);
bc_m{1}(t0,T,um_x,um_x,um_pars)
bc_m{2}(t0,T,um_x,um_x,um_pars)
bc_m{3}(t0,T,um_x,um_x,um_pars)
