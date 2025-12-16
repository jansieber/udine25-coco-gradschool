%% Demo AB reaction
% x(1)=concentration fraction
% x(2)=temperature
% parameters [alpha,beta,gamma]=reaction rate, exothermicity, cooling
%% set path
clear;
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
%% Define parameters and variables
pnames={'alpha','beta','gamma'};
xnames={'c','T'}; % concentration and temperature
[iv,ip]=deal(structind_from_names(xnames),structind_from_names(pnames));
%% define right-hand side
ab=@(c,T,alpha,beta,gamma)[...     % ODE r.h.s.
    -c+alpha.*(1-c).*exp(T);...
    -T+beta.*alpha.*(1-c).*exp(T)-gamma.*T];
rhs=@(x,p)ab(x(iv.c,:),x(iv.T,:),p(ip.alpha,:),p(ip.beta,:),p(ip.gamma,:));
%% initial guess for equilibrium
par0([ip.alpha,ip.beta,ip.gamma])=...
     [     0;       14;     2];
par0=par0(:);
x0=[0;0];
