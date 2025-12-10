%% Generate right-hand side for rate-induced tipping problem
% reference Thoraya Alharti PhD thesis
% <https://ore.exeter.ac.uk/ndownloader/files/56835545> 
%%
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'))
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
syms(pnames{:})
pars=cell2sym(pnames);
syms p q x
%% right-hand side
L=cm+(p^2+q^2)*(cp-cm)+beta*q*sqrt(p^2+q^2);
f=[(x+L)^2-a;...
    [r-r*(p^2+q^2), omega; -omega, r-r*(p^2+q^2)]*[p;q]];
sco_sym2funcs(f,{[x;p;q],pars},{'x','par'},'filename','r_tipping_rhs.m');
%% boundary condition for u_minus
syms t0 T p0 q0 x0 p1 q1 x1
bc_minus=[(x0+cm)^2-a; [p1;q1]-[cos(phi);sin(phi)]/sqrt(2)];
sco_sym2funcs(bc_minus,{T,[x0;p0;q0],[x1;p1;q1],pars},...
    {'T','x0','x1','par'},'vector',[0,1,1,1],'filename','r_tipping_bc_minus.m');
%% boundary condition for u_gamma
% 4 b.c., assigning phase of periodic forcing to phi
bc_gamma=[x0-x1;p0-p1;q0-q1;p1*sin(phi)-q1*cos(phi)];
sco_sym2funcs(bc_gamma,{T,[x0;p0;q0],[x1;p1;q1],pars},...
    {'T','x0','x1','par'},'vector',[0,1,1,1],'filename','r_tipping_bc_gamma.m');
%% boundary condition for u_plus
% only two b.c. are purely on the u_+ segment
bc_plus=[p0-cos(phi)/sqrt(2);q0-sin(phi)/sqrt(2)];
sco_sym2funcs(bc_plus,{T,[x0;p0;q0],[x1;p1;q1],pars},...
    {'T','x0','x1','par'},'vector',[0,1,1,1],'filename','r_tipping_bc_plus.m');
