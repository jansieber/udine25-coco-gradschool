%% use symbolic toolbox to generate r.h.s. and derivatives
clear
if sco_isoctave()
    pkg load symbolic
end
parnames={'nu','be','ga','r','a3','b3'};
%% Create symbols for parameters, states
syms(parnames{:});       % create symbols for parameters
par=cell2sym(parnames);  % now tau is par(1) etc
x=sym('x',[3,1]);
y= [( -(be+nu)*x(1) + be*x(2) - a3*x(1)^3 + b3*(x(2)-x(1))^3 )/r;
        be*x(1) - (be+ga)*x(2) - x(3) - b3*(x(2)-x(1))^3;
        x(2)];

F=sco_sym2funcs(...
    y,...
    {x,par(:)},...
    {'x','p'},...
    'filename','sym_tor');