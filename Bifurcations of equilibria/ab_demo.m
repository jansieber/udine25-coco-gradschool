%% Demo AB reaction
% x(1)=concentration fraction
% x(2)=temperature
% parameters [alpha,beta,gamma]=reaction rate, exothermicity, cooling
clear;
format compact
reaction=@(x,p)p(1,:).*(1-x(1,:)).*exp(x(2,:));
ab=@(x,p)[...     % ODE r.h.s.
    -x(1,:)+reaction(x,p);...
    -(1+p(3,:)).*x(2,:)+p(2,:).*reaction(x,p)];
p0=[0;14;2];
x0=[0;0];
