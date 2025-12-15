%% Solve ODE with Runge-Kutta method
% Equation is 
%
% $$ x'=f(t,x) $$ 
%
% on |[Tspan(1),Tspan(2)]|, initial value |x(Tspan(1))=x0|, solution
% makes |N| steps. This function is vectorised (not necessary for
% students).
%%
function [x,t,xtraj]=MyIVP(f,x0,Tspan,N,varargin)
%% optional arguments
default={'stop',@(x)false};
options=MySetOptions(default,varargin);
x=x0;
%% pre-allocate for trajectory if output is requested
if nargout>2
    xtraj=repmat(x,[1,1,N+1]);
end
h=(Tspan(2)-Tspan(1))/N;
t=linspace(Tspan(1),Tspan(2),N+1);
%% make N steps using DOPRI formula
% stop if optional stopping function returns true (useful for ODEs that
% blow up in finite time). 
for i=1:N
    x=dopri45(f,t(i),x,h);
    if nargout>2
        xtraj(:,:,i+1)=x;
    end
    if options.stop(x)
        t=t(1:i);
        if nargout>2
            xtraj=xtraj(:,:,1:i+1);
        end
        break
    end
end
if nargout>2
    xtraj=permute(xtraj,[1,3,2]);
end
end
%% Dormand & Prince formula (from ode45)
function x=dopri45(f,t,x1,h)
A = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];
tvec=t+A*h;
tvec=tvec(:,ones(1,size(x1,2)));
B = h*[
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    0           0       0       0               0               0
    ];
f1=f(tvec(1,:),x1);
x2=x1+B(1,1)*f1;
f2=f(tvec(2,:),x2);
x3=x1+B(1,2)*f1+B(2,2)*f2;
f3=f(tvec(3,:),x3);
x4=x1+B(1,3)*f1+B(2,3)*f2+B(3,3)*f3;
f4=f(tvec(4,:),x4);
x5=x1+B(1,4)*f1+B(2,4)*f2+B(3,4)*f3+B(4,4)*f4;
f5=f(tvec(5,:),x5);
x6=x1+B(1,5)*f1+B(2,5)*f2+B(3,5)*f3+B(4,5)*f4+B(5,5)*f5;
f6=f(tvec(6,:),x6);
x=x1+B(1,6)*f1+B(2,6)*f2+B(3,6)*f3+B(4,6)*f4+B(5,6)*f5+B(6,6)*f6;
end

%% classical Runge-Kutta formula (see Wikipedia)
% function x=rk34(f,x0,h)
% x=x0;
% k1=h*f(x);
% x=x0+0.5*k1;
% k2=h*f(x);
% x=x0+0.5*k2;
% k3=h*f(x);
% x=x0+k3;
% k4=h*f(x);
% x=x0+(k1+2*k2+2*k3+k4)/6;
% end
