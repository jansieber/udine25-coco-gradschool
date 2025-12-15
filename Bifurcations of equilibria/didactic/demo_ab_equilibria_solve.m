%% Find equilibria for AB reaction using a sequence of Newton iterations
clear;
addpath(fullfile('didactic_tools/'));
beta=14;
gamma=2;
ab=@(x,alpha)[...     % ODE r.h.s.
    -x(1)+alpha*(1-x(1))*exp(x(2));...
    -x(2)+beta*alpha*(1-x(1))*exp(x(2))-gamma*x(2)];
alpha=0.05;           % initial choice for alpha   
x0=[0;0];             % initial condition
n=length(x0);         % dimension of problem
f=@(x)ab(x,alpha);
%% arguments for Newton iteration
x_star=Solve(f,x0,'print',1); % find equilibrium with Newton iteration
%% sweep, trying t ofind equilibria for all  alpha
alpha_rg=0:0.0025:0.15;
na=length(alpha_rg);
xeq=NaN(n,na);
fval=xeq;
converged=false(1,na);
detJ=NaN(1,na);
xguess=x0;
for i=1:na
    alpha=alpha_rg(i);
    f=@(x)ab(x,alpha);
    [xeq(:,i),converged(i),J]=Solve(f,xguess,'print',1);
    %[xeq(:,i),fval(:,i),converged(i)]=fsolve(f,xguess,optimset('Display','off'));
    %J=Jacobian(f,xeq(:,i),1e-4);
    detJ(i)=det(J);
    fprintf('i=%2d, alpha=%5.3g, converged=%d, det(J)=%5.3g, x=(%6.4g,%6.4g)\n',...
        i,alpha,converged(i),detJ(i),xeq(1,i),xeq(2,i));
    if converged(i)
        xguess=xeq(:,i);
    end
end
%% plot
clf;
lw={'linewidth',2};
plot(alpha_rg(converged),xeq(1,converged),'o',...
    alpha_rg(~converged),xeq(1,~converged),'rx',lw{:});
set(gca,'fontweight','bold','ylim',[0,1]);
xlabel('alpha');
ylabel('concentration');
grid on
