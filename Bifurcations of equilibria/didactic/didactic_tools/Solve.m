function [x,converged,J]=Solve(f,x0,varargin)
defaults={'df',[],'hdev',1e-5,'maxit',10,'tolerance',1e-8,'print',1,'damping',1};
options=MySetOptions(defaults,varargin);
if isempty(options.df)
    options.df=@(x)Jacobian(f,x,options.hdev);
end
x=x0;
converged=0;
for i=1:options.maxit
    y=f(x);
    J=options.df(x);
    cor=-J\y;
    ynorm=norm(y,'inf');
    cornorm=norm(cor,'inf');
    x=x+options.damping*cor;
    if options.print>0
        fprintf('it=%d, |cor|=%g, |res|=%g\n',i,cornorm,ynorm);
    end
    if max(cornorm,ynorm)<options.tolerance
        converged=1;
        break
    end
end
end
