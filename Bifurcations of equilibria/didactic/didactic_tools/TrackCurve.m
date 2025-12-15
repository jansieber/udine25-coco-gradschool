%% Track curve of roots using pseudo arclength continuation
%
% Equation is
%
% $$ 0=f(y) $$
%
% with initial tangent |ytan|, and initial solution |y0|.
%%
function [ylist,ytan]=TrackCurve(userf,y0,ytan,varargin)
%% Initialization
defaults={'stepsize',0.01,'nmax',100,'print',1,'stop',@(y)(0),...
    'userdf',[],'hjac',1e-4,'maxstep',0.1,'minstep',1e-4,'stepinc',1.1};
[options,pass_on]=MySetOptions(defaults,varargin,'pass_on');
userargs=[pass_on,{'print',options.print-1}];
if isempty(options.userdf)
    options.userdf=@(y)Jacobian(userf,y,options.hjac);
end
dim=length(y0);
ylist=zeros(dim,options.nmax);
ytan0=[zeros(dim-1,1);1];
%% set up extended functions
fext=@(y,ytan,ypred)[userf(y);ytan'*(y-ypred)];
dfext=@(y,ytan)[options.userdf(y);ytan'];
%% correct initial guess
f=@(y)fext(y,ytan,y0);
df=@(y)dfext(y,ytan);
[y,conv,J]=Solve(f,y0,'df',df,userargs{:});
if conv
    npoints=1;
    ylist(:,npoints)=y;
    ytanold=ytan;
    ytan=J\ytan0;
    ytan=ytan/norm(ytan,'inf');
    ytan=ytan.*sign(ytanold'*ytan);
else
    fprintf('initial guess too far off or problem ill-posed\n');
    ylist=[];
    return;
end
%% loop to find branch of solutions
stepsize=options.stepsize;
for i=2:options.nmax
    yold=ylist(:,npoints);
    ypred=yold+stepsize*ytan;
    f=@(y)fext(y,ytan,ypred);
    df=@(y)dfext(y,ytan);
    [y,conv,J]=Solve(f,ypred,'df',df,userargs{:});
    if conv
        npoints=npoints+1;
        ylist(:,npoints)=y;
        ytanold=ytan;
        ytan=J\ytan0;
        ytan=ytan/norm(ytan,'inf');
        ytan=ytan.*sign(ytanold'*ytan);
        if options.print>0 
            fprintf('step %d, stepsize=%g, p=%g\n',i,stepsize,y(end));
        end
        stepsize=min(stepsize*options.stepinc,options.maxstep);
    elseif ~conv && stepsize>options.minstep
        stepsize=stepsize/2;
        if options.print>0
            fprintf('convergence failed at step %d, new stepsize=%g\n',...
                i,stepsize);
        end
    else
        if options.print>0
            fprintf('convergence failed with minimal stepsize, leaving\n');
        end
        break
    end
    if options.stop(ylist(:,npoints))
        break
    end        
end
ylist=ylist(:,1:npoints);
end
