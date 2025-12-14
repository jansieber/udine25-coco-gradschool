%% Iterate toward points with montor change along branch of periodic orbits
%% Input
%
% * funcs: problem functions
% * branch: solution branch
%
% * optional name-value pairs:
% * 'monitor' (default stability) function [val,p_update]=monitor(funcs,br,p,pref,...) returns
% discrete or continuous monitoring value val, with updated point p to
% p_update.
% * 'range' (default=|1:length(branch.point)|, which part of branch
% * 'type' (default 'sign'='continuous', alternative 'discrete') what kind
% of monitor to be detected
% * all other name-value pairs are passed on as fields for branch fields or
% method fields
%% Output
% * branch: updated branch with inserted special points
% * testrg: monitor values
% * bifpoints: 1 x ns cell array of points close to stability change
% * indices: 1 x ns array of pointers into branch.point at which stability
% change points are inserted
%%
function varargout=MonitorChange(funcs,branch,varargin)
default={'range',1:length(branch.point),'printlevel',1,...
    'type','discrete','output','branch','refine',true,'monitor',...
    @(f,b,p,pref,varargin)br_stability_monitor(f,b,p,varargin{:})};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch=replace_branch_pars(branch,[],pass_on);
ptrg=branch.point(options.range);
ptbefore=branch.point(1:options.range(1)-1);
ptafter=branch.point(options.range(end)+1:end);
br=setfield(branch,'point',ptrg); %#ok<SFLD>
flags=zeros(1,length(br.point));
getmon=@(p,pref)options.monitor(funcs,br,p,pref,pass_on{:});
ind=0;
changes=[];
while true
    [testrg,br]=br_check(br,getmon);
    changes=find_changes(testrg,options.type);
    unrefined=changes(flags(changes)==0);
    if isempty(unrefined) || ~options.refine
        break
    end
    curind=unrefined(1)+(0:1);
    curmon=testrg(curind);
    bif_mon=@(p,pref)iterate_mon(getmon,p,pref,curmon,options.type);
    [br,ibif,imap]=br_bisection(funcs,br,curind,bif_mon,...
        'bif_mon_reference',true,'printlevel',options.printlevel-1,pass_on{:});
    if isnan(ibif)
        flags(curind)=1;
    else
        flags(imap)=flags;
        flags(ibif+(-1:0))=1;
    end
    ind=ind+1;
    if options.printlevel>0
        if ~isnan(ibif)
            parbif=br.point(ibif).parameter(br.parameter.free);
            fprintf('Bif %d of %d at parameters (%s)\n',ind,length(changes),...
                num2str(parbif));
        else
            fprintf('Bif %d of %d, iteration unsuccessful\n',ind,length(changes));
        end
    end
end
indices=changes+length(ptbefore);
branchout=setfield(br,'point',[ptbefore,br.point,ptafter]); %#ok<SFLD>
bifpoints=branchout.point(indices);
testrg=[NaN(1,length(ptbefore)),testrg,NaN(1,length(ptafter))];
out=cell(1,4);
switch options.output
    case 'bifurcations'
        [out{:}]=deal(bifpoints,indices,branchout,testrg);
    case 'bifurcation_indices'
        [out{:}]=deal(indices,branchout,testrg,bifpoints);
    case 'tests'
        [out{:}]=deal(testrg,branchout,bifpoints,indices);
    case 'branch'
        [out{:}]=deal(branchout,testrg,bifpoints,indices);        
end
varargout=out;
end
%%
function [testrg,br]=br_check(br,mon)
testrg=NaN(1,length(br.point));
[testrg(1),br.point(1)]=mon(br.point(1),repmat(br.point(1),1,0));
for i=2:length(br.point)
    [testrg(i),br.point(i)]=mon(br.point(i),br.point(i-1));
end
end
%%
function sgn=find_changes(tests,type)
switch type
    case 'discrete'
        sgn=find(diff(tests));
    case {'sign','sgn','continuous','cont'}
        sgn=find(diff(sign(tests)));
end
end
function [res,pnew]=iterate_mon(mon,p,pref,bracket,type)
[val,pnew]=mon(p,pref);
switch type
    case 'discrete'
        res=val-mean(bracket)+0.1;
    case {'sign','sgn','continuous','cont'}
        res=val;
end
end
