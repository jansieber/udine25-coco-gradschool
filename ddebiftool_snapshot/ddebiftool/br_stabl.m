function [branch,nunst,dom,triv_defect]=br_stabl(funcs,branch,varargin)
%% compute stability information along branch
% function st_branch=br_stabl(funcs,branch,skip,recompute)
%
% *INPUT*
%  
%  * |funcs| problem function
%  * |branch| 
%  * |skip| (opt. default 0) number of points to skip between stability computations
%  * |recompute| (opt, false) if nonzero recompute stability info already
%  present
%
% Optional parameters can also be passed on as name-value pairs. Other
% possible name-value pairs
%
% * |'exclude_trivial'| (logical): exclude trivial eigenvalues (what they are
%                            depends on the type of points)
% * |'exclude'| (function): @(point,z)false(size(z)) which eigenvalues
% should be removed from consideration (eg high freuqencies for renewal
% equations), first argument is point structure, second is vector of
% eigenvalues, return is vector of logicals of same size as z
% * |'pointtype_list'| structure with instruction for stability criterion
% * |'locate_trivial'| (function): e.g. @(p)[1;-1] for removing 1 and -1
%                            for period doubling orbits (overwrites standard location)
%
% *OUTPUT*
%
%	* branch branch with stability information
%   * nunst number of unstable eigenvalues
%   * dom dominant eigenvalue
%   * triv_defect error of known eigenvalues
%
% (c) DDE-BIFTOOL v. 1.02, 02/11/2000
%
% $Id$
%
%% backwards compatibility for previous optional arguments
args=varargin;
if ~isempty(args) && isnumeric(args{1})
    skip=args{1};
    args=args(2:end);
else
    skip=0;
end
if ~isempty(args) && isnumeric(args{1})
    recompute=args{1};
    args=args(2:end);
else
    recompute=false;
end
%% process options
defaults={'exclude_trivial',true,'locate_trivial',[],'critical',false,...
    'recompute',recompute,'skip',skip,'nunst_sign',false,...
    'pointtype_list',@pointtype_list,'exclude',@(p,z)false(size(z))};
[options,pass_on]=dde_set_options(defaults,args,'pass_on');
%% return if empty
np=length(branch.point);
if np<1 
    nunst=[];
    dom=[];
    triv_defect=[];
    return
end
%% take care of case for psol bifurcations
% which are stored as simple psol
fullpoints=branch.point;
fullfuncs=funcs;
[kind,points,funcs]=dde_get_kind(fullpoints,fullfuncs,'solution_for_stability');
%% (re)compute stability where required
irg=[1:options.skip+1:np-1,np];
mth=branch.method.stability;
mth=dde_set_options(mth,pass_on,'pass_on');
for i=length(irg):-1:1
    if isempty(points(irg(i)).stability) || options.recompute
        points(irg(i)).stability=p_stabil(funcs,points(irg(i)),mth,pass_on{:});
    end
end
branch.method.stability=mth;
%% for psol bifurcations re-insert extended components
branch.point=arrayfun(@(x,y)setfield(x,'stability',y.stability),fullpoints,points); 
if nargout==1
    return
end
%% check number of unstable ev etc
nunst=NaN(np,1+double(options.nunst_sign));
dom=NaN(np,1);
triv_defect=NaN(np,1);
evfcn=options.pointtype_list(pass_on{:}); 
if ~isempty(options.locate_trivial)
    % overwrite location of trivial eigenvalues
    evfcn.triv.(kind)=options.locate_trivial;
end
for i=irg
    p=points(i);
    %% check if point is special point
    if isfield(p,'flag') && ~isempty(p.flag) && isfield(evfcn.triv,p.flag)
        type=p.flag;
    else
        type=kind;
    end
    %% exclude trivial eigenvalues and requested eigenvalues
    trivialev=evfcn.triv.(type)(p);
    ev=evfcn.getev.(type)(p);
    evdist=evfcn.dist.(type);
    to_keep=~options.exclude(p,ev);
    if options.exclude_trivial && ~isempty(trivialev) && ~isempty(ev)
        if length(trivialev)>=length(ev)
            ev=[];
            to_keep=[];
        else
            iex=dde_match_complex(trivialev(:),ev(:),evdist);
            triv_defect(i)=norm(reshape(ev(iex),[],1)-trivialev(:),'inf');
            ev(iex)=[];
            to_keep(iex)=[];
        end
    end
    ev=ev(to_keep);
    %% count unstable eigenvalues
    unstsel=evfcn.stab.(type)(ev)>=0;
    if ~options.nunst_sign
        nunstloc=sum(unstsel);
    else
        nunstloc=[sum(unstsel&real(ev)<0),sum(unstsel&real(ev)>0)];
    end
    if isempty(nunstloc)
        nunst(i)=0;
    else
        nunst(i,:)=nunstloc;
    end
    %% check which other (unstable) eigenvalue is dominant 
    % (closest to stability boundary) if options.critical is set, it
    % includes stable eigenvalues immediately
    if any(nunst(i,:)>0) && ~options.critical
        [dum,ind]=min(evfcn.stab.(type)(ev(unstsel))); %#ok<ASGLU>
        dom(i)=ev(ind);
    else
        [dum,ind]=min(abs(evfcn.stab.(type)(ev)));
        if ~isempty(dum)
            dom(i)=ev(ind);
        end
    end
    if ~isnan(dom(i)) && imag(dom(i))<0
        dom(i)=conj(dom(i));
    end
end
end
