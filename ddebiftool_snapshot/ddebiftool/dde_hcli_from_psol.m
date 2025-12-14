function [hcli,stst]=dde_hcli_from_psol(psol,varargin)
%% convert periodic orbit to connecting orbit
default={'connections',1,'eqtime_ind',[],'linear',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% double psol
ntst=length(psol.mesh)-1;
psol=dde_psol_create('point',psol,...
    'mesh',[psol.mesh/2,0.5+psol.mesh(2:end)/2],...
    'period',2*psol.period,...
    'profile',[psol.profile,psol.profile(:,2:end)]);
if isempty(options.eqtime_ind)
    %% find point(s) w min derivative
    npd=sqrt(sum(dde_coll_eva(psol,psol.mesh,'diff',1).^2,1));
    [dum, pos]=min(npd(1:ntst-1)); %#ok<*ASGLU>
    epos=pos+ntst;
    %% find other local minima, sorted according to value, if nconnections>1
    % ibounds are the indices at which the orbit is supposed to be cut
    imin=find(diff(sign(diff(npd(pos:epos))))>0)+pos;
    nmin=length(imin);
    assert(nmin>=options.connections-1,'dde_hcli_from_psol:connections',...
        ['dde_hcli_from_psol: requested %d connecting orbits but found only',...
        '%d local minima of derivative'],options.connections,nmin+1);
    [dum,is]=sort(npd(imin));
    imin=imin(is);
    imin=imin(1:options.connections-1);
    imin=sort(imin);
    ibounds=[pos,imin; imin,epos];
else
    %% user has provided hints
    ibounds=[options.eqtime_ind;...
        options.eqtime_ind(2:end),options.eqtime_ind(1)+ntst];
end
for i=size(ibounds,2):-1:1
    %% cut orbit at requested points
    coll=dde_psol_cut(psol,ibounds(:,i));
    hcli(i)=dde_hcli_create('point',coll,...
        'x1',coll.profile(:,1),'x2',coll.profile(:,end));
    %% correct equilibira and compute spectrum if requested
    for k=2:-1:1
        stst(k,i)=dde_stst_create('x',coll.profile(:,1),  'parameter',hcli(i).parameter);
        stst(k,i)=stst_stability(stst(k,i),pass_on{:});
        stst(k,i)=stst_correct(stst(k,i),sprintf('(%d,%d)',k,i),pass_on{:});
        stst(k,i)=stst_stability(stst(k,i),pass_on{:});
    end 
    %% add linear equilibrium information
    if options.linear 
        [hcli(i),stst(:,i)]=dde_hcli_from_hcli(hcli(i),pass_on{:});
    end
end
end
%%
function coll=dde_psol_cut(psol,bd)
[psol,submesh]=dde_coll_check(psol);
icoarse=1:psol.degree:length(psol.mesh);
[icoarse,~,ibd]=unique([bd(:)',icoarse]);
tcoarse=psol.mesh(icoarse(ibd(1):ibd(2)));
t=dde_coll_meshfill(tcoarse,1,'grid',submesh);
prof=dde_coll_eva(psol,t);
tscal=t-t(1);
T=tscal(end);
tscal=tscal/T;
coll=dde_coll_create('mesh',tscal,'parameter',psol.parameter,...
    'period',T*psol.period,'profile',prof,'degree',psol.degree);
end
%%
function stst=stst_correct(stst,label,varargin)
default={'free_par',[],'funcs',[],'stst_correct',true,...
    'method',getfield(df_mthod('stst'),'point'),'print',1};
[opts,pass_on]=dde_set_options(default,varargin,'pass_on');
if ~opts.stst_correct
    return
end
assert(~isempty(opts.funcs),'dde_hcli_from_psol:funcs',...
    'dde_hcli_from_psol: argument ''funcs'' is needed for stst correction');
[stst,suc]=p_correc(opts.funcs,stst,opts.free_par,[],opts.method,0,stst,pass_on{:});
if opts.print
    fprintf('stst %s correction: %d\n',label,suc);
end
assert(suc,'dde_hcli_from_psol:correct',...
    'dde_hcli_from_psol: correction failed for stst %s');
end
%%
function stst=stst_stability(stst,varargin)
default={'funcs',[],'stst_stability',true,...
    'method',getfield(df_mthod('stst'),'stability')};
[opts,pass_on]=dde_set_options(default,varargin,'pass_on');
if ~opts.stst_stability
    return
end
assert(~isempty(opts.funcs),'dde_hcli_from_psol:funcs',...
    'dde_hcli_from_psol: argument ''funcs'' is needed for stst stability');
stst.stability=p_stabil(opts.funcs,stst,opts.method,pass_on{:});
end
