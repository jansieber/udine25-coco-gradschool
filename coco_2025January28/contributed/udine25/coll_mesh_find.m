function [ix,x]=coll_mesh_find(tbp,x,varargin)
default={'assert',true};
options=sco_set_options(default,varargin,'pass_on');
x=x(:)';
if options.assert
    assert(all(x>=tbp(1))&&all(x<=tbp(end)));
end
ntst=size(tbp,2);
locdefault=ones(1,numel(x));
locdefault([x(1:end-1)==x(2:end),false])=-1;
%% options
default={'loc',locdefault};
options=sco_set_options(default,varargin,'pass_on');
%% fill in complete mesh
%% locate in which subinterval each x is (ix)
mcoarse=[tbp(1,:),tbp(end,end)];
ix=floor(interp1(mcoarse,1:ntst+1,x,'linear'));
%% adjust to lower interval if requested or required at final boundary
l_adj= (x==mcoarse(ix) & options.loc(:)'==-1 & ix>1) | ix==ntst+1;
ix(l_adj)=ix(l_adj)-1;
end