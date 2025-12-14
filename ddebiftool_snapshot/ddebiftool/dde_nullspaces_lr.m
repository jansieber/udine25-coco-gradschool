function [v,w]=dde_nullspaces_lr(A,varargin)
%% determine left and right null space of large sparse matrix (possibly non-square)
% This routine relies on the property of the LU factorization that it puts
% the zero elements at the bottom of U
%
% $Id: dde_nullspaces_lr.m 319 2019-01-31 02:14:56Z jansieber $
%% 
default={'nulltol',1e-8,'bordered',true,'check',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% obtain first (unsafe?) approximations using LU decomposition
% for left nullspace
v=dde_nullspace(A,'check',false,'nulltol',options.nulltol,pass_on{:});
[s1,s2]=size(A);
rightnulldim=size(v,2);
%% for right nullspace
leftnulldim=s1-s2+rightnulldim;
assert(leftnulldim>=0);
if options.bordered && ~dde_isoctave()
    wrn=warning('off','dde_nullspace:residual');
end
w=dde_nullspace(A',pass_on{:},'nulldim',leftnulldim);
assert(leftnulldim==size(w,2));
if options.bordered && ~dde_isoctave()
    warning(setfield(wrn,'state','on')); %#ok<SFLD>
end
if ~options.bordered
    return
end
%% obtain safer(?) nullspaces using bordered matrix (unnessecary?)
[v,w]=null_update(A,v,w,rightnulldim,leftnulldim);
[errvlarge,errwlarge]=check_null(A,v,w,options);
if ~isempty(errvlarge) || ~isempty(errwlarge)
    % repeat attempt
    [v,w]=null_update(A,v,w,rightnulldim,leftnulldim);
    [errvlarge,errwlarge,resv,resw]=check_null(A,v,w,options);    
    if ~isempty(errvlarge)
        [dum,ix]=max(abs(errvlarge)); %#ok<ASGLU>
        warning('dde_nullspaces_lr:residual',...
            'dde_nullspaces_lr: residual for right nullvector %d =%g\n',...
            errvlarge(ix),resv(errvlarge(ix)));
    end
    if ~isempty(errwlarge)
        [dum,ix]=max(abs(errwlarge)); %#ok<ASGLU>
        warning('dde_nullspaces_lr:residual',...
            'dde_nullspaces_lr: residual for left nullvector %d =%g\n',...
            errwlarge(ix),resw(errwlarge(ix)));
    end    
end
end
%%
function [v,w]=null_update(A,v,w,rightnulldim,leftnulldim)
Aext=[A,w;v',sparse(rightnulldim,leftnulldim)];
[s1,s2]=size(A);
assert(size(Aext,1)==size(Aext,2))
lrhs=cat(1,zeros(s1,rightnulldim),eye(rightnulldim));
v=Aext\lrhs;
[v,~]=qr(v(1:s2,:),0);
rrhs=cat(1,zeros(s2,leftnulldim),eye(leftnulldim));
w=Aext'\rrhs;
[w,~]=qr(w(1:s1,:),0);
end
%%
function [errvlarge,errwlarge,resv,resw]=check_null(A,v,w,options)
if options.check
    resv=max(abs(A*v),[],1);
    resw=max(abs(A'*w),[],1);
    An=norm(A,'inf');
    errvlarge=find(isnan(resv(:))|resv(:)>options.nulltol*An);
    errwlarge=find(isnan(resw(:))|resw(:)>options.nulltol*An);
end
end