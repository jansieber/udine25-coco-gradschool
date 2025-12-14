function [v,w,evs]=dde_svdspaces_lr(A,nulldim,varargin)
%% determine left and right null space of large sparse matrix (possibly non-square)
% This routine relies on the property of the LU factorization that it puts
% the zero elements at the bottom of U
%
% $Id: dde_svdspaces_lr.m 309 2018-10-28 19:02:42Z jansieber $
%% 
default={'nulltol',1e-8,'nullspaces',true,'use_svds',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if options.nullspaces
    [v,w]=dde_nullspaces_lr(A,'nulltol',options.nulltol,pass_on{:});
    tolnulldim=size(v,2);
    if tolnulldim==nulldim
        evs=zeros(nulldim,1);
        return
    end
else
    v=zeros(size(A,2),0);
    w=zeros(size(A,1),0);
    tolnulldim=0;
end
evdim=nulldim-tolnulldim;
J1=[A,w;v',sparse(tolnulldim,size(w,2))];
if ~options.use_svds
    [s1,s2]=size(J1);
    Je=[sparse(s1,s1),J1; J1',sparse(s2,s2)];
    Je=Je+speye(size(Je))*1e-15*1i;
    [L,U,P,Q,R]=lu(Je);
    %[Lt,Ut,Pt,Qt,Rt]=lu(J1');
    opts.issym=0;
    opts.isreal=0;
    [V,D]=eigs(@(v)Q*(U\(L\(P*(R\v)))),...
        s1+s2,2*evdim,'sm',opts);
    %[W,~]=eigs(@(v)(Qt*(Ut\(Lt\(Pt*(Rt\(Q*(U\(L\(P*(R\v)))))))))),...
    %    size(J1,1),evdim,-1e-6,opts);
    sv=diag(D)-1e-15*1i;
    [~,is]=sort(abs(sv));
    sv=real(sv(is));
    we=real(V(1:s1,is));
    ve=real(V(s1+1:end,is));
    sel=sv>=0;
    evs=sv(sel);
    ve=ve(1:size(A,2),sel);
    we=we(1:size(A,1),sel);
else
    if dde_isoctave()
        [we,D,ve]=svds(J1,evdim,options.nulltol);
    else
        [we,D,ve]=svds(J1,evdim,'smallest');
    end
    evs=diag(D);
    [evs,ix]=sort(evs);
    we=we(:,ix);
    ve=ve(:,ix);
end
[v,~]=qr([v,ve],0);%qr([v,V(1:size(A,2),:)],0);
[w,~]=qr([w,we],0);%qr([w,W(1:size(A,1),:)],0);
end
