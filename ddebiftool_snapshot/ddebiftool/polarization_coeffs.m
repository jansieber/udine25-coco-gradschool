%% Determine coefficients for polarization formula to compute higher-order derivatives
% The combinations in |mats| and |factors| help to reduce an arbitrary-order
% mixed derivative to a linear combination of single-directional
% derivatives.
%%
function [mats,factors,used]=polarization_coeffs(order,groups)
%% input
%
% * |order|: coefficients for mixed derivatives for degree |n=order| are
% computed
% * |groups|: cell array, forming a partition of 1:order (eg, {1:2,3} to
% indicate that linear deviation arguments 1 and 2 are equal), or keyword
% 'complex' indicating that |mats| and |factors| are for splitting up
% linear functional evaluation of complex argument into complex combination
% of real arguments (S=factors, P=mats): 
% flin[a+ib]^k=sum S_j flin(P(j,1)a+P(j,2)b)^k
%
%% Output
%
% * |mats|: nd x nterms matrix P of 1s and 0s or -1s
% * |factors|: nd x 1 array S
% * |used|: logical n x 1, if nterms<n: which of the factors should be used
% (if groups is non-trivial, eg, groups={1:2,3}, then
% used=[true;false;true];
%
% if groups is missing then nd=2^(n-1), otherwise, it can be shorter.
%
%% Usage of output
% (S=factors, P=mats) to obtain multilinear form B(x1,x2,...,xn), sum up
% S(j)*B(yj,yj,...,yj) over j where j=1..length(S),  and
% yj=P(j,1)*x1+...+P(j,1)*xn.
mats=(dec2bin(2^(order-1):2^order-1)-'0')*2-1;
factors=prod(mats,2)/factorial(order)/2^(order-1);
used=true(order,1);
if nargin==1||isempty(groups)||...
        (iscell(groups)&&length(groups)==order)
    return
end
%% sum coefficients for equal terms
if ischar(groups) && strcmp(groups,'complex')
    tr=tril(ones(order,order+1));
    factors=reshape(repmat(factors,1,order+1)*diag(arrayfun(@(k)(1i)^k*nchoosek(order,k),0:order)),[],1);
    mats=reshape(cat(3,mats*tr,mats*(1-tr)),[],2);
    used=true(2,1);
else
    mats=cell2mat(arrayfun(@(i){sum(mats(:,groups{i}),2)},1:length(groups)));
    sel=cellfun(@(g)g(1),groups);
    used=~used;
    used(sel)=true;
end
%% remove terms that result in all-zero coefficients
isnon0=any(mats~=0,2);
[mats,factors]=deal(mats(isnon0,:),factors(isnon0));
%% factor out common divisors
if size(mats,2)==1
    gc=mats(:,1);
else
    gc=gcd(mats(:,1),mats(:,2));
    for i=3:size(mats,2)
        gc=gcd(gc,mats(:,i));
    end
end
mats=mats./gc(:,ones(1,size(mats,2)));
factors=factors.*gc.^order;
%% unify signs: make sign first non-zero entry in mats positive
% and adjust factors
[dum,ix]=max(double(mats~=0),[],2); %#ok<ASGLU>
sgn=sign(mats(sub2ind(size(mats),(1:size(mats,1))',ix)));
mats=mats.*sgn(:,ones(size(mats,2),1));
factors=factors.*sgn.^order;
%% collect unique columns
[mats,dum,ind]=unique(mats,'rows'); %#ok<ASGLU>
factors=accumarray(ind,factors);
%% finally some factors may turn out zero, so remove them
isnon0=factors~=0;
[mats,factors]=deal(mats(isnon0,:),factors(isnon0));
end
