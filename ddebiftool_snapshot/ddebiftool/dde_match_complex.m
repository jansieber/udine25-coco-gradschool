function ind=dde_match_complex(z1,z2,costfcn)
if nargin<3
    costfcn=@(x,y)abs(x-y);
end
cost=costfcn(z1(:,ones(length(z2),1)),z2(:,ones(length(z1),1)).');
ind=munkres(cost);
end
