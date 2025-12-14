function xpattern=dde_join_xpattern(newpattern,xpattern)
if nargin<2
    xpattern=zeros(2,0);
end
if size(xpattern,1)==3 && size(newpattern,1)==2
    newpattern=cat(1,newpattern,ones(1,size(newpattern,2)));
end
if size(newpattern,1)==3 && size(xpattern,1)==2
    xpattern=cat(1,xpattern,ones(1,size(xpattern,2)));
end
xpattern=cat(2,xpattern,newpattern);
xpattern=unique(xpattern(end:-1:1,:)','rows')';
xpattern=xpattern(end:-1:1,:);
end
