function c=find_chain(nbh)
entries=unique(nbh(:));
loc=zeros(max(entries),1);
for i=1:size(nbh,1)
    if loc(nbh(i,1))==1 || loc(nbh(i,2))==2
        nbh(i,:)=nbh(i,[2,1]);
    end
    loc(nbh(i,1:2))=1:2;
end
nbh=sortrows(nbh);
c=ones(1,length(entries))*nbh(1,1);
for i=2:length(entries)
    c(i)=nbh(c(i-1),2);
end
end