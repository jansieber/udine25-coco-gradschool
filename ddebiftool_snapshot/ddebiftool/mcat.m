function P=mcat(dim,varargin)
%% concatenate sparse matrices in hogher dimensions
sz=cellfun(@(A){mshape(A,'size')},varargin);
sz=cell2mat(reshape(sz,[],1));
ini=[0;cumsum(sz(:,dim))];
szdim=ini(end);
ini=ini(1:end-1);
assert(norm(diff(sz(:,[1:dim-1,dim+1:end]),[],1))==0);
irc=cellfun(@(A){mshape(A,'irc')},varargin);
vals=cellfun(@(A){mshape(A,'vals')},varargin);
for i=1:size(sz,1)
    irc{i}(:,dim)=irc{i}(:,dim)+ini(i);
end
P=mshape(zeros(0,1),'struct','skip');
ircP=num2cell(cat(1,irc{:}),1);
P.sz=sz(1,:);
P.sz(dim)=szdim;
P.ir=sub2ind(P.sz,ircP{:});
P.vals=cat(1,vals{:});
P=mshape(P,'state','sparse','mat');
end