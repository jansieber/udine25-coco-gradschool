function y=fun_vec(f,dims,isvec,nf,varargin)
%% wrapper to pseudo-vectorize f if needed
if isnumeric(dims)
    argdims=dims(:);
    outdims=1;
else
    argdims=dims{2}(:);
    outdims=dims{1};
end
orig_args=arg_array_expand(argdims,varargin{:});
sz=cellfun(@(x){size(x)},orig_args(:));
vecdim=[sz{1}(argdims(1)+1:end),1];
nvec=prod(vecdim);
%% empty arguments and output dimension provided
if nvec==0 
    y=zeros([nf,vecdim]);
    return
end
szbase=cellfun(@(s,d){s(1:d)},sz,num2cell(argdims));
lenbase=cellfun(@prod,szbase);
lenbase_c=num2cell(lenbase);
args=cellfun(@(x,len){reshape(x,[len,nvec])},orig_args(:),lenbase_c);
if any(isvec)
    len_nonvecblocks=diff([0,find(any(diff(cat(1,args{~isvec}),[],2)~=0,1)),nvec]);
else
    len_nonvecblocks=ones(1,nvec);
end
n_nonvecblocks=length(len_nonvecblocks);
if n_nonvecblocks==1 % fully vectorized
    y=f(orig_args{:});
else % not or partially vectorized
    argblocks=mat2cell(cat(1,args{:}),lenbase,len_nonvecblocks);
    argblocks(~isvec,:)=cellfun(@(a){a(:,1)},argblocks(~isvec,:));
    szrep=repmat(szbase,1,n_nonvecblocks);
    argblocks(~isvec,:)=cellfun(@(a,s){reshape(a,[s,1])},   argblocks(~isvec,:),szrep(~isvec,:));
    argblocks( isvec,:)=cellfun(@(a,s,len){reshape(a,[s,numel(a)/len])},...
        argblocks( isvec,:),szrep( isvec,:),lenbase_c(isvec,ones(1,n_nonvecblocks)));
    yc=arrayfun(@(i){f(argblocks{:,i})},1:n_nonvecblocks);
    dimcat=max(cellfun(@ndims,yc));
    y=cat(dimcat,yc{:});
end
nf=size(y,1:outdims);
y=reshape(y,[nf,vecdim]);
end

