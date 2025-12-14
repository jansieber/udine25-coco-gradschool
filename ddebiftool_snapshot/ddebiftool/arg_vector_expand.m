function [out,pos_isc,ind_isc]=arg_vector_expand(n,varargin)
%% expand n column vector arguments for directional derivatives
%
% varargin are n*(order+1) arrays or cells, and will be reshaped to
% args=reshape(n,varargin). First column of args are base args, others
% deviations. Each args(j,k) is dim(j) x [vectorizations] or expansion
% instructions for deviations expand for vectorization as necessary, and
% dim(j)=size(args{j,1},1).
%
% deviation arguments can be 1x2 cells {k,rg} with integer k as first entry
% and range rg. The integer k determines at which vectorized dimension the
% range gets expanded. The symbol 'I' is a shortcut for range 1: dim(j) 
% 
% The range is interpreted as index vector, e.g. 2:3 is replaced by a
% collection of unit deviations of dimension dim(j) in the indicated
% directions (so, eg, [e2,e3], where e2, e3 may have length 4)
%
% 0's in deviations get expanded to full zero deviations
%
%% Example
% n=2, length(varargin)=8, so
% args={...
% 6x1   {1,'I'}    0         6x1
% 4x2       0    {2,2:3}     4x1x3
%
% base arguments have length (first dim) 6 and 4 -> deviations must have
% same first dim format. Two vectorization expansion are required directly
% by arguments: one (x2) in the provisionally 1st vectorization dimension
% for arg{2,1}, and one (x3) in the provisionally 2nd vectorized dimension
% for arg{2,4}. There are two expansion instruction arguments. Their first
% entries are 1 and 2, such that vecdim will be m1 x m2 x 2 x 3 where m1
% and m2 are the lengths of the ranges. The symbol 'I' is replaced by
% 1:dim(1) such that m1=6 and m2=2 (determined by 2:3 in this case.
%
% All arguments are initially expanded to shape dim x 2 x 3, so 6x2x3 for
% arg{1,:}, 4x2x3 for args{2,:}. Cell arguments (with 'I' or a range in
% second place) are initially NaN's of this size, Arguments with single 0
% are expanded to zeros of size dim x vecdim. Others are expanded by
% repmat.
%
% The argument with 'I' is then expanded to (here 6 x 6) identity, repeated
% with repmat to accommodate vectorization and placed at first dimension of
% vecdim (determined by first entry in cell). In the same way the cell
% argument with {2,2:3} is expanded into deviation [e2,e3] of 4x1 unit
% vector columns, repeated with repmat to accommodate all vectorizations
% (including those created prior by 'I'). Finally the dimensions created by
% the 'I' and range expansion get reshaped into the format given by the
% format or the {..}. So, for the above the outputs will be of shapes
%
% 6 x 6 x 2 x 2 x 3 (dim(1), 'I' with dim(1), range 2:3, vectorizations)
% 
% and 
% 
% 4 x 6 x 2 x 2 x 3 (dim(2), 'I' with dim(1), range 2:3, vectorizations)
%
% Cell arguments and 'I' arguments should be combined with 0's in other
% argument.
args=reshape(varargin,n,[]);
dim=repmat(cellfun(@(a){size(a,1)},args(:,1)),1,size(args,2));
%% replace 'I' with appropriate range argument
isid=cellfun(@(a)iscell(a)&&ischar(a{2})&&strcmp(a{2},'I'),args);
if any(isid(:))
    [ir,ic]=find(isid); %#ok<ASGLU>
    args(isid)=cellfun(@(a,d){[a(1),{1:d}]},args(isid),dim(ir,1));
end
is0=cellfun(@(a)isnumeric(a)&&numel(a)==1&&a==0,args);
isc=cellfun(@(a)iscell(a)&&isnumeric(a{2}),args);
pos_isc=reshape(cellfun(@(a)a{1},args(isc)),1,[]); % which args are cells
fmt_isc=reshape(cellfun(@(a)a(2),args(isc)),1,[]); % save format
if isempty(fmt_isc) && isnumeric(fmt_isc)
    fmt_isc=num2cell(fmt_isc);
end
len_isc=cellfun(@length,fmt_isc); % number of columns to be expanded
args(isc)=cellfun(@(f){NaN(f,1)},  dim(isc)); % fill in cell args with NaN of right shape
args(is0)=cellfun(@(f){zeros(f,1)},dim(is0)); % expand singleton zeros as needed
maxdim=max(cellfun(@ndims,args(:))); % maximal dimension of all args
argdims=cell2mat(cellfun(@(a){[size(a),ones(1,maxdim-ndims(a))]},args(:))); % fill up all dimensinos with 1's to turn into array
n_args=size(argdims,1);
nvec=max(argdims,[],1); % vectorizations
ind_isc=find(isc);
n_isc=length(ind_isc);
curdims=[len_isc,nvec(2:end)];
n_vecdims=length(curdims);
argdims=[argdims(:,1),ones(n_args,n_isc),argdims(:,2:end)];
%% take into account possibility of zero-length arguments
is0dim=any(argdims(:,2:end)==0,1);
curdims(is0dim)=0;
isnon0dim=true(1,n_vecdims);
isnon0dim(curdims==0)=false;
reps=ones(n_args,n_vecdims);
reps(:,curdims==0)=0;
reps(:,isnon0dim)=curdims(ones(n_args,1),isnon0dim)./argdims(:,[false,isnon0dim]);
reps=num2cell([ones(n_args,1),reps],2);
argdimc=num2cell(argdims,2);
args(~isc)=cellfun(@(a,d,r){repmat(reshape(a,d),r)},...
    reshape(args(~isc),size(argdimc(~isc))),argdimc(~isc),reps(~isc));
%% expand all ranges
% expanding   all other arguments simultaneously, initially put expansions
% at the beginning
args(ind_isc)=cellfun(@(pos,rg,a){rg_create(pos,rg,size(a,1),curdims)},...
    num2cell((1:n_isc)'),fmt_isc(:),reshape(args(ind_isc),n_isc,1));
%% reorder to put expansions in the right place
dperm=NaN(1,n_vecdims);
dperm(pos_isc)=1:n_isc;
dperm(isnan(dperm))=n_isc+1:n_vecdims;
out=cellfun(@(a){permute(reshape(a,[size(a,1),curdims]),[1,dperm+1])},args);
end
%%
function ar=rg_create(pos,rg,sz,dims)
reps=length(rg);
mat=full(sparse(rg(:),(1:reps)',ones(reps,1),sz,reps));
fmt=ones(size(dims));
fmt(pos)=reps;
rep=dims;
rep(pos)=1;
mat=reshape(mat,[sz,fmt]);
ar=repmat(mat,[1,rep]);
end