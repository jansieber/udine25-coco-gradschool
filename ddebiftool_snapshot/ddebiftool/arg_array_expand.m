function [out,fmtrep]=arg_array_expand(dims,varargin)
%% expand args assuming format fmt for directional derivatives
%
% args is narg x (order+1) dim cell, first column are base args, others
% deviations each arg{k} is size(arg{k},1) x..x size(arg{k},dims(k)) x
% [vectorizations] or special symbols for deviations expand for
% vectorization as necessary.
%
% deviation arguments can be 1x2 cells {k,s} with integer k as first entry
% and special symbols s. The integer k determines before which vectorized
% dimension the symbol expansion gets placed. The symbol determines the
% type of expansion. Symbol s='I' replaces by num x num
% 
% if s is tself cell then it is interpreted as index vector, e.g. {2:3} is
% replaced by a collection of unit deviations in the indicated directions
% (so, eg, [e2,e3])
%
% 0's in deviations get expanded to full zero deviations
%
%% Example
%
% say format dims=[2,2] and format of 
% varargin={...
% 2x3,   {1,'I'}, 0,         2x3;
% 1x4x2,  0,    {2,{1,2:3}}, 1x4x1x3}
%
% base arguments have format 2x3 and 1x4 -> deviations must have same
% format. Two vectorization expansion are required directly by arguments:
% one (x2) in the provisionally 1st vectorization dimension for arg{2,1},
% and one (x3) in the provisionally 2nd vectorized dimension for arg{2,4}.
% Call the initial vectorized dimensions vecdim= 2x3, such that argument j
% has initially dimension fmt{j} x vecdim=fmt{j} x 2 x 3. There are two
% cell deviation arguments. Their first entries are 1 and 2, such that
% vecdim will be m1 x m2 x 2 x 3 where m1 is the size determined by 'I'
% (equal to 2*3=6) and m2 the dimension determined by {1,2:3} (equal to 2
% in this case).
%
% All arguments are initially expanded to shape fmt x 2 x 3, so 2x3x2x3 for
% arg{1,:}, 1x4x2x3 for args{2,:}. Cell arguments (with 'I' or {..} in
% second place) are initially NaN's of this size, Arguments with single 0
% are expanded to zeros of size fmt{j} x vecdim. Others are expanded by
% repmat.
%
% The argument with 'I' is then expanded to (here 2*3 x 2*3) identity,
% repeated with repmat to accommodate vectorization and placed at first
% dimension of vecdim (determined by first entry in cell). In the same way
% the cell argument with {2,{1,2:3}} is expanded into deviation
% reshape([e2,e3],1,4,2), repeated with repmat to accommodate all
% vectorizations (including those created prior by 'I'). Finally the
% dimensions created by the 'I' and {...} expansion get reshaped into the
% format given by the format or the {..}. So, for the above the outputs
% will be of shapes 
% 
%
% 2x3 x 2x3 x 1x2 x 2x3 (fmt{1}, 'I' with fmt{1}, dev {1,2:3},
% vectorizations)
% 
% and 1x4 x 2x3 x 1x2 x 2x3 (fmt{2}, 'I' with fmt{1}, dev {1,2:3},
% vectorizations)
%
% if the s in a cell argument for a deviation is a cell with more than one
% element, the dimensions will equal the shape of the cell. cell arguments
% and 'I' arguments should be cobmbined with 0's in other argument.
n_args=length(dims);
args=reshape(varargin,n_args,[]);
fmt=cellfun(@(a,f){size(a,1:f)},args(:,1),num2cell(dims(:)));
fmtrep=repmat(fmt(:),1,size(args,2));
%% shortcut if args are already ok
if check_shortcut(dims,args)
    out=args;
    return
end
is0=cellfun(@(a)isnumeric(a)&&numel(a)==1&&a==0,args);
iscf=cellfun(@(a)iscell(a)&&iscell(a{2}),args);
isid=cellfun(@(a)iscell(a)&&ischar(a{2})&&strcmp(a{2},'I'),args);
other=~isid&~iscf&~is0;
if any(other(:))
    sz=cellfun(@(a){[size(a),1]},args(other));
    args(other)=cellfun(@(a,f,s){reshape(a,[prod(f),s(numel(f)+1:end)])},...
        args(other),fmtrep(other),sz);
end
subi=@(fmt,ind)subsref(reshape(1:prod(fmt),[fmt,1]),struct('type','()','subs',{ind}));
cargs_sz=cell(size(args));
cargs_sz(isid)=fmtrep(isid);
if any(iscf(:))
    cargs_fmt=cellfun(@(f,a){subi(f,a{2})},fmtrep(iscf),args(iscf));
    cargs_sz(iscf)=cellfun(@(f){cellfun(@numel,f{2})},args(iscf));
    args(iscf)=cellfun(@(p,a){{p{1},reshape(a,[],1)}},args(iscf),cargs_fmt);
end
[out,pos_isc,ind_isc]=arg_vector_expand(size(args,1),args{:});
vecdims=size(out{end,end});
vecdimc=num2cell(vecdims(2:end));
vecdimc(pos_isc)=cargs_sz(ind_isc);
vecdims=cat(2,vecdimc{:});
out=cellfun(@(a,f){reshape(a,[f,vecdims,1])},out,fmtrep);
end
%%
function shortcut=check_shortcut(dims,args)
shortcut=~any(cellfun(@iscell,args(:)));
if ~shortcut
    return
end
dims=repmat(dims(:),1,size(args,2));
numa=cellfun(@numel,args);
diffnum=diff(numa,[],2);
shortcut=isempty(diffnum)||~any(diffnum(:));
if ~shortcut
    return
end
vecdims=arrayfun(@(i){[size(args{i},dims(i)+1:ndims(args{i})),1]},1:numel(args));
lendimdiff=diff(cellfun(@numel,vecdims));
shortcut=isempty(lendimdiff)|| ~any(lendimdiff);
if ~shortcut
    return
end
dimdiff=diff(cat(1,vecdims{:}),[],1);
shortcut=isempty(dimdiff) || all(dimdiff(:)==0);
end