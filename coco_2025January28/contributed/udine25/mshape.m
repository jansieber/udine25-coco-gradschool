function Ar=mshape(A,varargin)
%% reshape, permute, repmat, slice potentially sparse matrix A 
% with intermittent higher-dimensional object
altc={'r','r';   'rs','r';       'reshape','r';...
    'p','p';     'permute','p';...
    'rep','rep'; 'repmat','rep';...
    'sl','sl';   'slice','sl';    'assign','sl'; 'asgn','sl';...
    'sm','sm';   'submatrix','sm';'smasgn','sm';...
    'mat','mat';  'matrix','mat';...
    'str','str';  'struct','str';...
    'vals','vals'; 'values','vals';...
    'irc','irc';   'apply','apply';  'skip','skip';  'state','state';...
    'blkdiag','blkdiag'; 'diagblk','diagblk';...
    'size','size'; 'ir','ir'; 'nrmlz','normalize'; 'normalize','normalize'};
altc=altc';
fcn=struct(altc{:});
[Ar,trafo]=sparse2struct(A);
state=getstate(Ar);
cmds.sparse=struct('r',{{@sparse_reshape,1}},'p',{{@sparse_permute,1}},...
    'rep',{{@sparse_repmat,1}},'size',{{@(A)A.sz,0}},...
    'mat',{{@(M)struct2sparse(M,1),0}},'str',{{@(M)sparse2struct(M,1),0}},...
    'irc',{{@A2irc,0}},'sm',{{@subsparse,1}},...
    'vals',{{@(M)getfield(sparse2struct(M),'vals'),0}},...
    'ir',{{@(M)getfield(sparse2struct(M),'ir'),0}},...
    'normalize',{{@sparse_normalize,0}});
cmds.full=struct('r',{{@(M,a)reshape(M,a{:}),1}},'p',{{@mat_permute,1}},...
    'rep',{{@repmat,1}},'size',{{@size,0}},'mat',{{@struct2full,0}},...
    'str',{{@mat2struct,0}},'irc',{{@(M)A2irc(mat2struct(M)),0}},...
    'sm',{{@mat_subfull,1}},...
    'vals',{{@(M)getfield(mat2struct(M),'vals'),0}},...
    'ir',{{@(M)getfield(mat2struct(M),'ir'),0}},...
    'normalize',{{@(M)M,0}});
cmds.full.sl={@(M,a)mat_slice(M,a,cmds.full),1};
cmds.sparse.sl={@(M,a)mat_slice(M,a,cmds.sparse),1};
[cmds.full.apply,cmds.sparse.apply]=deal({@(M,fun)fun(M),1});
[cmds.full.skip,cmds.sparse.skip]=deal({@(M)M,0});
[cmds.full.blkdiag,cmds.sparse.blkdiag]=deal({@mat_blkdiag,0});
[cmds.full.diagblk,cmds.sparse.diagblk]=deal({@mat_diagblk,1});
i=0;
while i<length(varargin)
    i=i+1;
    cmd=varargin{i};
    if strcmp(cmd,'state')
        state=varargin{i+1};
        i=i+1;
        continue
    elseif i==length(varargin) && strcmp(cmd,'skip')
        return
    end
    fc=cmds.(state).(fcn.(cmd));
    f=fc{1};
    nargs=fc{2};
    argc=varargin(i+(1:nargs));
    Ar=f(Ar,argc{:});
    i=i+nargs;
end
Ar=struct2sparse(Ar,trafo);
end
function Ar=nsparse_create(varargin)
default={'sz',[1,1],'vals',zeros(0,1),'ir',ones(0,1)};
Ar=struct(default{:});
user=reshape(varargin,2,[]);
for i=1:size(user,2)
    Ar.(user{1,i})=user{2,i};
end
if isfield(Ar,'size')
    Ar.sz=Ar.size;
    Ar=rmfield(Ar,'size');
end
end
%%
function Ar=sparse_reshape(Ar,arg)
empty=find(cellfun(@isempty,arg));
if ~isempty(empty)
    arg{empty}=prod(Ar.sz)/prod(cat(2,arg{:}));
end
sz2=cat(2,arg{:});
assert(prod(sz2)==prod(Ar.sz)&&all(sz2==round(sz2)));
Ar.sz=sz2;
end
%%
function Ar=sparse_permute(Ar,arg)
irc=A2irc(Ar);
irc=irc(:,arg);
Ar.sz=Ar.sz(arg);
Ar=irc2A(Ar,irc);
end
%%
function Anew=mat_permute(Am,arg)
len=length(arg);
if len==length(unique(arg))
    Anew=permute(Am,arg);
else
    sz=[size(Am),ones(1,len-ndims(Am))];
    irc=mat_s2s(numel(Am),sz,(1:numel(Am))');
    irc=irc(:,arg);
    sz=sz(arg);
    ir=mat_s2s(sz,numel(Am),irc);
    Anew=zeros(sz);
    Anew(ir)=Am(:);
end
end
%%
function Ar=sparse_repmat(Ar,arg)
ndadd=length(arg)-length(Ar.sz);
sz1=[Ar.sz,ones(1,ndadd)];
sz2=sz1.*arg;
nmult=prod(arg);
Ar.vals=repmat(Ar.vals,[nmult,1]);
irc=A2irc(Ar);
irc=[irc,ones(size(irc,1),ndadd)];
for i=1:length(sz2)
    nirc_orig=size(irc,1);
    irc=repmat(irc,arg(i),1);
    nirc=size(irc,1);
    sziadd=reshape(repmat((1:arg(i)-1)*sz1(i),nirc_orig,1),nirc-nirc_orig,1);
    irc(nirc_orig+1:end,i)=irc(nirc_orig+1:end,i)+sziadd;
end
Ar.sz=sz2;
Ar=irc2A(Ar,irc);
end
%%
function Ar=mat_slice(Ar,arg,cmds)
sz=cmds.size{1}(Ar);
[asgn,rg,len,rhsval]=check_asgn(sz,arg); %#ok<ASGLU>
lch=cellfun(@ischar,arg);
lsl=~lch;
ich=find(lch);
isl=find(lsl);
invperm([ich,isl])=1:length(sz);
sub=cellfun(@(i)i(:),arg(lsl),'UniformOutput',false);
ind=sub2ind([sz(isl),1],sub{:});
Ar=cmds.p{1}(Ar,[ich,isl]);
Ar=cmds.r{1}(Ar,{prod(sz(ich)),prod(sz(isl))});
Am=struct2sparse(Ar);
if ~asgn
    Am=Am(:,ind);
    sznew=sz([ich,isl]);
    sznew(length(ich)+2:end)=1;
    sznew(length(ich)+1)=length(ind);
else
    Am(:,ind)=reshape(rhsval,[],numel(ind));
    sznew=sz([ich,isl]);
end
Ar=sparse2struct(Am);
Ar=cmds.r{1}(Ar,num2cell(sznew));
Ar=cmds.p{1}(Ar,invperm);
end
%%
function Am=mat_subfull(Am,arg)
sz=size(Am);
[asgn,rg,len,rhsval,cleanrg]=check_asgn(sz,arg);
nsz=length(sz);
mshc=cell(1,nsz);
[mshc{:}]=ndgrid(rg{:});
msh=reshape(cat(nsz+1,mshc{:}),[],nsz);
[sz,Am]=enlarge_zeros(asgn,msh,sz,Am);
ir=mat_s2s(sz,prod(sz),msh);
if asgn
    assert(all(cleanrg));
    Am(ir)=rhsval;
else
    Am=reshape(Am(ir),len);
end
end
%%
function As=subsparse(As,arg)
sz=As.sz; 
[asgn,rg,len,rhsval,cleanrg,is_ch]=check_asgn(sz,arg);
% ~cleanrg are index ranges that contain duplicates, requiring allocation
i_rgmv=find(~is_ch); % which index ranges are changing?
if asgn
    assert(all(cleanrg));
    irca=A2irc(As);
    asel=true(size(irca,1),1);
    for i=i_rgmv
        asel=asel&ismember(irca(:,i),rg{i});
    end
    B=mshape(rhsval,'state','sparse','struct','r',num2cell(len));
    ircb=A2irc(B);
    [sz,As]=enlarge_zeros(asgn,cellfun(@numel,rg),sz,As);
    for i=i_rgmv
        ircb(:,i)=rg{i}(ircb(:,i));
    end
    irb=mat_s2s(sz,prod(sz),ircb);
    As.ir=[As.ir(~asel);irb];
    As.vals=[As.vals(~asel);B.vals];
else
    As=sparse_fill_rg(As,rg,len,find(~cleanrg));
    irca=A2irc(As);
    asel=true(numel(As.ir),1);
    for i=i_rgmv
        asel=asel&ismember(irca(:,i),rg{i});
    end
    As.vals=As.vals(asel);
    irca=irca(asel,:);
    for i=i_rgmv
        trans(rg{i})=1:len(i); %#ok<AGROW>
        irca(:,i)=trans(irca(:,i));
    end
    As.sz=len;
    As=irc2A(As,irca);
end
end
%%
function irc2=mat_s2s(sz1,sz2,irc1)
irc1c=num2cell(irc1,1);
lin=sub2ind([sz1,1],irc1c{:});
nout=length(sz2);
irc2c=cell(1,nout);
[irc2c{:}]=ind2sub(sz2,lin);
irc2=cat(2,irc2c{:});
end
%%
function irc=A2irc(A)
irc=mat_s2s([prod(A.sz),1],A.sz,A.ir);
end
%%
function A=irc2A(A,irc)
ircc=num2cell(irc,1);
A.ir=sub2ind(A.sz,ircc{:});
end
%%
function Am=struct2sparse(As,trafo)
if nargin<2
    trafo=true;
end
if trafo && isstruct(As) && (length(As.sz)<=2||all(As.sz(3:end)==1))
    irc=A2irc(As);
    Am=sparse(irc(:,1),irc(:,2),As.vals,As.sz(1),As.sz(2));
else
    Am=As;
end
end
%%
function Am=struct2full(As)
if isnumeric(As)
    Am=full(As);
else
    Am=zeros(As.sz);
    Am(As.ir)=As.vals;
end
end
%%
function [As,trafo]=sparse2struct(Am,force)
if nargin<2
    force=false;
end
if (force&&isnumeric(Am))||issparse(Am)
    [ir,ic,vals]=find(Am);
    As=nsparse_create('sz',size(Am),'vals',vals(:),...
        'ir',sub2ind(size(Am),ir(:),ic(:)));
    trafo=true;
elseif iscell(Am)
    As=nsparse_create(Am{:});
    trafo=false;
else
    As=Am;
    trafo=false;
end
end
%%
function As=mat2struct(Am)
As=nsparse_create('vals',Am(:),'sz',size(Am),'ir',(1:numel(Am))');
end
%%
function state=getstate(Ar)
if isnumeric(Ar) && ~issparse(Ar)
    state='full';
else
    state='sparse';
end
end
%% normalize sparse structure
function As=sparse_normalize(As)
if isstruct(As)
    [As.ir,ix]=sort(As.ir);
    As.vals=As.vals(ix);
    if length(As.sz)>2 && all(As.sz(3:end)==1)
        As.sz=As.sz(1:2);
    end
end
end
%% check if assignment or selection
function [asgn,rg,len,rhsval,isclean,is_ch]=check_asgn(sz,arg)
%% check if argument is provided for assignment
if length(sz)<length(arg) && length(arg)>2 && strcmp(arg{end-1},'=')
    [arg,rhsval]=deal(arg(1:end-2),arg{end});
    asgn=true;
else
    asgn=false;
    rhsval=[];
end
[rg,len,isclean,is_ch] = char2rg(arg,sz);
end
%% enlarge matrix if required by submatrix assignment or range duplication
function [sz,A]=enlarge_zeros(asgn,ircnew,sz,A)
if ~asgn || all(max(ircnew,[],1)<=sz) % enlargement not necessary
    return
end
sparse=isstruct(A);
newsz=max([ircnew;sz],[],1);
if sparse % sparse case
    ir=A.ir;
else
    ir=(1:prod(sz))';
end
ircold=mat_s2s(prod(sz),sz,ir);
irnew=mat_s2s(newsz,prod(newsz),ircold);
if sparse
    A=nsparse_create('sz',newsz,'ir',irnew,'vals',A.vals);
else
    Anew=zeros(newsz);
    Anew(irnew)=A(:);
    A=Anew;
end
sz=newsz;
end
%%
function [rg,len,isclean,l_ischar] = char2rg(arg,sz)
sz=[sz,ones(1,max(0,length(arg)-length(sz)))];
l_isrg=cellfun(@isnumeric,arg);
l_ischar=~l_isrg;
rg=cell(1,length(sz));
rg(l_isrg)=cellfun(@(i)i(:),arg(l_isrg),'UniformOutput',false);
rg(~l_isrg)=arrayfun(@(i)(1:i)',sz(~l_isrg),'UniformOutput',false);
len=cellfun(@length,rg);
isdirty=find(l_isrg);
clean_num=arrayfun(@(i)all(accumarray(rg{i},1)<=1),isdirty);
isdirty(clean_num)=[];
isclean=true(1,length(sz));
isclean(isdirty)=false;
end
%%
function As=sparse_fill_rg(As,rg,len,dirtyrg)
perm=(1:length(len));
for i=1:length(dirtyrg)
    perm=[setdiff(perm,dirtyrg(i)),dirtyrg(i)];
    invperm(perm)=1:length(len); %#ok<AGROW>
    sz=As.sz(perm);
    Am=mshape(As,'p',perm,'r',{[],sz(end)},'mat');
    Am=Am(:,rg{dirtyrg(i)});
    sz(end)=len(dirtyrg(i));
    As=mshape(Am,'r',num2cell(sz),'p',invperm,'struct');
end
end
%%
function A=mat_blkdiag(A)
isfull=isnumeric(A);
A=mshape(A,'state','sparse','struct');
nsz=length(A.sz);
ircA=num2cell(A2irc(A),1);
if nsz==2 && A.sz(2)>1
    A.sz=[A.sz,ones(1,3-nsz)];
    ircA=[ircA,num2cell(ones(size(A.vals,1),3-nsz),1)];
end
ircA=ircA([1,3,2,3]);
A.sz=A.sz([1,3,2,3]);
A.ir=sub2ind(A.sz,ircA{:});
A=mshape(A,'r',{A.sz(1)*A.sz(2),A.sz(3)*A.sz(4)});
if isfull
    A=mshape(A,'state','full','mat');
end
end
%%
function A=mat_diagblk(A,dims)
isfull=isnumeric(A);
if length(dims)==2
    sz=mshape(A,'size');
    dims=[dims,sqrt(prod(sz)/prod(dims))];
end
A=mshape(A,'r',num2cell(dims([1,3,2,3])),'p',[1,3,2,4],'state','sparse','struct');
ircA=num2cell(A2irc(A),1);
ilast=min(3,length(ircA));
ircA=ircA(1:ilast);
A.sz=A.sz(1:ilast);
A.ir=sub2ind(A.sz,ircA{:});
if isfull
    A=mshape(A,'state','full','mat');
end
end
