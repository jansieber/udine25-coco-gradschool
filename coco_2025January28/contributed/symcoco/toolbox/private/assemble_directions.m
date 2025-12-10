function [args,ndirs]=assemble_directions(src,names)
%% provide information that helps convert symbols into directional derivatives
%
% fun: function containing symbolic derivatives and format inforamtion
% names: arguments given by user to identify quested derivatives
symbolic=~iscell(src);
if symbolic
    argnames=reshape(fieldnames(src('vector')),1,[]);
else
    argnames=reshape(src,1,[]);
end
nargs=length(argnames);
cnames=cat(1,argnames,num2cell(1:nargs));
inm=struct(cnames{:});
if ischar(names)
    names={names};
end
deg=length(names);
if symbolic
    max_order=src('maxorder');
    assert(deg<=max_order,'sco_gen:degree',...
        'derivatives only computed up to order %d, but order %d requested',max_order,deg);
end
stars=strfind(names,'*');
isdirectional=~cellfun(@isempty,stars);
names(isdirectional)=cellfun(@(s,pos)s(1:pos-1),...
    names(isdirectional),stars(isdirectional),'UniformOutput',false);
unknown=setdiff(names,argnames);
assert(isempty(unknown),'sco_gen:args',...
    'sco_gen: arguments %s not known',strjoin(unknown,','));
args=num2cell(zeros(nargs,deg));
na=0;
ndirs=NaN(1,sum(isdirectional));
for i=1:deg
    if isdirectional(i)
        na=na+1;
        args{inm.(names{i}),i}=['a',num2str(na)];
        ndirs(na)=sub2ind(size(args),inm.(names{i}),i);
    else
        args{inm.(names{i}),i}='I';
    end
end
end
