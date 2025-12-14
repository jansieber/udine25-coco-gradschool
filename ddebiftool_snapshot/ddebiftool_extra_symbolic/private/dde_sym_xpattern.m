function xpattern=dde_sym_xpattern(fs,xxs,xpattern)
if nargin<3
    xpattern=zeros(2,0);
end
fnames=dde_names_from_sym(symvar(fs(:)));
xnames=dde_names_from_sym(xxs(:));
[dum1,dum2,ixxs]=intersect(fnames(:),xnames(:)); %#ok<ASGLU>
[ir,ic]=ind2sub(size(xxs),ixxs);
xpattern=dde_join_xpattern([ir(:),ic(:)]',xpattern);
end
