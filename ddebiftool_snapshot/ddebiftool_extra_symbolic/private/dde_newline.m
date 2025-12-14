%% newline depending on platform
function nl=dde_newline()
if dde_isoctave()
    nl=sprintf('\n'); %#ok<SPRINTFN>
else
    nlfun=builtin('newline');
    nl=nlfun();
end
end