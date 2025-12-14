%% Script removing point.nvec.extracolumns from all variables
% as these take a lot of space. This uses eval(...) to check all variables
% in current workspaceand br_remove_extracolumns function.
tmp=who();
for i=1:length(tmp)
    eval([tmp{i},'=br_remove_extracolumns(',tmp{i},');']);
end
