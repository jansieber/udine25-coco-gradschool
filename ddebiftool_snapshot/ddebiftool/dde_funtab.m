function tab=dde_funtab(tab)
%% create empty table containing function structures and their references
% if tab is not present or empty
if nargin<1 || isempty(tab)
    tab=struct('ntau',0,'taumap',zeros(0,2),'xtest',zeros(0,1),'partest',zeros(1,0),'nf',0,...
        'fun_ref',repmat(dde_funformat_create(),0,1),...
        'fun_names',struct(),'fun_namelist',{{}});
end
end
