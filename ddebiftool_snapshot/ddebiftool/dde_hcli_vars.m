function [varfields,xtype,kind]=dde_hcli_vars()
varfields=struct('profile',1,'period',1,'parameter',1,...
    'x1',1,'x2',1,'v',1,'lambda_v',1,'w',1,'lambda_w',1,'alpha',1);
kind='coll';
xtype={'profile','x1','x2','v','w'};
end