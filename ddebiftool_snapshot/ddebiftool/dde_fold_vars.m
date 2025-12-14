function [varfields,xtype,kind]=dde_fold_vars()
varfields=struct('x',1,'v',1,'parameter',1);
kind='stst';
xtype={'x','v'};
end