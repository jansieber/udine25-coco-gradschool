function [nunst,pout]=br_stability_monitor(funcs,branch,pt,varargin) 
b0=setfield(branch,'point',pt); %#ok<SFLD>
[b0,nunst]=br_stabl(funcs,b0,'exclude_trivial',true,'critical',true,varargin{:});
pout=b0.point;
end