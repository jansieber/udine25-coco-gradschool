function [eigval,eigprofile]=mult_crit(funcs,point,method,remove,closest,mult)
%% find Floquet multiplier closest to unit circle and its mode
% (wrapper around dde_psol_eig kept for backward compatibility)
% function eig_v=mult_crit(period,point)
% INPUT: 
%   funcs problem functions
%   point periodic solution structure
%   method method parameters for computing stability
%   remove (optional, default=1), array of trivial Floquet multipliers to
%         remove
%   closest ([]) choose mult eigenpairs closest to value
%   mult (1) select mult eigenpairs
% OUTPUT:
%   eigval: critical non-trivial eigenvalue (critical=closest to unit
%   circle)
%	eigprofile: profile of eigenvector (mesh equals point.mesh)
%
% $Id: mult_crit.m 309 2018-10-28 19:02:42Z jansieber $
%

%% initialise arguments and functions
if nargin<6
    mult=1;
end
if nargin<5
    closest=[]; 
end
if nargin<4
    remove=[];
end
stability=dde_psol_eig(funcs,point,method,'closest',closest,'geteigenfuncs',true);
if ~isempty(remove)
    iremove=dde_match_complex(remove(:),stability.mu(:));
    stability.mu(iremove)=[];
    stability.eigenfuncs(iremove)=[];
end
%% find critical eigenvalue
if isempty(closest)
    [dum,ix]=sort(abs(log(abs(stability.mu)))); %#ok<ASGLU>
else
    [dum,ix]=sort(abs(stability.mu-closest)); %#ok<ASGLU>
end
eigval=stability.mu(ix(1:mult));
eigprofile=stability.eigenfuncs(ix(1:mult)).profile;
end
