function branch=br_rvers(branch)

% function t_branch=br_rvers(branch)
% INPUT:
%       branch 
% OUTPUT:
%       t_branch branch with points in reversed order

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

branch.point=branch.point(length(branch.point):-1:1);
if branch.method.continuation.use_tangent
    for i=1:length(branch.point)
        if isfield(branch.point(i),'nvec') && ...
                isfield(branch.point(i).nvec,'tangent')
            branch.point(i).nvec.tangent.vector=...
                p_axpy(-1,branch.point(i).nvec.tangent.vector,[]);
        end
    end
end
end
