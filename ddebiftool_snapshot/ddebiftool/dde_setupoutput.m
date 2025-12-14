function out=dde_setupoutput(varargin)
%% typical output permutations for Setup... and dde_..._jac_res routines
switch varargin{1}
    case {'Setup','setup'}
        if length(varargin)==5 && varargin{5}
            % original: branch,success,fcnstruc (if outputfuncs present &
            % true
            out=varargin(1+(1:3));
        else
            out=varargin(1+[2,3,1]);
        end
    case 'jac_res'
        % arguments J,res,ieq
    switch varargin{2}
        case {'res','res+J','res+J+ieq'}
            out=varargin(2+[2,1,3]);
        case {'J','J+res','J+res+ieq'}
            out=varargin(2+(1:3));
        case 'ieq'
            out=varargin(2+[3,1,2]);
    end
    case 'jac_res:reverse'
        % arguments J,res,ieq
    switch varargin{2}
        case {'res','res+J','res+J+ieq'}
            out=varargin(2+[2,1,3]);
        case {'J','J+res','J+res+ieq'}
            out=varargin(2+(1:3));
        case 'ieq'
            out=varargin(2+[2,3,1]);
    end
    case {'perm','permute','permutation'}
        % permute outputs
        out=varargin(2+varargin{2});
end
end
