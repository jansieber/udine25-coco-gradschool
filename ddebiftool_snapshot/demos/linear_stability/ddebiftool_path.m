function ddebiftool_path(base,flag)
if nargin<1
    base=[pwd(),filesep,'..',filesep,'..',filesep];
end
if nargin<2
    flag='base';
end
basefolders={'ddebiftool','ddebiftool_extra_psol',...
            'ddebiftool_utilities','ddebiftool_extra_nmfm'};
switch flag
    case 'base'
        folders=basefolders;
    case 'rotsym'
        folders=[basefolders,{'ddebiftool_extra_rotsym'}];
    case 'coco'
        folders=[basefolders,{'ddebiftool_coco'}];
    case {'sym','symbolic'}
        folders={'ddebiftool','ddebiftool_extra_symbolic'};
    case 'all'
        folders=[basefolders,{'ddebiftool_extra_rotsym',...
            'ddebiftool_coco','ddebiftool_extra_symbolic'}];
end
pth=cellfun(@(f){[base,filesep,f]},folders);
addpath(pth{:});
end