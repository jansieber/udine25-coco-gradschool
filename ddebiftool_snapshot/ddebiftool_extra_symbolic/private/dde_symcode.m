function [str,iscollected]=dde_symcode(df,x,p,dx,dp,varargin)
%% convert symbolic expressions to matlab code using matlabFunction (uses temporary files!)
%
% output is string containing function body with name chooseable by
% optional input funcname.
%
% set optional input scalar to true to return every component as separate
% function (for sys_tau and its derivatives)
%
% $Id: dde_symcode.m 309 2018-10-28 19:02:42Z jansieber $
%%
%#ok<*AGROW>
default={'splitexpr',{1:length(df{1})},'funcname','sys','keeptemp',false,'outname','out',...
    'directional_derivative',true,'multifile',false,'rename',true,'inname','in',...
    'print_progress',false};
options=dde_set_options(default,varargin,'pass_on');
%% generate code with matlabFunction
% matlabFunction can only write to files, so generate temporary files where
% we save the function then read the file as string, concatenated 
nl=dde_newline();
if options.multifile
    str={};
else
    str='';
end
iscollected=false;
folder=tempname;
gotfolder=mkdir(folder);
if ~gotfolder
    error('dde_symcode:perm','dde_symcode: could not create temp folder %s',folder);
end
if dde_isoctave()
    iscollected=true;
end
%% break up output array into sequence of outputs 
% to enable scalar expansion after function call
[vars,v_orig_names]=generate_vars(df,x,p,dx,dp,options);
%% rename inputargs if requested (default to avoid clashes)
[df,vars]=rename_vars(df,vars,options);
vnames=check_names(vars,v_orig_names,'input');
%% generate code w matlab function
strcount=0;
if dde_isoctave()
end
for i=1:length(options.splitexpr)
    dfc=cellfun(@(c)c(options.splitexpr{i}),df,'UniformOutput',false);
    %% generate code in files tempname/funcname_ind_(order-1).m and re-read into string str
    for k=1:size(dfc,2)
        dfcell=num2cell(dfc{k});
        outnames=arrayfun(@(k)sprintf('%s_%d',options.outname,k),1:length(dfcell),'uniformoutput',false);
        %%
        check_names(outnames,[vnames{:}],'output');
        fname=sprintf('%s_%d_%d',options.funcname,i,k-1);
        %if dde_isoctave()
        %    strnew=dde_octave_code(dfcell,fname,cat(2,vars{:,k}));
        %else
            % write code to temporary file
            filename=fullfile(folder,[fname,'.m']);
            if iscollected
                inputs=horzcat(vars{:,k});
                outargs={};
            else
                outargs={'outputs',outnames};
                inputs=vars(:,k);
            end
            w_orig=warning;
            warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
            matlabFunction(dfcell{:},'file',filename,'vars',inputs,outargs{:});
            warning(w_orig);
            % read code back into string str
            fid=fopen(filename,'r');
            strnew=fread(fid,inf);
            fclose(fid);
        %end
        strcount=strcount+1;
        if options.multifile
            str(strcount,1:2)={fname,char(strnew(:)')};
        else
            str=[str,nl,char(strnew(:)'),nl];
        end
        if options.print_progress>0
            fprintf('dde_symcode: %s of %d %d generated\n',fname,size(dfc,2),length(options.splitexpr));
        end
    end
end
%% remove folder (unless optionally prevented)
if ~dde_isoctave && ~options.keeptemp
    rmdir(folder,'s')
end
end
%%
function [vars,vnames]=generate_vars(df,x,p,dx,dp,options)
vars={};
if ~options.directional_derivative
    vars{1}=[x(:).',p(:).'];
    for i=2:numel(df)
        dxd=horzcat(dx{1:i-1});
        dpd=horzcat(dp{1:i-1});
        vars{i}=[x(:).',p(:).',dxd(:).',dpd(:).'];
    end
    vnames=cellfun(@(v){dde_names_from_sym(v)},vars);
else
    vars={x(:).';p(:).';dx{1}(:).';dp{1}(:).'};
    vnames=cellfun(@(v){dde_names_from_sym(v)},vars);
    vars=repmat(vars,1,numel(df));
    vnames=repmat(vnames,1,numel(df));
end
end
%% rename variables to options.inname
function [df,vars]=rename_vars(df,vars,options)
if ~options.rename
    return
end
for i=numel(df):-1:1
    for k=1:size(vars,1)
        newname=sprintf('%s_%d_n',options.inname,k);
        newvar=sym(newname,[1,numel(vars{k,i})]);
        df{i}=subs(df{i},vars{k,i},newvar);
        vars{k,i}=newvar;
    end
end
end
%% check if renamed names overlap
function vnames=check_names(vars,v_orig_names,label)
vnames=vars;
if ~isempty(vars) && ~iscell(vars{1})
    for i=numel(vars):-1:1
        vnames{i}=dde_names_from_sym(vars{i});
    end
end
if ~isempty(intersect([vnames{:}],[v_orig_names{:}]))
    error('dde_symcode:names',...
        'dde_symcode: name clash between %s names and variables',...
        label);
end
end