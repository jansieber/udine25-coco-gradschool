function [funcstr,derivatives]=dde_sym2funcs(fs,xxs,ps,varargin)
%% create matlab code for r.h.s (and delays if needed) from symbolic expressions
%
%   function [funcstr,derivatives]=dde_sym2funcs(fs,xxs,ps,...) 
%
% It converts the input fs (an array of symbolc expressions) into a
% right-hand side function and its derivatives that can be used for
% DDE-Biftool. The wrapper |set_symfuncs| around |set_funcs| can be used to
% create a structure containing right-hand sides and derivatives.
%
%% Inputs:
%
% * |fs|:  n x 1 array of symbolic expressions, defining the right-hand side
% * |xxs|: n x (ntau+1) array of symbols for states, |xxs(j,k)| is the
% $x_j(t-\tau_{k-1})$, and |xxs(j,1)|=$x_j(t)$.
% * |ps|: 1 x np (or np x 1) array of symbols for parameters, |ps(j)| is
% the |j|th parameter.
%
%% Common optional name-value input pairs
%
% * |'sd_delay'| (default |sym([]))|): ntau x 1 (or 1 x ntau) array tau of
% symbolic expressions depending on symbols in |xxs| and |par|. |tau(k)| is
% delay number |k| ($\tau_{k}(x,p)$) an dmay depend on |xxs(j,l)| for
% |l=1..k|.
% * |'write'| (default |true|): write output to file?
% * |'keepnames'| (default |false|): function arguments are renamed to
% avoid bugs of matlabfunction when symbols clash with buitlins (eg beta),
% controllable with argument |'inname'| (default |'in'|, so variables are
% 'in1','in2',...)
% * |'filename'| (default |'sys'|): results are written to function file
% with this name.
% * |'multifile'| (default |false|): write each function into its own file?
% * |'folder'| (default pwd()) is used as the base folder where to put the file(s).
% * |'maxorder'| (default 5): maximal order of derivatives to be computed.
% * |'sd_delay_seq'| (default {1,...,ntau}): do delays need to be called in
% sequence? If sd_delay_seq={1,[2,3],4} then sys_tau(1,x(t),p),
% sys_tau(2:3,[x(t),x(t-tau1)],p) and
% sys_tau(4,[x(t),x(t-tau1),...x(t-tau3)],p) can be called.
%
%% Outputs (often not needed)
%
% * |fstr|: strong containing the code for writing the function.
% * |derivatives|: symbolic expressions for all dervatives (a structure,
% containing the fields |df|, |xx|, |parameter|, |dx| and |dp|. |df|
% contains the expressions, |xx|, |parameter|, |dx| and |dp| the symbols
% used in these expressions. 
%
% The routine will create a vectorized matlab
% function |y=f(xx,p)| from the expressions in |fs|, and its directional
% derivatives up to a maximum order (default 5). For state-dependent delays
% (optional argument |sd_delay| non-empty), it also computes the values and
% directional derivatives of all delays.
%% Warning
% The file will write to temporary files, since matlabFunction can only
% write to files.
%
% $Id: dde_sym2funcs.m 320 2019-02-01 00:38:43Z jansieber $
%%
default={'svn_id','','write',true,'filename','sys',{'sys_tau','sd_delay'},sym([]),'directional_derivative',true,...
    {'sys_tau_seq','sd_delay_seq'},{},'multifile',false,'folder',pwd()};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% use directional derivative if delay is state dependent
% because multilinear derivatives are not yet implemented
if ~isempty(options.sys_tau)
    options.directional_derivative=true;
end
pass_on=[{'directional_derivative',options.directional_derivative},pass_on];
%% convert possibly cell arrays safely into sym arrays
fs=dde_sym_from_cell(fs);
xxs=dde_sym_from_cell(xxs);
ps=dde_sym_from_cell(ps);
options.sys_tau=dde_sym_from_cell(options.sys_tau);
fs=fs(:);
xpattern=dde_sym_xpattern(fs,xxs);
%% extract function name from (optional) filename
[pth,funcname,ext]=fileparts(options.filename);
[fstr{1},iscollected,df{1},v{1},q{1}]=dde_symdericode(fs,xxs,ps,[funcname,'_rhs'],...
    'multifile',options.multifile,'splitexpr',{1:length(fs)},pass_on{:});
maxorder=length(df{1})-1;
npar=length(ps);
sd_delay_in=[];
nx=size(xxs,1);
nf=length(fs);
if ~isempty(options.sys_tau)
    tp_del=true;
    tau=options.sys_tau;
    xpattern=dde_sym_xpattern(tau,xxs,xpattern);
    ntau=length(tau);
    if isempty(options.sys_tau_seq)
        taugroup=num2cell(1:ntau);
        options.sys_tau_seq=taugroup;
    else
        taugroup=options.sys_tau_seq;
    end        
    sd_delay_in=dde_tauin_from_tauseq(taugroup);
    [fstr{2},iscollected,df{2},v{2},q{2}]=dde_symdericode(tau,xxs,ps,[funcname,'_tau'],...
        'splitexpr',taugroup,'multifile',options.multifile,pass_on{:});
    %[fstr{3},df{3},v{3}]=dde_sdmf_symdericode(fs,tau,xxs,ps,[funcname,'_combined'],pass_on{:});
    %mf_dxlength=numel(v{3});
else
    ntau=size(xxs,2)-1;
    tp_del=false;
    %mf_dxlength=0;
end
if ~options.multifile
    str=[fstr{:}];
else
    str=cat(1,fstr{:});
end
derivatives=struct('df',df,'xx',xxs,'parameter',ps,'dx',v,'dp',q);
%% octave's output ends with 'end', matlab's does not
if dde_test_for_end()
    function_end='end';
else
    function_end='';
end
%% create full function (still string)
nl=dde_newline();
header=sprintf('function varargout=%s(action,varargin)',funcname);
comment=[...
    '%% Automatically generated with matlabFunction',nl,...
    '% ',options.svn_id,nl,...
    '%#ok<*DEFNU,*INUSD,*INUSL>',nl];
body=[...
    'switch action',nl,...
    '  case ''ntau''',nl,...
    '   varargout{1}=',num2str(ntau),';',nl,...
    '   return',nl,...
    '  case ''npar''',nl,...
    '   varargout{1}=',num2str(npar),';',nl,...
    '   return',nl,...
    '  case ''nf''',nl,...
    '   varargout{1}=',num2str(nf),';',nl,...
    '   return',nl,...
    '  case ''nx''',nl,...
    '   varargout{1}=',num2str(nx),';',nl,...
    '   return',nl,...
    '  case ''tp_del''',nl,...
    '   varargout{1}=',num2str(tp_del),';',nl,...
    '   return',nl,...    
    '  case ''maxorder''',nl,...
    '   varargout{1}=',num2str(maxorder),';',nl,...
    '   return',nl,... 
    '  case ''iscollected''',nl,...
    '   varargout{1}=',num2str(iscollected),';',nl,...
    '   return',nl,... 
    '  case ''directional_derivative''',nl,...
    '   varargout{1}=',num2str(options.directional_derivative),';',nl,...
    '   return',nl,... 
    '  case ''sys_tau_seq''',nl,...
    '   varargout{1}=',cell2string(options.sys_tau_seq),';',nl,...
    '   return',nl,... 
    '  case ''sys_tau_in''',nl,...
    '   varargout{1}=[',num2str(sd_delay_in),'];',nl,...
    '   return',nl,... 
    '  case ''xpattern''',nl,...
    '   varargout{1}=[',num2str(xpattern(:)'),'];',nl,...
    '   return',nl,... 
    'end',nl,...
    'ind=varargin{1};',nl,...
    'order=varargin{2};',nl,...
    'nout=varargin{3};',nl,...
    'f=str2func(sprintf(''',funcname,'_%s_%d_%d'',action,ind,order));',nl,...
    'varargout=cell(nout,1);',nl,...
    '[varargout{:}]=f(varargin{4:end});',nl,...
    function_end,nl,...
    nl];
if ~options.multifile
    funcstr=[header,nl,comment,nl,body,nl,str];
else
    funcstr=cat(1,{funcname,[header,nl,comment,nl,body,nl]},str);
end
%% write to file if requested
if options.write
    if isempty(ext)
        ext='.m';
    end
    if ~exist(options.folder,'dir')
        mkdir(options.folder);
    end
    if ~options.multifile
        filename=fullfile(options.folder,pth,[funcname,ext]);
        fid=fopen(filename,'w');
        if fid<=0
            error('dde_sym2funcs:perm','dde_sym2funcs: could not create function file %s',filename);
        end
        fwrite(fid,funcstr);
        fclose(fid);
    else
        for i=1:size(funcstr,1)
            filename=fullfile(options.folder,pth,[funcstr{i,1},ext]);
            fid=fopen(filename,'w');
            if fid<=0
                error('dde_sym2funcs:perm','dde_sym2funcs: could not create function file %s',filename);
            end
            fwrite(fid,funcstr{i,2});
            fclose(fid);
        end
    end
end
end
%% test if output creates fucntions with 'end' keyword
function funend=dde_test_for_end()
x=sym('x');
p=sym('p');
dx={sym('dx')};
dp={sym('dp')};
str=strsplit(dde_symcode({x*dx+p*dp},x,p,dx,dp),dde_newline());
funend=['',str{strncmp(str,'end',3)}];
end
%% split up cell array of arrays into string
function cstr=cell2string(c)
cstr='{';
for i=1:length(c)
    cstr=[cstr,'[',num2str(c{i}),']']; %#ok<AGROW>
    if i<length(c)
        cstr=[cstr,', ']; %#ok<AGROW>
    end
end
cstr=[cstr,'}'];
end