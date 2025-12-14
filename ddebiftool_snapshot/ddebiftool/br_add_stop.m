function out=br_add_stop(mth,varargin)
stop=create_stop(varargin{:});
if isempty(varargin)
    stop=repmat(stop,0,1);
end
isbranch=false;
if isfield(mth,'method')
    br=mth;
    isbranch=true;
    mth=br.method.continuation;
end
mth=stop_struct(mth);    
if ~isempty(stop) && isempty(stop.name)
    stop.name=sprintf('stop%d',length(mth.stops)+1);
end
mth.stops(end+(1:length(stop)))=stop;
out=mth;
names={mth.stops.name};
[dum,ix]=unique(names); %#ok<ASGLU>
mth.stops=mth.stops(ix);
if isbranch
    br.method.continuation=mth;
    out=br;
end
end
%%
function mth=stop_struct(mth)
if ~isfield(mth,'stops')
    mth.stops=repmat(create_stop(),0,1);
    return
end
switch class(mth.stops)
    case 'struct'
        stops=mth.stops;
    case 'cell'
        for i=length(mth.stops):-1:1
            stops(i)=create_stop('online',mth.stops{i},'name',sprintf('stop%d',i));
        end
    case 'function_handle'
        stops=create_stop('detect',mth.stops,sprintf('stop%d',1));
    otherwise
        stops=repmat(create_stop(),0,1);
end
mth.stops=stops;
end
%%
function stop=create_stop(varargin)
default={'state','corrector','online',[],'final',[],'name','','argument_format','points',...
    'force',false,'order',0};
stop=dde_set_options(default,varargin,'pass_on');
end