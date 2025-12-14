function br_plot(branch,x_m,y_m,varargin)
%% plot branch
% function br_plot(branch,x_measure,y_measure,line_type)
% INPUT:
%	branch branch of points
%	x_measure scalar measure for x-coordinate
%	y_measure scalar measure for y_coordinate
%	line_type line type to plot with

% (c) DDE-BIFTOOL v. 2.00, 22/12/2000
%
%
if isempty(varargin)
  line_type='';
elseif ischar(varargin{1}) && ~isvarname(varargin{1})
  line_type=varargin{1};
  varargin=varargin(2:end);
end
default={'ax',gca};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ax=options.ax;
ish=ishold(ax);
if isempty(line_type)
  line_type='';
end
ll=length(branch.point);
if ll<=1
  error('BR_PLOT: ll=%d, branch contains too few points.',ll);
end
[ffx, llx] = loc_extract_from_banch(x_m, branch);
[ffy, lly] = loc_extract_from_banch(y_m, branch);
    
mx=max(llx);
my=max(lly);

if mx~=min(llx)
  for i=1:mx
    j_end=0;
    while 1
      llx(ll+1)=mx+1;
      j_start=j_end+1;
      while llx(j_start)<i
       j_start=j_start+1;
      end
      if j_start==ll+1
        break;
      end
      j_end=j_start+1;
      llx(ll+1)=0;
      while llx(j_end)>=i
        j_end=j_end+1;
      end
      j_end=j_end-1;
      plot(ax,ffx(j_start:j_end,min(i,size(ffx,2))),ffy(j_start:j_end),line_type);
      hold(ax,'on');
    end
  end
elseif my~=min(lly)
  for i=1:my
    j_end=0;
    while 1
      lly(ll+1)=my+1;
      j_start=j_end+1;
      while lly(j_start)<i
       j_start=j_start+1;
      end
      if j_start==ll+1
        break;
      end
      j_end=j_start+1;
      lly(ll+1)=0;
      while lly(j_end)>=i
        j_end=j_end+1;
      end
      j_end=j_end-1;
      plot(ax,ffx(j_start:j_end),ffy(j_start:j_end,min(i,size(ffy,2))),line_type);
      hold(ax,'on')
    end
  end
else
  plot(ax,ffx,ffy,line_type); 
  hold(ax,'on');
end
if ~ish
    hold(ax,'off');
end
end

function [ffx, llx] = loc_extract_from_banch(x_m, branch)
ll=length(branch.point);
if isstruct(x_m)
    [ffx,llx]=br_measr(branch,x_m);
elseif isempty(x_m)
    ffx=1:ll;
    llx=ones(1,ll);
elseif isa(x_m,'function_handle')
    ffx=arrayfun(x_m,branch.point).';
    llx=ones(1,ll);
elseif iscell(x_m)
    ffx=arrayfun(x_m{1},branch.point).';
    if length(x_m)>1
        llx=arrayfun(x_m{2},branch.point);
    else
        llx=ones(1,ll);
    end
else
    ffx=1:ll;
    llx=ones(1,ll);
end
end