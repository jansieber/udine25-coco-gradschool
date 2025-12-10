function prob=loc_add_adjt(prob,seg,adj,varargin)
if ~adj.do
    return
end
if ismember(adj.type,{'cont','switch'})
    chart = coco_read_adjoint(seg, adj.run,adj.lab, 'chart');
end
switch adj.type
    case  'init'
        tl0={};
    case 'cont'
        tl0={'l0',chart.x};
    case 'switch'
        tl0={'l0',chart.x,'tl0',chart.t};
end
prob=coco_add_adjt(prob,seg,varargin{:},tl0{:});
end
