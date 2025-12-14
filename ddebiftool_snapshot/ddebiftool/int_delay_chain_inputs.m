function xinput=int_delay_chain_inputs(histories,ixarg)
if isnumeric(ixarg)
    assert(~isnan(ixarg));
    ixarg={'x',ixarg,'id'};
else
    assert(iscell(ixarg));
    assert(size(ixarg,2)==2||size(ixarg,2)==3);
    if size(ixarg,2)==2
        ixarg(:,3)=repmat({'int'},size(ixarg,1),1);
    end
end
ninputs=cumsum([0,cellfun(@(x)numel(x),ixarg(:,2)).']);
xinput=NaN(1,ninputs(end));
for i=1:size(ixarg,1)
    [igxname,igxind,id_or_int]=deal(ixarg{i,:});
    xargnum=histories.names.(igxname);
    xinput(ninputs(i)+1:ninputs(i+1))=...
        histories.xvals{xargnum+1}{histories.(id_or_int)}(igxind(:)');
end
end