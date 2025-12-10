bd_r_of_phi=coco_bd_table('r_of_phi');
FPs=bd_r_of_phi{strcmp('FP',bd_r_of_phi{:,'TYPE'}),'LAB'};
seglist={'u_gamma','u_plus','u_minus'};
prob=coco_prob();
for i=1:length(seglist)
    prob=ode_bvp2bvp(prob,seglist{i},'r_of_phi',FPs{end});
    [data.(seglist{i}),uidx.(seglist{i}),u0.(seglist{i})]=...
        coco_get_func_data(prob,[seglist{i},'.bvp.seg1.coll'],'data','uidx','u0');
    maps.(seglist{i})=data.(seglist{i}).coll_seg.maps;
    ind=uidx.(seglist{i});
    prob=coco_add_pars(prob,['T',seglist{i}],ind(maps.(seglist{i}).T_idx),...
        ['T',seglist{i}]);
end
