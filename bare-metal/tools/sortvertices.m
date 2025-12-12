function vs = sortvertices(P)

temp = cell2mat(P.faceV');
vs   = temp(1,:);
lookfor   = temp(1,2);
temp(1,:) = [0, 0];
for k=2:P.nFaces-1
  i = find(lookfor==temp(:,1));
  if isempty(i)
    i = find(lookfor==temp(:,2));
    lookfor = temp(i(1),1);
  else
    lookfor = temp(i(1),2);
  end
  vs(k+1) = lookfor;
  temp(i(1),:) = [0, 0];
end

end
