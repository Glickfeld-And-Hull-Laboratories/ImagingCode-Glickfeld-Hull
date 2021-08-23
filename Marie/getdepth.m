function depth = getdepth(depths, descent, unit);

index = find(depths(:,1) == unit);
depth = depths(index, 2);
depth = descent - depth;
end