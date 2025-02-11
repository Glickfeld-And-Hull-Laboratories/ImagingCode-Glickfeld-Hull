function unitNumber = GetUnitVector(unit, struct)
for i = 1:length(struct)
    if struct(i).unitID == unit
        index = i;
    end
end
unitNumber = struct(index).timestamps;
end
