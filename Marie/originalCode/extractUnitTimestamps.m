function unit = extractUnitTimestamps(UnivStruct, unit)

    strUnit = string(unit);
    unitID = [UnivStruct.unitID];
    unitID = unitID.';
    x = find(unitID == unit);
    %eval(strUnit '=GoodUnitStruct(x).timestamps');
    unit =  UnivStruct(x).timestamps;
end
    