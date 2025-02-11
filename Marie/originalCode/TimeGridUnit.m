function TimeGridUnit = TimeGridUnit(TimeGridA, TimeGridB, unit)
TimeGridUnit = [];
for i = 1:length(TimeGridB)
    AddThis1 = unit(unit<TimeGridB(i));
    AddThis2 = AddThis1(AddThis1 > TimeGridA(i));
    TimeGridUnit = [TimeGridUnit; AddThis2];
end
end