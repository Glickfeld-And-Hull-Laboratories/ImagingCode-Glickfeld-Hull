
function c = MapValue2Color(Value, MinMax, cmap)
%function c = MapValue2Color(Value, MinMax, cmap)
%
%maps value to rgb color
Value = Value(1);

%this supports backwards MinMax
if MinMax(1)>MinMax(2)
    cmap = cmap(end:-1:1, :);
    MinMax = MinMax(end:-1:1);
end

if Value<=MinMax(1)
    c = cmap(1,:);
elseif Value>=MinMax(2)
    c = cmap(end,:);
else
    numC = size(cmap,1);
    c = cmap(round(1 + (numC-1)*(Value-MinMax(1))/(MinMax(2)-MinMax(1))) ,:);
end

