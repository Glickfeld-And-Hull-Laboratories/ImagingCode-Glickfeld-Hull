function FirstJuice = ExtractFirstJuice(JuiceLicks);
FirstJuice = zeros(length(JuiceLicks(:,1)),1);
for i =1:length(JuiceLicks(:,1))
    licks = JuiceLicks{i,2};
    afterjuice = find((licks > JuiceLicks{i,1}), 1);
    if(isempty(afterjuice))
        FirstJuice(i)= NaN;
    else
        FirstJuice(i) = licks(afterjuice);
        
    end
end
end
