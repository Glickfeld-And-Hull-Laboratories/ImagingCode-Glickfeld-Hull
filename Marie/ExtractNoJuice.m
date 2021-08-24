function NoJuice = ExtractNoJuice(JuiceLicks);
NoJuice = [];
for i =1:length(JuiceLicks(:,1))
    licks = JuiceLicks{i,2};
    beforejuice = licks(find((licks < JuiceLicks{i,1})));
    if(isempty(beforejuice))
       
    else
        NoJuice = [NoJuice beforejuice];
        
    end
end
end
