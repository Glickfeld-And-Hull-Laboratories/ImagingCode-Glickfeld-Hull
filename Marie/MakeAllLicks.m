function AllLicks = MakeAllLicks(JuiceLicks)
AllLicks= [];
for i = 1:length(JuiceLicks(:,1))
    AllLicks = [AllLicks, JuiceLicks{i, 2}];
end
AllLicks = AllLicks.';
end