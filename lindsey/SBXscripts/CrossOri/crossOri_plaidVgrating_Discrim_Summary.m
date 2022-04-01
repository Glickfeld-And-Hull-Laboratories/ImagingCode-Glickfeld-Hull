mouse_mat = {'i484','i485','i490','i492','i784','i785'};
expt = 'plaid_test';
nm = length(mouse_mat);
gratingPctLeft_all = [];
plaidPctLeft_all = [];
gratingCI_all = [];
plaidCI_all = [];
selectGrating = [];
selectPlaid = [];

for i = 1:nm
    mouse = mouse_mat{i};
    load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '\' expt '_data.mat'])
    fprintf([mouse '\n'])
    gratingPctLeft_all = [gratingPctLeft_all gratingPctLeft];
    plaidPctLeft_all = [plaidPctLeft_all plaidPctLeft];
    gratingCI_all = cat(3,gratingCI_all, gratingCI);
    plaidCI_all = cat(3,plaidCI_all, plaidCI);
    selectGratingCI = [sqrt(sum((gratingPctLeft-gratingCI(:,1)).^2)) sqrt(sum((gratingCI(:,2)-gratingPctLeft).^2))]
    selectGrating = [selectGrating (gratingPctLeft(2) - gratingPctLeft(1)); ];
    selectPlaid = [selectPlaid (plaidPctLeft(2) - plaidPctLeft(1))];
end

figure;

for i = 1:nm
    subplot(2,2,1)
    errorbar([1 2],gratingPctLeft_all(:,i),gratingPctLeft_all(:,i)-gratingCI_all(:,1,i),gratingCI_all(:,2,i)-gratingPctLeft_all(:,i),'-ok')
    hold on
    subplot(2,2,2)
    errorbar([1 2],plaidPctLeft_all(:,i),plaidPctLeft_all(:,i)-plaidCI_all(:,1,i),plaidCI_all(:,2,i)-plaidPctLeft_all(:,i),'-ok')
    hold on
    subplot(2,2,3)
    plot([1 2],[selectGrating(i) selectPlaid(i)],'-oc') 
    hold on
end
subplot(2,2,1)
title('Grating')
xlim([0 3])
ylabel('Fraction left choice')
set(gca,'XTick',[1:2],'XTickLabels',{'0','90'})
xlabel('Direction (deg)')
subplot(2,2,2)
title('Plaid')
xlim([0 3])
ylabel('Fraction left choice')
set(gca,'XTick',[1:2],'XTickLabels',{'0','90'})
xlabel('Direction (deg)')
subplot(2,2,3)
errorbar([1 2],mean([selectGrating; selectPlaid],2), std([selectGrating; selectPlaid],[],2)./sqrt(nm),'-ok') 
xlim([0 3])
ylim([0 1])
set(gca,'XTick',[1:2],'XTickLabels',{'grating','plaid'})
ylabel('Selectivity')

print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\AllMice_PlaidTest_Summary.pdf'],'-dpdf','-fillpage');
