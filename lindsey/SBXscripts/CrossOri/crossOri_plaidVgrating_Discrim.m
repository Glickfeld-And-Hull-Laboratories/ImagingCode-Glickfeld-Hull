% mouse = 'i484';
% date = {'210502','210503','210504','210505'};
% time = {'1326', '1417', '1408', '1314'};
% tr_range = {190, 295, 275, 340};
% useFB = 0;
% expt = 'plaid_test';

% mouse = 'i485';
% date = {'210506','210507','210509','210510'};
% time = {'1320','1325',{'1357', '1441', '1451','1521'},'1301'};
% tr_range = {150,200,296,200};
% useFB = 0;
% expt = 'plaid_test';
% 
% mouse = 'i484';
% date = {'210506','210507','210509','210510','210511','210512'};
% time = {'1321','1326','1355','1258','1255','1258'};
% tr_range = {444, 290,350,300,290,240};
% useFB = 1;
% expt = 'plaid_train';

mouse = 'i485';
date = {'210511','210512','210513','210514','210517','210518'};
time = {'1258','1302','1258','1336','1250','1332'};
tr_range = {270,235,328,340,235,210};
useFB = 1;
expt = 'plaid_train';
 
% mouse = 'i484';
% date = {'210513'};
% time = {'1255'};
% tr_range = {375};
% useFB = 0;
% expt = 'posttrain_plaid_test';

% mouse = 'i484';
% date = {'210517','210518'};
% time = {'1243','1323'};
% tr_range = {350,383};
% useFB = 0;
% expt = 'SF_plaid_test';


%%
nd = length(date);
start = 0;
clear temp
close all
for i = 1:nd
    clear input
    if iscell(time{i})
        clear tt
        for ii = 1:length(time{i})
            t = time{i}{ii};
            load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date{i} '-' t '.mat'])
            tt(ii) = input;
        end
        input = concatenateStructuresLG(temp);
    else
        load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date{i} '-' time{i} '.mat'])
    end

SIx = strcmp(input.trialOutcomeCell,'success');
FIx = strcmp(input.trialOutcomeCell,'incorrect');
MIx = strcmp(input.trialOutcomeCell,'ignore');

doFB = celleqel2mat_padded(input.tDoFeedbackMotion);
if useFB
    trial_use = 1:tr_range{i};
else
    trial_use = intersect(1:tr_range{i},find(~doFB));
end

B2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
nb = length(unique(B2Ix));

tMask = celleqel2mat_padded(input.tDoMask);
% for ib = 1:nb
% [gratingPctCorrect(ib) gratingCI(ib,:)] = binofit(length(intersect(find(B2Ix==ib-1),intersect(find(~tMask&SIx),trial_use))),sum(tMask(intersect(find(B2Ix==ib-1),intersect(find(SIx+FIx),trial_use)))==0));
% [plaidPctCorrect(ib) plaidCI(ib,:)] = binofit(length(intersect(find(B2Ix==ib-1),intersect(find(tMask&SIx),trial_use))),sum(tMask(intersect(find(B2Ix==ib-1),intersect(find(SIx+FIx),trial_use)))));
% end
tGratingDir = celleqel2mat_padded(input.tGratingDirectionStart);
gratingDirs = unique(tGratingDir);

tPlaidDir = celleqel2mat_padded(input.tPlaidDirectionStart);
plaidDirs = unique(tPlaidDir);

tDir = zeros(size(tGratingDir));
tDir(find(tMask==0)) = tGratingDir(find(tMask==0));
tDir(find(tMask)) = tPlaidDir(find(tMask));

tLeftChoice = (SIx & tDir == 90) + (FIx & tDir == 0);
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);

RTs = celleqel2mat_padded(input.tDecisionTimeMs);
gratingPctLeft = zeros(2,1,nb);
plaidPctLeft = zeros(2,1,nb);
ntrials_grating = zeros(2,1,nb);
ntrials_plaid = zeros(2,1,nb);
gratingCI = zeros(2,2,nb);
plaidCI = zeros(2,2,nb);
gratingRT = zeros(2,2,nb);
plaidRT = zeros(2,2,nb);
for ii = 1:2
    for ib = 1:nb
        grating_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask==0 & tLeftTrial==ii-1)));
        plaid_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask & tLeftTrial==ii-1)));
        [gratingPctLeft(ii,:,ib) gratingCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),grating_ind)),length(grating_ind));
        [plaidPctLeft(ii,:,ib) plaidCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid_ind)),length(plaid_ind));
        gratingRT(ii,1,ib) = mean(RTs(grating_ind));
        gratingRT(ii,2,ib) = std(RTs(grating_ind),[],2)./sqrt(length(grating_ind));
        plaidRT(ii,1,ib) = mean(RTs(plaid_ind));
        plaidRT(ii,2,ib) = std(RTs(plaid_ind),[],2)./sqrt(length(plaid_ind));
        ntrials_grating(ii,1,ib) = length(grating_ind);
        ntrials_plaid(ii,1,ib) = length(plaid_ind);
    end
end
% figure;
for ib = 1:nb
subplot(nd+1,2,1+start)
errorbar([0 90], gratingPctLeft(:,:,ib), gratingPctLeft(:,:,ib)-gratingCI(:,1,ib), gratingCI(:,2,ib)-gratingPctLeft(:,:,ib))
hold on
errorbar([0 90], plaidPctLeft(:,:,ib), plaidPctLeft(:,:,ib)-plaidCI(:,1,ib), plaidCI(:,2,ib)-plaidPctLeft(:,:,ib))
xlim([-10 100])
ylim([0 1])
xlabel('Direction (deg)')
ylabel('Fraction left')
if start==0
legend('grating','plaid','location','southeast')
end
pctCorr_grating = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&~tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(~tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask(intersect(find(B2Ix==ib-1),trial_use)));
title([date{i} ' Plaid %Corr- ' num2str(chop(pctCorr_plaid,2))])
start = start+1;
end
movegui('center')



if nb<2
subplot(nd+1,2,1+start)
errorbar([0 90], gratingRT(:,1), gratingRT(:,2))
hold on
errorbar([0 90], plaidRT(:,1), plaidRT(:,2))
xlim([-10 100])
ylim([0 8000])
xlabel('Direction (deg)')
ylabel('RT (ms)')
start = start+1;
end
% suptitle([mouse ' ' date{i}])
% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '_' date{i} '_GratingPlaidComp_' expt '.pdf'],'-dpdf','-fillpage');
if i>1
    [isDiff onlyIn1 onlyIn2] = compareStructures(temp(1),input);
    if isDiff & length(onlyIn1)
        for ifield = 1:length(onlyIn1)
        	input = rmfield(input,onlyIn1{ifield});
        end
    end
    if isDiff & length(onlyIn2)
        for ifield = 1:length(onlyIn2)
        	input = rmfield(input,onlyIn2{ifield});
        end
    end  
end
temp(i) = trialChopper(input,[1 tr_range{i}]);

end

%%
input = concatenateStructuresLG(temp);
SIx = strcmp(input.trialOutcomeCell,'success');
FIx = strcmp(input.trialOutcomeCell,'incorrect');
MIx = strcmp(input.trialOutcomeCell,'ignore');

doFB = celleqel2mat_padded(input.tDoFeedbackMotion);
if useFB
    trial_use = 1:length(SIx);
else
    trial_use = find(~doFB);
end
tMask = celleqel2mat_padded(input.tDoMask);
B2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
nb = length(unique(B2Ix));

tGratingDir = celleqel2mat_padded(input.tGratingDirectionStart);
gratingDirs = unique(tGratingDir);

tPlaidDir = celleqel2mat_padded(input.tPlaidDirectionStart);
plaidDirs = unique(tPlaidDir);

tDir = zeros(size(tGratingDir));
tDir(find(tMask==0)) = tGratingDir(find(tMask==0));
tDir(find(tMask)) = tPlaidDir(find(tMask));

tLeftChoice = (SIx & tDir == 90) + (FIx & tDir == 0);
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);

RTs = celleqel2mat_padded(input.tDecisionTimeMs);
gratingPctLeft = zeros(2,1,nb);
plaidPctLeft = zeros(2,1,nb);
gratingCI = zeros(2,2,nb);
plaidCI = zeros(2,2,nb);
gratingRT = zeros(2,2,nb);
plaidRT = zeros(2,2,nb);
for ii = 1:2
    for ib = 1:nb
        grating_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask==0 & tLeftTrial==ii-1)));
        plaid_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask & tLeftTrial==ii-1)));
        [gratingPctLeft(ii,:,ib) gratingCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),grating_ind)),length(grating_ind));
        [plaidPctLeft(ii,:,ib) plaidCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid_ind)),length(plaid_ind));
        gratingRT(ii,1,ib) = mean(RTs(grating_ind));
        gratingRT(ii,2,ib) = std(RTs(grating_ind),[],2)./sqrt(length(grating_ind));
        plaidRT(ii,1,ib) = mean(RTs(plaid_ind));
        plaidRT(ii,2,ib) = std(RTs(plaid_ind),[],2)./sqrt(length(plaid_ind));
        ntrials_grating(ii,1,ib) = length(grating_ind);
        ntrials_plaid(ii,1,ib) = length(plaid_ind);
    end
end
%figure;
for ib = 1:nb
subplot(nd+1,2,1+start)
errorbar([0 90], gratingPctLeft(:,:,ib), gratingPctLeft(:,:,ib)-gratingCI(:,1,ib), gratingCI(:,2,ib)-gratingPctLeft(:,:,ib))
hold on
errorbar([0 90], plaidPctLeft(:,:,ib), plaidPctLeft(:,:,ib)-plaidCI(:,1,ib), plaidCI(:,2,ib)-plaidPctLeft(:,:,ib))
xlim([-10 100])
ylim([0 1])
xlabel('Direction (deg)')
ylabel('Fraction left')
if start==0
legend('grating','plaid','location','southeast')
end
movegui('center')
pctCorr_grating = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&~tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(~tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask(intersect(find(B2Ix==ib-1),trial_use)));
title(['%Corr: Grating- ' num2str(chop(pctCorr_grating,2)) '; Plaid- ' num2str(chop(pctCorr_plaid,2))])
start = start+1;
end



if nb<2
title([date{i} '- Plaid %Correct: ' num2str(chop(pctCorr_plaid,2))])
subplot(nd+1,2,1+start)
errorbar([0 90], gratingRT(:,1), gratingRT(:,2))
hold on
errorbar([0 90], plaidRT(:,1), plaidRT(:,2))
xlim([-10 100])
ylim([0 8000])
xlabel('Direction (deg)')
ylabel('RT (ms)')
start = start+1;
end

suptitle([mouse '- ' num2str(nd) ' sessions: ' num2str(sum(ntrials_grating)) ' grating trials; ' num2str(sum(ntrials_plaid)) ' plaid trials'])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '_GratingPlaidComp_' expt '.pdf'],'-dpdf','-fillpage');
