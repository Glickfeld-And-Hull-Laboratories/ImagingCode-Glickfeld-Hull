% mouse = 'i484';
% date = {'210502','210503','210504','210505'};
% time = {'1326', '1417', '1408', '1314'};
% tr_range = {190, 295, 275, 340};
% useFB = 0;
% expt = 'plaid_test';
%  
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
% 
% mouse = 'i485';
% %date = {'210511','210512','210513','210514','210517','210518','210519'};
% date = {'210520','210521','210524','210525','210526','210527','210528','210530','210601','210602','210603','210604'};
% %time = {'1258','1302','1258','1336','1250','1332','1312'};
% time = {'1318','1226','1203','1253','1144','1210','1110','1221','1130','1201','1150','1021'};
% %tr_range = {270,235,328,340,235,210,220};
% tr_range = {200,225,250,250,150,200,175,125,240,150,250,225};
% useFB = 1;
% expt = 'plaid_train';
%  
% mouse = 'i484';
% date = {'210513'};
% time = {'1255'};
% tr_range = {375};
% useFB = 0;
% expt = 'posttrain_plaid_test';
% 
% mouse = 'i484';
% date = {'210517','210518','210519','210520','210521','210524','210525'};
% time = {'1243','1323','1259','1314','1220','1151','1246'};
% tr_range = {350,383,200,414,325,437,310};
% useFB = 0;
% expt = 'SF_plaid_test';
 
% mouse = 'i484';
% date = {'210526','210527','210528','210530','210601','210602','210603','210604'};
% time = {'1136','1206','1102','1218','1126','1152','1152','1017',};
% tr_range = {250,287,125,275,225,210,240,125};
% useFB = 1;
% expt = 'low_SF_train';
 
% mouse = 'i484';
% date = {'210819','210820','210823','210825','210826'};
% time = {'1138','1114','1120','1000','1011'};
% tr_range = {170,200,250,324,260};
% useFB = 0;
% expt = 'type2_test';
 
% mouse = 'i485';
% date = {'210820','210823','210824','210825','210826'};
% time = {'1116','1123','1047','1003','1009'};
% tr_range = {225,263,252,278,177};
% useFB = 0;
% expt = 'type2_test';
 
% mouse = 'i490';
% date = {'210806','210809','210810','210811','210812','210813','210816','210818'};
% time = {'1000','1042','1125','1029','1032','1012','1113','1039'};
% tr_range = {358,358,359,375,371,347,358,378};
% useFB = 0;
% expt = 'plaid_test';
 
% mouse = 'i484';
% date = {'211011','211012','211013','211014','211015','211018','211019','211022','211025','211026','211027','211028','211029'};
% time = {'1053','1036','1008','1023','1015','1016','1227','0956','1022','1104','0957','1118','1027'};
% tr_range = {[],[],[],[],[],[],[],[],[],[],[],[],[]};
% useFB = 0;
% expt = 'plaid_image_FF_pt1CPD';

% mouse = 'i484';
% date = {'211101','211102','211103','211104','211105','211116','211117','211118'};
% time = {'1035','1156','1104','1135','1200','1134','1005','0947'};
% tr_range = {[],[],[],[],[],[],[],[]};
% useFB = 0;
% expt = 'plaid_image_50deg_pt1CPD';

% mouse = 'i484';
% date = {'211108','211109'};
% time = {'1034','1156'};
% tr_range = {[],[],};
% useFB = 0;
% expt = 'plaid_image_50deg_pt05CPD';

% mouse = 'i484';
% date = {'211110','211111'};
% time = {'1133','1011'};
% tr_range = {[],[],};
% useFB = 0;
% expt = 'plaid_image_FF_pt05CPD';

% mouse = 'i484';
% date = {'211110','211111'};
% time = {'1133','1011'};
% tr_range = {[],[],};
% useFB = 0;
% expt = 'plaid_image_FF_pt05CPD';

% mouse = 'i484';
% date = {'211119','211122','211123'};
% time = {'0950','1033','1134'};
% tr_range = {[],[],[]};
% useFB = 0;
% expt = 'plaid_image_FF_pt2CPD';
%  

% mouse = 'i492';
% date = {'211122','211123','211126','211129','211201','211203','211206','211207'};
% time = {'1123','1309','1146','1147','1128','1130','1142','1045'};
% tr_range = {[],[],220,[],215,100,330,420};
% useFB = 0;
% expt = 'plaid_test';
% 
% mouse = 'i784';
% date = {'211208','211210','211213','211214','211215','211216','211217','211220','211221','211222'};
% time = {'1200','1140','1122','1237','1001','1138','1059','1102','1117','1112'};
% tr_range = {180,150,[],[],240,[],170,[],190,[]};
% useFB = 0;
% expt = 'plaid_test';
% 
mouse = 'i785';
date = {'211208','211210','211213','211214','211215','211216','211217','211220','211221','211222'};
time = {'1139','1143','1124','1240','1003','1140','1101','1104','1120','1114'};
tr_range = {190,180,170,[],225,260,180,210,190,230};
useFB = 0;
expt = 'plaid_test';
%%
nd = length(date);
start = 0;
clear temp
close all
all_plaid_pctcorr = zeros(1,nd);
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
if isempty(tr_range{i})
    tr_range{i} = length(doFB);
end
if useFB
    trial_use = 1:tr_range{i};
else
    trial_use = intersect(1:tr_range{i},find(~doFB));
end
 
B2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
nb = length(unique(B2Ix));
tLeft = celleqel2mat_padded(input.tLeftResponse); 
tRight = celleqel2mat_padded(input.tRightResponse); 
tMask = celleqel2mat_padded(input.tDoMask);
if isfield (input,'doType2Mask')
    if input.doType2Mask
        tMask2 = celleqel2mat_padded(input.tDoType2Mask);
        tMask(find(tMask2)) = 0;
    else
        tMask2 = zeros(size(tMask));
    end
else
    tMask2 = zeros(size(tMask));
end
% for ib = 1:nb
% [gratingPctCorrect(ib) gratingCI(ib,:)] = binofit(length(intersect(find(B2Ix==ib-1),intersect(find(~tMask&SIx),trial_use))),sum(tMask(intersect(find(B2Ix==ib-1),intersect(find(SIx+FIx),trial_use)))==0));
% [plaidPctCorrect(ib) plaidCI(ib,:)] = binofit(length(intersect(find(B2Ix==ib-1),intersect(find(tMask&SIx),trial_use))),sum(tMask(intersect(find(B2Ix==ib-1),intersect(find(SIx+FIx),trial_use)))));
% end
 
tGratingDir = chop(celleqel2mat_padded(input.tGratingDirectionStart),3);
tGratingDir([find(tMask) find(tMask2)]) = nan;
gratingDirs = unique(tGratingDir);
gratingDirs(isnan(gratingDirs)) = [];
 
tPlaidDir = chop(celleqel2mat_padded(input.tPlaidDirectionStart),3);
tPlaidDir(find(tMask==0)) = nan;
plaidDirs = unique(tPlaidDir);
plaidDirs(isnan(plaidDirs)) = [];
 
tPlaid2Dir = chop(celleqel2mat_padded(input.tPlaidDirectionStart),3);
tPlaid2Dir(find(tMask2==0)) = nan;
plaid2Dirs = unique(tPlaid2Dir);
plaid2Dirs(isnan(plaid2Dirs)) = [];
 
tDir = zeros(size(tGratingDir));
tDir(find(tMask==0)) = tGratingDir(find(tMask==0));
tDir(find(tMask)) = tPlaidDir(find(tMask));
tDir(find(tMask2)) = tPlaid2Dir(find(tMask2));
 
tLeftChoice = calcChoice2AFC(input);
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
 
RTs = celleqel2mat_padded(input.tDecisionTimeMs);
gratingPctLeft = zeros(length(gratingDirs),1,nb);
plaidPctLeft = zeros(length(plaidDirs),1,nb);
plaid2PctLeft = zeros(length(plaid2Dirs),1,nb);
ntrials_grating = zeros(length(gratingDirs),1,nb);
ntrials_plaid = zeros(length(plaidDirs),1,nb);
ntrials_plaid2 = zeros(length(plaid2Dirs),1,nb);
gratingCI = zeros(length(gratingDirs),2,nb);
plaidCI = zeros(length(plaidDirs),2,nb);
plaid2CI = zeros(length(plaid2Dirs),2,nb);
gratingRT = zeros(length(gratingDirs),2,nb);
plaidRT = zeros(length(plaidDirs),2,nb);
plaid2RT = zeros(length(plaid2Dirs),2,nb);
for ib = 1:nb
    for ii = 1:length(gratingDirs)
        grating_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask==0 & tGratingDir == gratingDirs(ii))));
        [gratingPctLeft(ii,:,ib) gratingCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),grating_ind)),length(grating_ind));
        gratingRT(ii,1,ib) = mean(RTs(grating_ind));
        gratingRT(ii,2,ib) = std(RTs(grating_ind),[],2)./sqrt(length(grating_ind));
        ntrials_grating(ii,1,ib) = length(grating_ind);
    end
    for ii = 1:length(plaidDirs)
        plaid_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask & tPlaidDir == plaidDirs(ii))));
        [plaidPctLeft(ii,:,ib) plaidCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid_ind)),length(plaid_ind));
        plaidRT(ii,1,ib) = mean(RTs(plaid_ind));
        plaidRT(ii,2,ib) = std(RTs(plaid_ind),[],2)./sqrt(length(plaid_ind));
        ntrials_plaid(ii,1,ib) = length(plaid_ind);
    end
    if sum(tMask2)
        for ii = 1:length(plaid2Dirs)
            plaid2_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask2 & tPlaid2Dir == plaid2Dirs(ii))));
            [plaid2PctLeft(ii,:,ib) plaid2CI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid2_ind)),length(plaid2_ind));
            plaid2RT(ii,1,ib) = mean(RTs(plaid2_ind));
            plaid2RT(ii,2,ib) = std(RTs(plaid2_ind),[],2)./sqrt(length(plaid2_ind));
            ntrials_plaid2(ii,1,ib) = length(plaid2_ind);
        end
    end
end
 
% figure;
for ib = 1:nb
subplot(nd+1,2,1+start)
errorbar(gratingDirs, gratingPctLeft(:,:,ib), gratingPctLeft(:,:,ib)-gratingCI(:,1,ib), gratingCI(:,2,ib)-gratingPctLeft(:,:,ib))
hold on
errorbar(plaidDirs, plaidPctLeft(:,:,ib), plaidPctLeft(:,:,ib)-plaidCI(:,1,ib), plaidCI(:,2,ib)-plaidPctLeft(:,:,ib))
if sum(tMask2)
    errorbar(plaid2Dirs, plaid2PctLeft(:,:,ib), plaid2PctLeft(:,:,ib)-plaid2CI(:,1,ib), plaid2CI(:,2,ib)-plaid2PctLeft(:,:,ib))
end
xlim([-10 100])
ylim([0 1])
xlabel('Direction (deg)')
ylabel('Fraction left')
if start==0
    if sum(tMask2)
        legend('grating','plaid','plaid2','location','southeast')
    else
        legend('grating','plaid','location','southeast')
    end
end
pctCorr_grating = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&~tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(~tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid2 = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask2(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask2(intersect(find(B2Ix==ib-1),trial_use)));
if sum(tMask2)
    title([date{i} ' Plaid %Corr- ' num2str(chop(pctCorr_plaid,2)) ' Type2 Plaid %Corr- ' num2str(chop(pctCorr_plaid2,2))])
else
    title([date{i} ' Plaid %Corr- ' num2str(chop(pctCorr_plaid,2))])
end
start = start+1;
end
movegui('center')
all_plaid_pctcorr(1,i) = pctCorr_plaid;
all_plaid2_pctcorr(1,i) = pctCorr_plaid2;
 
 
if nb<2
subplot(nd+1,2,1+start)
errorbar(gratingDirs, gratingRT(:,1), gratingRT(:,2))
hold on
errorbar(plaidDirs, plaidRT(:,1), plaidRT(:,2))
errorbar(plaid2Dirs, plaid2RT(:,1), plaid2RT(:,2))
xlim([-10 100])
ylim([0 8000])
xlabel('Direction (deg)')
ylabel('RT (ms)')
start = start+1;
end
% sgtitle([mouse ' ' date{i}])
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
if isfield(input,'doCatchTrial')
    input = rmfield(input,'doCatchTrial');
    input = rmfield(input,'fractionCatchTrial');
    input = rmfield(input,'tDoCatchTrial');
end
if isfield(input,'firstTrConsts')
    input = rmfield(input,'firstTrConsts');
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
if isfield (input,'doType2Mask')
    if input.doType2Mask
        tMask2 = celleqel2mat_padded(input.tDoType2Mask);
        tMask(find(tMask2)) = 0;
    else
        tMask2 = zeros(size(tMask));
    end
else
    tMask2 = zeros(size(tMask));
end
% for ib = 1:nbB2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
B2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
nb = length(unique(B2Ix));
tLeft = celleqel2mat_padded(input.tLeftResponse); 
tRight = celleqel2mat_padded(input.tRightResponse); 

tGratingDir = chop(celleqel2mat_padded(input.tGratingDirectionStart),3);
tGratingDir([find(tMask) find(tMask2)]) = nan;
gratingDirs = unique(tGratingDir);
gratingDirs(isnan(gratingDirs)) = [];
 
tPlaidDir = chop(celleqel2mat_padded(input.tPlaidDirectionStart),3);
tPlaidDir(find(tMask==0)) = nan;
plaidDirs = unique(tPlaidDir);
plaidDirs(isnan(plaidDirs)) = [];
 
tPlaid2Dir = chop(celleqel2mat_padded(input.tPlaidDirectionStart),3);
tPlaid2Dir(find(tMask2==0)) = nan;
plaid2Dirs = unique(tPlaid2Dir);
plaid2Dirs(isnan(plaid2Dirs)) = [];
 
tDir = zeros(size(tGratingDir));
tDir(find(tMask==0)) = tGratingDir(find(tMask==0));
tDir(find(tMask)) = tPlaidDir(find(tMask));
tDir(find(tMask2)) = tPlaid2Dir(find(tMask2));
 
tLeftChoice = calcChoice2AFC(input);
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
 
RTs = celleqel2mat_padded(input.tDecisionTimeMs);
gratingPctLeft = zeros(length(gratingDirs),1,nb);
plaidPctLeft = zeros(length(plaidDirs),1,nb);
plaid2PctLeft = zeros(length(plaid2Dirs),1,nb);
ntrials_grating = zeros(length(gratingDirs),1,nb);
ntrials_plaid = zeros(length(plaidDirs),1,nb);
ntrials_plaid2 = zeros(length(plaid2Dirs),1,nb);
gratingCI = zeros(length(gratingDirs),2,nb);
plaidCI = zeros(length(plaidDirs),2,nb);
plaid2CI = zeros(length(plaid2Dirs),2,nb);
gratingRT = zeros(length(gratingDirs),2,nb);
plaidRT = zeros(length(plaidDirs),2,nb);
plaid2RT = zeros(length(plaid2Dirs),2,nb);
for ib = 1:nb
    for ii = 1:length(gratingDirs)
        grating_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask==0 & tGratingDir == gratingDirs(ii))));
        [gratingPctLeft(ii,:,ib) gratingCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),grating_ind)),length(grating_ind));
        gratingRT(ii,1,ib) = mean(RTs(grating_ind));
        gratingRT(ii,2,ib) = std(RTs(grating_ind),[],2)./sqrt(length(grating_ind));
        ntrials_grating(ii,1,ib) = length(grating_ind);
    end
    for ii = 1:length(plaidDirs)
        plaid_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask & tPlaidDir == plaidDirs(ii))));
        [plaidPctLeft(ii,:,ib) plaidCI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid_ind)),length(plaid_ind));
        plaidRT(ii,1,ib) = mean(RTs(plaid_ind));
        plaidRT(ii,2,ib) = std(RTs(plaid_ind),[],2)./sqrt(length(plaid_ind));
        ntrials_plaid(ii,1,ib) = length(plaid_ind);
    end
    if sum(tMask2)
        for ii = 1:length(plaid2Dirs)
            plaid2_ind = intersect(find(B2Ix==ib-1),intersect(trial_use, find(tMask2 & tPlaid2Dir == plaid2Dirs(ii))));
            [plaid2PctLeft(ii,:,ib) plaid2CI(ii,:,ib)] = binofit(length(intersect(find(tLeftChoice),plaid2_ind)),length(plaid2_ind));
            plaid2RT(ii,1,ib) = mean(RTs(plaid2_ind));
            plaid2RT(ii,2,ib) = std(RTs(plaid2_ind),[],2)./sqrt(length(plaid2_ind));
            ntrials_plaid2(ii,1,ib) = length(plaid2_ind);
        end
    end
end
 
%figure;
for ib = 1:nb
subplot(nd+1,2,1+start)
errorbar(gratingDirs, gratingPctLeft(:,:,ib), gratingPctLeft(:,:,ib)-gratingCI(:,1,ib), gratingCI(:,2,ib)-gratingPctLeft(:,:,ib))
hold on
errorbar(plaidDirs, plaidPctLeft(:,:,ib), plaidPctLeft(:,:,ib)-plaidCI(:,1,ib), plaidCI(:,2,ib)-plaidPctLeft(:,:,ib))
errorbar(plaid2Dirs, plaid2PctLeft(:,:,ib), plaid2PctLeft(:,:,ib)-plaid2CI(:,1,ib), plaid2CI(:,2,ib)-plaid2PctLeft(:,:,ib))
xlim([-10 100])
ylim([0 1])
xlabel('Direction (deg)')
ylabel('Fraction left')
if start==0
    if sum(tMask2)
        legend('grating','plaid','plaid2','location','southeast')
    else
        legend('grating','plaid','location','southeast')
    end
end
movegui('center')
pctCorr_grating = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&~tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(~tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid2 = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask2(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask2(intersect(find(B2Ix==ib-1),trial_use)));
if sum(tMask2)
    title(['%Corr: Grating- ' num2str(chop(pctCorr_grating,2)) '; Plaid- ' num2str(chop(pctCorr_plaid,2)) '; Type2Plaid- ' num2str(chop(pctCorr_plaid2,2))])
else
    title(['%Corr: Grating- ' num2str(chop(pctCorr_grating,2)) '; Plaid- ' num2str(chop(pctCorr_plaid,2))])
end    
start = start+1;
end
 
 
 
if nb<2
title(['All days- Plaid %Correct: ' num2str(chop(pctCorr_plaid,2)) 'Type2 Plaid %Correct: ' num2str(chop(pctCorr_plaid2,2))])
subplot(nd+1,2,1+start)
errorbar(gratingDirs, gratingRT(:,1), gratingRT(:,2))
hold on
errorbar(plaidDirs, plaidRT(:,1), plaidRT(:,2))
errorbar(plaid2Dirs, plaid2RT(:,1), plaid2RT(:,2))
xlim([-10 100])
ylim([0 8000])
xlabel('Direction (deg)')
ylabel('RT (ms)')
start = start+1;
end
if sum(tMask2)
    sgtitle([mouse '- ' num2str(nd) ' sessions: ' num2str(sum(ntrials_grating)) ' grating trials; ' num2str(sum(ntrials_plaid)) ' plaid trials; ' num2str(sum(ntrials_plaid2)) ' type2 plaid trials'])
else
    sgtitle([mouse '- ' num2str(nd) ' sessions: ' num2str(sum(ntrials_grating)) ' grating trials; ' num2str(sum(ntrials_plaid)) ' plaid trials; ' num2str(sum(ntrials_plaid2)) ' type2 plaid trials'])
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '_GratingPlaidComp_' expt '.pdf'],'-dpdf','-fillpage');
 
figure;
for ib = 1:nb
subplot(2,2,ib)
errorbar(gratingDirs, gratingPctLeft(:,:,ib), gratingPctLeft(:,:,ib)-gratingCI(:,1,ib), gratingCI(:,2,ib)-gratingPctLeft(:,:,ib))
hold on
errorbar(plaidDirs, plaidPctLeft(:,:,ib), plaidPctLeft(:,:,ib)-plaidCI(:,1,ib), plaidCI(:,2,ib)-plaidPctLeft(:,:,ib))
errorbar(plaid2Dirs, plaid2PctLeft(:,:,ib), plaid2PctLeft(:,:,ib)-plaid2CI(:,1,ib), plaid2CI(:,2,ib)-plaid2PctLeft(:,:,ib))
xlim([-10 100])
ylim([0 1])
xlabel('Direction (deg)')
ylabel('Fraction left')
if ib==1
    if sum(tMask2)
        legend('grating','plaid','plaid2','location','southeast')
    else
        legend('grating','plaid','location','southeast')
    end
end
movegui('center')
pctCorr_grating = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&~tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(~tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask(intersect(find(B2Ix==ib-1),trial_use)));
pctCorr_plaid2 = sum(SIx(intersect(find(B2Ix==ib-1),trial_use))&tMask2(intersect(find(B2Ix==ib-1),trial_use)))./sum(tMask2(intersect(find(B2Ix==ib-1),trial_use)));
if sum(tMask2)
    title(['%Corr: Grating- ' num2str(chop(pctCorr_grating,2)) '; Plaid- ' num2str(chop(pctCorr_plaid,2)) '; Type2 Plaid- ' num2str(chop(pctCorr_plaid2,2))])
else
    title(['%Corr: Grating- ' num2str(chop(pctCorr_grating,2)) '; Plaid- ' num2str(chop(pctCorr_plaid,2))])
end
start = start+1;
end
 
if nb<2
if sum(tMask2)
    title(['%Correct- Plaid: ' num2str(chop(pctCorr_plaid,2)) '- Type2 Plaid: ' num2str(chop(pctCorr_plaid2,2))])
else
    title(['%Correct- Plaid: ' num2str(chop(pctCorr_plaid,2))])
end
subplot(2,2,2)
errorbar(gratingDirs, gratingRT(:,1), gratingRT(:,2))
hold on
errorbar(plaidDirs, plaidRT(:,1), plaidRT(:,2))
errorbar(plaid2Dirs, plaid2RT(:,1), plaid2RT(:,2))
xlim([-10 100])
ylim([0 8000])
xlabel('Direction (deg)')
ylabel('RT (ms)')
start = start+1;
else
    subplot(2,1,2)
    plot(all_plaid_pctcorr)
    xlabel('Session')
    ylabel('Plaid %Corr')
    ylim([0 1])
end
 
 
if sum(tMask2)
    sgtitle([mouse '- ' num2str(nd) ' sessions: ' num2str(sum(ntrials_grating)) ' grating trials; ' num2str(sum(ntrials_plaid)) ' plaid trials; ' num2str(sum(ntrials_plaid2)) ' type2 plaid trials'])
else
    sgtitle([mouse '- ' num2str(nd) ' sessions: ' num2str(sum(ntrials_grating)) ' grating trials; ' num2str(sum(ntrials_plaid)) ' plaid trials'])
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '_GratingPlaidComp_' expt '_allData.pdf'],'-dpdf','-fillpage');

if ~exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse])
    mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse])
end
save(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\Behavior\CrossOri\' mouse '\' expt '_data.mat'],'gratingPctLeft', 'plaidPctLeft','plaid2PctLeft','gratingCI','plaidCI','plaid2CI','ntrials_plaid','ntrials_plaid2','ntrials_grating','input');