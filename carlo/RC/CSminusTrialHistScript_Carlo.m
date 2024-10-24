%expt info
% mouse = '1081';
% date = '2102122';
% time = '1755';
RCExptListCarlo_Inter;

whichSession = input('Look at experiment list and give value to session that you wish to analyze: ');
frameRateHz = 30;
prewin_frames = round(1500./frameRateHz);
postwin_frames = round(3000./frameRateHz);
tt = (-prewin_frames:postwin_frames-1).*1000/frameRateHz;
%load mworks data
fName = dir(['A:\home\carlo\rawData\behavioral\data-i' expt(whichSession).mouse '-' expt(whichSession).date '*']);
load([fName.folder '\' fName.name]);
b2Ix = celleqel2mat_padded(input.tBlock2TrialNumber);
%load events tc
fData = ['A:\home\carlo\analysis\2P\' expt(whichSession).date '_img' expt(whichSession).mouse '\getTC_' expt(whichSession).run '\figDataUntrimmed_' expt(whichSession).date '_img' expt(whichSession).mouse '.mat'];
load(fData);
targetAlign_events = targetAlign_events./(1/frameRateHz);
nFrames = size(targetAlign_events,1);
nCells = size(targetAlign_events,2);
nTrials = size(targetAlign_events,3);
%correct for trial mismatch
b2Ix = b2Ix(2:end);

%find Cs+Cs+ trials
b1b1Ix = zeros(size(b2Ix));
for itrial = 2:nTrials
    if b2Ix(itrial)==0 & b2Ix(itrial-1)==0
        b1b1Ix(itrial) = 1;
    end
end
%find Cs-Cs+ trials
b2b1Ix = zeros(size(b2Ix));
for itrial = 2:nTrials
    if b2Ix(itrial)==0 & b2Ix(itrial-1)==1
        b2b1Ix(itrial) = 1;
    end
end
%find Cs+Cs- trials
b1b2Ix = zeros(size(b2Ix));
for itrial = 2:nTrials
    if b2Ix(itrial)==1 & b2Ix(itrial-1)==0
        b1b2Ix(itrial) = 1;
    end
end
%find Cs-Cs- trials
b2b2Ix = zeros(size(b2Ix));
for itrial = 2:nTrials
    if b2Ix(itrial)==1 & b2Ix(itrial-1)==1
        b2b2Ix(itrial) = 1;
    end
end
%find Cs+Cs-Cs- trials
b1b2b2Ix = zeros(size(b2Ix));
for itrial = 3:nTrials
    if b2Ix(itrial)==1 & b2Ix(itrial-1)==1 & b2Ix(itrial-2) == 0
        b1b2b2Ix(itrial) = 1;
    end
end
%find Cs-Cs-Cs- trials
b2b2b2Ix = zeros(size(b2Ix));
for itrial = 3:nTrials
    if b2Ix(itrial)==1 & b2Ix(itrial-1)==1 & b2Ix(itrial-2) == 1
        b2b2b2Ix(itrial) = 1;
    end
end


%plot all CS+ and CS-
figure;
subplot(3,1,1)
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==0)),3),2), std(mean(targetAlign_events(:,:,find(b2Ix==0)),3),[],2)./sqrt(nCells),{'k'});
hold on
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==1)),3),2), std(mean(targetAlign_events(:,:,find(b2Ix==1)),3),[],2)./sqrt(nCells),{'r'});
vline(0,'--k')
vline(770,'k')
xlim([-500 1500])
%xlim([-1000 2000])

%plot Cs+Cs-, Cs+Cs-Cs-, Cs-Cs-Cs-
% figure;
subplot(3,1,2)
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b1b2Ix)),3),2), std(mean(targetAlign_events(:,:,find(b1b2Ix)),3),[],2)./sqrt(nCells),{'k'});
hold on
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b1b2b2Ix)),3),2), std(mean(targetAlign_events(:,:,find(b1b2b2Ix)),3),[],2)./sqrt(nCells),{'r'});
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2b2b2Ix)),3),2), std(mean(targetAlign_events(:,:,find(b2b2b2Ix)),3),[],2)./sqrt(nCells),{'m'});
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
%plot without errorbar
% figure;
subplot(3,1,3)
plot(tt,mean(mean(targetAlign_events(:,:,find(b1b2Ix)),3),2),'k');
hold on
plot(tt,mean(mean(targetAlign_events(:,:,find(b1b2b2Ix)),3),2),'r');
plot(tt,mean(mean(targetAlign_events(:,:,find(b2b2b2Ix)),3),2),'m');
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
legend({['Cs+>Cs- n = ' num2str(sum(b1b2Ix))], ['Cs+>Cs->Cs- n = ' num2str(sum(b1b2b2Ix))], ['Cs->Cs->Cs- n = ' num2str(sum(b2b2b2Ix))]})
suptitle([date ' ' mouse])
%print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\for Court\' date '_' mouse '_trialHistory.pdf'],'-dpdf')

%plot all Cs+ and Cs-, Cs+Cs- and Cs-Cs-, Cs+Cs+ and Cs-Cs+ 
figure;
subplot(3,1,1)
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==0)),3),2), std(mean(targetAlign_events(:,:,find(b2Ix==0)),3),[],2)./sqrt(nCells),{'k'});
hold on
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==1)),3),2), std(mean(targetAlign_events(:,:,find(b2Ix==1)),3),[],2)./sqrt(nCells),{'r'});
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+ (black) n = ' num2str(sum(b2Ix == 0)) '; Cs- (red) n = ' num2str(sum(b2Ix))])

subplot(3,1,2)
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b1b2Ix)),3),2), std(mean(targetAlign_events(:,:,find(b1b2Ix)),3),[],2)./sqrt(nCells),{'k'});
hold on
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2b2Ix)),3),2), std(mean(targetAlign_events(:,:,find(b2b2Ix)),3),[],2)./sqrt(nCells),{'r'});
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+Cs- (black) n = ' num2str(sum(b1b2Ix)) '; Cs-Cs- (red) n = ' num2str(sum(b2b2Ix))])

subplot(3,1,3)
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b1b1Ix)),3),2), std(mean(targetAlign_events(:,:,find(b1b1Ix)),3),[],2)./sqrt(nCells),{'k'});
hold on
shadedErrorBar(tt,mean(mean(targetAlign_events(:,:,find(b2b1Ix)),3),2), std(mean(targetAlign_events(:,:,find(b2b1Ix)),3),[],2)./sqrt(nCells),{'r'});
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+Cs+ (black) n = ' num2str(sum(b1b1Ix)) '; Cs-Cs+ (red) n = ' num2str(sum(b2b1Ix))])

suptitle([date ' ' mouse])
%print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\for Court\' date '_' mouse '_trialHistory2.pdf'],'-dpdf')



figure;
subplot(3,1,1)
plot(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==0)),3),2),'k');
hold on
plot(tt,mean(mean(targetAlign_events(:,:,find(b2Ix==1)),3),2),'r');
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+ (black) n = ' num2str(sum(b2Ix == 0)) '; Cs- (red) n = ' num2str(sum(b2Ix))])

subplot(3,1,2)
plot(tt,mean(mean(targetAlign_events(:,:,find(b1b2Ix)),3),2),'k');
hold on
plot(tt,mean(mean(targetAlign_events(:,:,find(b2b2Ix)),3),2),'r');
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+Cs- (black) n = ' num2str(sum(b1b2Ix)) '; Cs-Cs- (red) n = ' num2str(sum(b2b2Ix))])

subplot(3,1,3)
plot(tt,mean(mean(targetAlign_events(:,:,find(b1b1Ix)),3),2),'k');
hold on
plot(tt,mean(mean(targetAlign_events(:,:,find(b2b1Ix)),3),2),'r');
ylim([0 6])
vline(0, '--k')
vline(770, 'k')
xlim([-500 1500])
title(['Cs+Cs+ (black) n = ' num2str(sum(b1b1Ix)) '; Cs-Cs+ (red) n = ' num2str(sum(b2b1Ix))])

suptitle([date ' ' mouse])


%}
%create matrix that gives number of consecutive Cs- trials
b2all = zeros(size(b2Ix));
t = 0;
for itrial = 1:nTrials
    if b2Ix(itrial) == 1
        t = t+1;
        b2all(itrial) = t;
    else
        t = 0;
        b2all(itrial) = t;
    end
end

%bin and plot 1xCs-, 2xCs-, 3xCs-, and >4xCs-
figure;
[n edges bin] = histcounts(b2all);
for i = 2:4
    ind = find(bin == i);
    plot(tt,mean(mean(targetAlign_events(:,:,ind),3),2));
    hold on
    ind_n(i-1) = n(i);
end
ind = find(bin > i);
ind_n(i) = sum(n(i+1:end));
plot(tt,mean(mean(targetAlign_events(:,:,ind),3),2));
%ylim([0 0.15])
vline(0, '--k')
vline(770, 'k')
xlim([-1000 2000])
legend(num2str(ind_n'))
suptitle([date ' ' mouse])
%print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\for Court\' date '_' mouse '_trialHistory3.pdf'],'-dpdf')
