%setup 
clear all; clear global; close all;
load('FigureS7_C_D_data.mat')
%% Figure S7C - Time Course Entropy

% Pre- and Post- DART time course from example mouse
i3309_ITI_tc_pre = i3309_ITI_tc_pre - mean(i3309_ITI_tc_pre);
i3309_ITI_tc_post = i3309_ITI_tc_post - mean(i3309_ITI_tc_post);
ymax = max(vertcat(i3309_ITI_tc_pre,i3309_ITI_tc_post))+100;
ymin = min(vertcat(i3309_ITI_tc_pre,i3309_ITI_tc_post))-100;
f1 = figure;
f1.Name = 'i3309_ITI_TC';
sgtitle('i3309 ITI Timecourses')
hold on
subplot(2,1,1)
plot(i3309_ITI_tc_pre);
title('Control')
ylim([ymin ymax])
ylabel('F')
subplot(2,1,2)
plot(i3309_ITI_tc_post);
title('Post-DART')
ylim([ymin ymax])
ylabel('F')
xlabel('Frames')
hold off

% PV entropy calculation

nMice = size(TCs_PV,1);
mice = TCs_PV(:,2);
ent_PV = cell(nMice,2);
ent_PV(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_PV{mouse,1};
    for row = 1:length(mouse_data)
        mouse_data{row,1} = mouse_data{row,1} - mean(mouse_data{row,1});
    end
    mouse_id = TCs_PV{mouse,2};
    nSesh = size(mouse_data,1);
    ent_mat = NaN(nSesh,1);
    for sesh = 1:size(mouse_data,1)
        curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        ent_mat(sesh) = SampleEn_TH(curr_data,0.2,3);
    end
    ent_PV{mouse,1} = ent_mat;
end

% SST entropy calculation
nMice = size(TCs_SST,1);
mice = TCs_SST(:,2);
ent_SST = cell(nMice,2);
ent_SST(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_SST{mouse,1};
    for row = 1:length(mouse_data)
        mouse_data{row,1} = mouse_data{row,1} - mean(mouse_data{row,1});
    end
    mouse_id = TCs_SST{mouse,2};
    nSesh = size(mouse_data,1);
    ent_mat = NaN(nSesh,1);
    for sesh = 1:size(mouse_data,1)
        curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        ent_mat(sesh) = SampleEn_TH(curr_data,0.2,3);
    end
    ent_SST{mouse,1} = ent_mat;
end

% sample entropy comparison plots
ent_PV_mat = reshape(cell2mat(ent_PV(:,1)),2,[])';
ent_SST_mat = reshape(cell2mat(ent_SST(:,1)),2,[])';
ent_PV_DART = vertcat(ent_PV_mat(1:2,:),ent_PV_mat(4,:));

EI_PV = ent_PV_mat(:,2) ./ ent_PV_mat(:,1);
EI_PV_DART = vertcat(EI_PV(1:2,:),EI_PV(4,:));
EI_PV_PEG = vertcat(EI_PV(3,:),EI_PV(5,:));
EI_SST = ent_SST_mat(:,2) ./ ent_SST_mat(:,1);

f2 = figure;
f2.Name = 'SampleEntropy';
sgtitle('FullSessionSample Entropy');
hold on
plot(0.7,EI_PV_DART(:,1),'o','Color',"#0047AB")
plot(1.4,EI_PV_PEG(:,1),'o','Color',"#89CFF0")
plot(2.1,EI_SST(:,1),'o','Color',"#D95319")
ylim([0.5 2.5])
ylabel('Sample Entropy Post/Pre')
xlim([0 2.8])
xticks([0.7 1.4 2.1])
xticklabels({'PV+DART','PV+PEG','SST+DART'})
hold off

EI_PV_PEG(3) = NaN;
EI_anovamat = [EI_PV_DART EI_PV_PEG EI_SST];
[p,tbl,stats] = anova1(EI_anovamat);
[c,m,h,gnames] = multcompare(stats); % variable "c" contains the result of post-hoc pairwise tests for the between-DART-condition ANOVA

%% Figure S7D - Neuronal Activity Synchrony

% plots percent active cells as a measure of synchrony for the example
% mouse
i3309_iti_sync = [iti_synch_PV{1,1} iti_synch_PV{1,2}];
nDay = 2;
f3 = figure;
f3.Name = 'Synch_hist';
sgtitle('i3309 Percent Active Histogram');
hold on
subplot(2,1,1)
histogram(i3309_iti_sync(:,1),[0:0.05:1]);
xlim([0 1])
ylim([0 100])
title('Baseline')
ylabel('Frames')
subplot(2,1,2)
histogram(i3309_iti_sync(:,2),[0:0.05:1]);
xlim([0 1])
ylim([0 100])
title('Post-DART')
ylabel('Frames')
xlabel('Percent Active')
hold off


iti_synch_PV_DART = vertcat(iti_synch_PV(1:2,:),iti_synch_PV(4,:));
[nrow ncol] = size(iti_synch_PV_DART);
synch_value = nan(nrow,ncol);
for row = 1:nrow
    for col = 1:ncol
        synch_value(row,col) = mean(iti_synch_PV_DART{row,col});
    end
end

nMice_PV = size(iti_synch_PV,1);
nDay = 2;
start = 1;
SI_PV = NaN(nMice_PV,1);
for iM = 1:nMice_PV
    for iD = 1:nDay
        mean_synch(iM,iD) = mean(iti_synch_PV{iM,iD});
    end
    SI_PV(iM) = (mean_synch(iM,2)./mean_synch(iM,1));
end

nMice_SST = size(iti_synch_SST,1);
nDay = 2;
start = 1;
SI_SST = NaN(nMice_SST,1);
for iM = 1:nMice_SST
    for iD = 1:nDay
        start = start+1;
        mean_synch(iM,iD) = mean(iti_synch_SST{iM,iD});
    end
    SI_SST(iM) = (mean_synch(iM,2)./mean_synch(iM,1));
end

SI_PV_DART = vertcat(SI_PV(1:2),SI_PV(4));
SI_PV_PEG = vertcat(SI_PV(3),SI_PV(5));

% compare synchrony across all mice
f4 = figure;
f4.Name = 'Synchrony';
sgtitle('Synchrony');
hold on
plot(0.7,SI_PV_DART(:,1),'o','Color',"#0047AB")
plot(1.4,SI_PV_PEG(:,1),'o','Color',"#89CFF0")
plot(2.1,SI_SST(:,1),'o','Color',"#D95319")
ylim([0 3])
ylabel('Synchrony Index')
xlim([0 2.8])
xticks([0.7 1.4 2.1])
xticklabels({'PV+DART','PV+PEG','SST+DART'})
hold off

SI_PV_PEG(3) = NaN;

% stats
anova_mat = [SI_PV_DART SI_PV_PEG SI_SST];
[p,tbl,stats] = anova1(anova_mat);
[c,m,h,gnames] = multcompare(stats); % variable "c" contains the result of post-hoc pairwise tests for the between-DART-condition ANOVA