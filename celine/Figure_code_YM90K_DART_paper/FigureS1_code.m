%% Figure S1 
close all
clear all
load('FigureS1_data.mat')
%% Figure S1A - Example traces
figure;
subplot(2,1,1)
plot(processedTraces(:,1),'k')
hold on
ylim([-40, 5])
for i = 1:length(eventTimes_control)
    if eventTimes_control(i) >= 10000 && eventTimes_control(i) <= 20000  % Only add markers within the visible range
        hold on
        plot([eventTimes_control(i) eventTimes_control(i)], [-40], 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
    end
end
box off
xlim([10000,20000])
hold on
axis off
% Add scale bars
% Horizontal scale bar (200 ms)
samplesIn200ms = 2000;
x_start = 19000;  % Position near the end of the visible range
y_position = -35;
plot([x_start, x_start + samplesIn200ms], [y_position, y_position], 'k', 'LineWidth', 2)
text(x_start + samplesIn200ms/2, y_position - 3, '200 ms', 'HorizontalAlignment', 'center')

% Vertical scale bar (10 pA)
x_position = 19800;  % Position near the end of the visible range
y_start = -35;
plot([x_position, x_position], [y_start, y_start + 10], 'k', 'LineWidth', 2)
text(x_position + 200, y_start + 5, '10 pA', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')

subplot(2,1,2)
plot(processedTraces(:,2),'b')
hold on
ylim([-40, 5])
for i = 1:length(eventTimes_DART)
    if eventTimes_DART(i) >= 10000 && eventTimes_DART(i) <= 20000  % Only add markers within the visible range
        hold on
        plot([eventTimes_DART(i) eventTimes_DART(i)], [-40], 'r^', 'MarkerFaceColor', 'r', 'MarkerSize', 6)
    end
end
box off
axis off
xlim([10000,20000])
%% Figure S1B - Frequency plot
%load('FigureS1_B-C_data.mat')
combinedControlFreq = [control_data.baseline_frequency;control_data.drug_frequency]';
meanControlFreq = mean(combinedControlFreq,'omitmissing');
stdControlFreq = std(combinedControlFreq,'omitmissing');
nControCells = size(combinedControlFreq,1);
seControlFreq=stdControlFreq/nControCells;

combinedExpFreq = [DART_data.baseline_frequency;DART_data.drug_frequency]';
meanExpFreq = mean(combinedExpFreq,'omitmissing');
stdExpFreq = std(combinedExpFreq,'omitmissing');
nExperimentalCells = size(combinedExpFreq,1);
seExperimentalFreq=stdExpFreq/nExperimentalCells;


figure
subplot(1,2,1)
plot(combinedControlFreq','k','LineWidth',.25)
hold on
scatter([1,2],combinedControlFreq',20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar([1,2],meanControlFreq, seControlFreq,'.-k','MarkerSize',15)
xticks([1,2])
xticklabels({'Ctrl','+NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 30])
xlim([0.5 2.5])
ylabel('Hz')
subplot(1,2,2)
plot(combinedExpFreq','b','LineWidth',.25)
hold on
scatter([1,2],combinedExpFreq',20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",'b',"MarkerFaceAlpha",.25)
errorbar([1,2],meanExpFreq, seExperimentalFreq,'.-b','MarkerSize',15)
xticks([1,2])
xticklabels({'DART','+NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 30])
xlim([0.5 2.5])
x0=5;
y0=5;
width=1.75;
height=1.2;
set(gcf,'units','inches','position',[x0,y0,width,height])

%% Figure S1B - ANOVA and pairwise t-tests for frequency
lme_freq = fitlme(data_long, 'frequency~drug_condition*group+(1|cellID)');
anova(lme_freq)

[h,p1,~,stats1] =ttest2(control_data.baseline_frequency,DART_data.baseline_frequency);
%two-sample t-test to test whether control cells have sig different frequency from DART in the baseline (ACSF)
[h,p2,~,stats2] =ttest2(control_data.drug_frequency,DART_data.drug_frequency);
%two-sample t-test to test whether control cells have sig different frequency from DART in the drug (NBQX)

[h,p3,~,stats3] =ttest(control_data.baseline_frequency,control_data.drug_frequency); 
%paired test whether control cells have sig different frequency in baseline vs NBQX
[h,p4,~,stats4] =ttest(DART_data.baseline_frequency,DART_data.drug_frequency); 
%paired test whether DART cells have sig different frequency in baseline vs NBQX

[p1*4, p2*4, p3*4, p4*4]
%% Figure S1C - Amplitude plot
combinedControlAmp = [control_data.baseline_amplitude;control_data.drug_amplitude]';
meanControlAmp = mean(combinedControlAmp,'omitmissing');
stdControlAmp = std(combinedControlAmp,'omitmissing');
nControCells = size(combinedControlAmp,1);
seControlAmp=stdControlAmp/nControCells;

combinedExpAmp = [DART_data.baseline_amplitude;DART_data.drug_amplitude]';
meanExpAmp = mean(combinedExpAmp,'omitmissing');
stdExpAmp = std(combinedExpAmp,'omitmissing');
nExperimentalCells = size(combinedExpAmp,1);
seExperimentalAmp=stdExpAmp/nExperimentalCells;


figure
scatter([1],control_data.baseline_amplitude,20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
hold on
errorbar([1],meanControlAmp(1), seControlAmp(1),'.-k','MarkerSize',15)
scatter([2],DART_data.baseline_amplitude,20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",'b',"MarkerFaceAlpha",.25)
hold on
errorbar([2],meanExpAmp(1), seExperimentalAmp(1),'.-b','MarkerSize',15)
xticks([1,2])
xticklabels({'Ctrl','DART'})
set(gca,'TickDir','out')
box off
ylim([0 25])
xlim([0.5 2.5])
ylabel('pA')

x0=5;
y0=5;
width=1;
height=1.2;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% Figure S1C - t-test for amplitude
[h,p1,~,stats1] =ttest2(control_data.baseline_amplitude,DART_data.baseline_amplitude)

%% Figure S1E - Example traces and amplitude ratio plot

nFiles_HT = numel(group_data_HT.SOM); 
nFiles_ddHT = numel(group_data_ddHT.SOM);

group_data_HT.ratio = [group_data_HT.SOM.peakAmp]./[group_data_HT.Pyr.peakAmp];
group_data_ddHT.ratio = [group_data_ddHT.SOM.peakAmp]./[group_data_ddHT.Pyr.peakAmp];

group_data_HT.mean = mean(group_data_HT.ratio);
group_data_ddHT.mean = mean(group_data_ddHT.ratio);

group_data_HT.SEM = std(group_data_HT.ratio)/nFiles_HT;
group_data_ddHT.SEM = std(group_data_ddHT.ratio)/nFiles_ddHT;


figure; 
timeVector = make_time(group_data_ddHT.Pyr(1).subtract_Trace,20,1);

ddHT_cell = 2;
ax(1) = subplot(1,3,1); hold on;
plot(timeVector,smooth(group_data_ddHT.Pyr(ddHT_cell).subtract_Trace)); 
plot(timeVector,smooth(group_data_ddHT.SOM(ddHT_cell).subtract_Trace)); 
fix_axes(gcf,10,'time (ms)','pA'); axis square;
% ylim([-Inf 10]);

HT_cell = 4;
ax(2) = subplot(1,3,2); hold on;
plot(timeVector,smooth(group_data_HT.Pyr(HT_cell).subtract_Trace)); 
plot(timeVector,smooth(group_data_HT.SOM(HT_cell).subtract_Trace)); 
fix_axes(gcf,10,'time (ms)','pA');axis square;
% ylim([-Inf 10]);

linkaxes(ax);

% Ratio plot, ddHT vs HT 
subplot(1,3,3); hold on;
plot(1,group_data_ddHT.ratio,'o','MarkerSize',4,'Color',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor','none','MarkerSize',7);
errorbar(1,group_data_ddHT.mean,group_data_ddHT.SEM,'ko','MarkerSize',7,'MarkerFaceColor','k');
plot(2,group_data_HT.ratio,'o','MarkerSize',4,'Color',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor','none','MarkerSize',7);
errorbar(2,group_data_HT.mean,group_data_HT.SEM,'ko','MarkerSize',7,'MarkerFaceColor','k');
xticks([1 2]); xticklabels({['ddHT'];['HT']}); ylim([0 1]); xlim([0.5 2.5]);
fix_axes(gcf,10,'virus type','SST:Pyr EPSC ratio'); axis square;
%% Figure S1G - example trace
%load('FigureS1_G_data.mat')
figure; 
subplot(2,1,1)
plot(data(:,2)-mean(data(2500:2700,2),1))
hold on
plot(data(:,4)-mean(data(2500:2700,4),1))
xlim([2600 5000])
ylim([-900 100])
subplot(2,1,2)
plot(-(data(:,2)-mean(data(2500:2700,2),1))./data_resp(2,1))
hold on
plot(-(data(:,4)-mean(data(2500:2700,4),1))./data_resp(4,1))
xlim([2600 5000])
ylim([-7 1])
%%  Figure S1I - Normalized EPCS plots
% load('FigureS1_H-I_data.mat')
sz = size(data_norm2base_all);
figure;
for i = 2
    for ii = 1:2
        errorbar(1:10,mean(data_norm2base_all(match(i,ii),:,:),3),std(data_norm2base_all(match(i,ii),:,:),[],3)./sqrt(sz(3)-1))
        hold on
    end
    anova_mat = reshape(permute(data_norm2base_all([match(i,1) match(i,2)],:,:), [2 3 1]),[sz(2) sz(3)*2])';
    p_n2B = anova2(anova_mat,sz(3),'off');
    title([num2str(freq_list(match(i,1))) ' Hz; pStim = ' num2str(chop(p_n2B(1),2)) ' pDart = ' num2str(chop(p_n2B(2),2)) ' pInt = ' num2str(chop(p_n2B(3),2))])
    ylim([0 6])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
legend({'Control','YM90K-DART'})
ylabel('Norm to baseline EPSC amplitude')

figure;
for i = 2
    errorbar(1:10,mean(data_norm2con_all(i,:,:),3),std(data_norm2con_all(i,:,:),[],3)./sqrt(sz(3)-1))
    hold on
    ylim([0 1])
    xlim([0 11])
    ylabel('Norm EPSC')
    xlabel('Stimulus #')
    p_n2C = anova1(squeeze(data_norm2con_all(i,:,:))',[],'off');
    title([num2str(freq_list(match(i,1))) ' Hz; p = ' num2str(chop(p_n2C,2))])
    hold on
    plot(1:10,squeeze(data_norm2con_all(i,:,:)),'k')
end
ylabel('Norm to control EPSC amplitude')
xlabel('Stimulus #')


%% Figure S1K - histogram of capture index values

data = readtable('FigureS1_K_data.csv');

% Extract data for each drug
dart_data = data.Capture(strcmp(data.Drug, 'DART'));
peg_data = data.Capture(strcmp(data.Drug, 'PEG'));
edges = linspace(1, 2, 10); % Create 20 bins.

% Create figure
figure;
hold on;

% Create histograms with transparency
h1 = histogram(dart_data, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinEdges',edges);
h2 = histogram(peg_data, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BinEdges',edges);

% Add labels and title
xlabel('Capture Index');
ylabel('# Mice');

% Add legend
legend('DART', 'PEG', 'Location', 'northeast');

% Apply preferred figure settings
set(gca, 'TickDir', 'out');
grid off;
box off;

