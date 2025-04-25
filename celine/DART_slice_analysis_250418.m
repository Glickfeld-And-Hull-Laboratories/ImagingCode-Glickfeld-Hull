% this script goes through a list of files for control cells and a list for
% DART cells, finding the mean EPCS frequency in baseline and drug
% sweeps. This is based on the output of the EPSC_detect_fxn script.

close all;clear all
cd('\home\celine\Analysis\ePhys\YM90K_DART')
d=string(datetime('today'));
outpath = fullfile('\home\celine\Analysis\ePhys\YM90K_DART',d);
mkdir(outpath);

plottingEnabled = 0;

file_config = {
    % Control group
    {'23615001_epsc_results.mat',[1:30],[45:70],[],'control'},
    {'23616000_epsc_results.mat',[1:29],[47:70],[],'control'},
    {'23616001_epsc_results.mat',[1:47],[62:87],[],'control'},
    {'23616003_epsc_results.mat',[1:40],[48:62],[],'control'},
    {'23621001_epsc_results.mat',[1:20],[35:56],[],'control'},
    {'23621002_epsc_results.mat',[1:21],[],[],'control'},
    {'23724000_epsc_results.mat',[3:25],[40:54],[],'control'},
    {'24110000_epsc_results.mat',[1:53],[76:93],[],'control'},
    {'24104000_epsc_results.mat',[1:65],[89:125],[],'control'},
    
    % Experimental group
    {'23623002_epsc_results.mat',[2:22],[],[42],'DART'},
    {'23623003_epsc_results.mat',[1:20],[35:58],[41,46,48],'DART'},
    {'23719004_epsc_results.mat',[1:24],[39:50],[],'DART'},
    {'23720003_epsc_results.mat',[1:23],[],[],'DART'},
    {'23724003_epsc_results.mat',[1:23],[38:59],[],'DART'},
    {'23615003_epsc_results.mat',[1:20],[37:60],[],'DART'},
    %{'24117001_epsc_results.mat',[1:35],[50:65],'DART'},
    {'24117003_epsc_results.mat',[1:73],[],[],'DART'},
    {'24115001_epsc_results.mat',[1:48],[72:85],[79],'DART'}

};


% Initialize data structures
control_data = struct();
control_data.cell_ids = {};
control_data.baseline_frequency = [];
control_data.drug_frequency = [];
control_data.baseline_amplitude = [];
control_data.drug_amplitude = [];
control_data.all_baseline_amps = [];
control_data.total_baseline_durations = 0;
control_data.all_drug_amps=[];    


DART_data = struct();
DART_data.cell_ids = {};
DART_data.baseline_frequency = [];
DART_data.drug_frequency = [];
DART_data.baseline_amplitude = [];
DART_data.drug_amplitude = [];
DART_data.all_baseline_amps = [];
DART_data.total_baseline_durations =0;
DART_data.all_drug_amps=[]; 

for iFile = 1:length(file_config)
    
    config = file_config{iFile};
    filename = config{1};
    baseline_start = config{2}(1);
    
    reject_sweeps = config{4};
    group = config{5};
    cell_id = char(filename(1:8));

    
    load(filename); %load the current file
    duration_temp = allResults.effectiveDuration;
    duration = duration_temp(1);
    baseline_freq_temp = [];
    baseline_amp_temp = [];
    all_baseline_amps_temp = [];

    if duration > 5
        baseline_sweeps = [baseline_start:baseline_start+10];
    else
        baseline_sweeps = [baseline_start:baseline_start+20];
    end

    for iSweep = baseline_sweeps
        if allResults(iSweep).baselineStable && ~ismember(iSweep, reject_sweeps) %only use good sweeps
            baseline_amp_temp = [baseline_amp_temp, allResults(iSweep).stats.meanAmplitude]; %append the mean amplitude for events in this sweep
            baseline_freq_temp = [baseline_freq_temp, allResults(iSweep).stats.frequency]; %append the frequency for events in this sweep
            all_baseline_amps_temp = [all_baseline_amps_temp,allResults(iSweep).epscs.amplitude];%append all the amplitudes for all the events in this sweep 
            all_baseline_durations_temp = length(baseline_sweeps)*duration;
        end
    end

    if ~isempty(config{3})
        drug_start = config{3}(1);
        drug_freq_temp = [];
        drug_amp_temp = [];
        all_drug_amps_temp = [];

        if duration > 5
            drug_sweeps = [drug_start:drug_start+10];
        else
            drug_sweeps = [drug_start:drug_start+20];
        end
        for iSweep = drug_sweeps
            if allResults(iSweep).baselineStable  && ~ismember(iSweep, reject_sweeps)%only use good sweeps
                drug_amp_temp = [drug_amp_temp, allResults(iSweep).stats.meanAmplitude];
                drug_freq_temp = [drug_freq_temp, allResults(iSweep).stats.frequency];
                all_drug_amps_temp = [all_drug_amps_temp, allResults(iSweep).epscs.amplitude];%append all the amplitudes for all the events in this sween 
            end
        end
    else
        
        drug_freq_temp = NaN;
        drug_amp_temp = NaN;
        all_drug_amps_temp = NaN;
    end
    
if strcmpi(group, 'control')
            control_data.cell_ids{end+1} = cell_id;
            control_data.baseline_frequency(end+1) = mean(baseline_freq_temp, 'omitmissing');
            control_data.drug_frequency(end+1) = mean(drug_freq_temp, 'omitmissing');
            control_data.baseline_amplitude(end+1) = mean(baseline_amp_temp, 'omitmissing');
            control_data.drug_amplitude(end+1) = mean(drug_amp_temp, 'omitmissing');
            control_data.all_baseline_amps = [control_data.all_baseline_amps, all_baseline_amps_temp];
            control_data.total_baseline_durations = control_data.total_baseline_durations+all_baseline_durations_temp;
            control_data.all_drug_amps = [control_data.all_drug_amps, all_drug_amps_temp];  
elseif strcmpi(group, 'DART')
            DART_data.cell_ids{end+1} = cell_id;
            DART_data.baseline_frequency(end+1) = mean(baseline_freq_temp, 'omitmissing');
            DART_data.drug_frequency(end+1) = mean(drug_freq_temp, 'omitmissing');
            DART_data.baseline_amplitude(end+1) = mean(baseline_amp_temp, 'omitmissing');
            DART_data.drug_amplitude(end+1) = mean(drug_amp_temp, 'omitmissing');
            DART_data.all_baseline_amps = [DART_data.all_baseline_amps, all_baseline_amps_temp];
            DART_data.total_baseline_durations = DART_data.total_baseline_durations+all_baseline_durations_temp;
            DART_data.all_drug_amps = [DART_data.all_drug_amps, all_drug_amps_temp];    
else
    error('neither control nor experimetnal group')
end

% Extract frequency values from the allResults structure
sweepNumbers = 1:length(allResults);
frequencies = zeros(size(sweepNumbers));

for iSweep = 1:length(allResults)
    if  ismember(iSweep, reject_sweeps)
        frequencies(iSweep) = NaN;
    else
        frequencies(iSweep) = allResults(iSweep).stats.frequency;
    end
end
if plottingEnabled
    % Create the plot
    figure;
    plot(sweepNumbers, frequencies, 'o-', 'LineWidth', 1.5);
    ylim([0 30])
    vline(baseline_sweeps(1),'k')
    vline(baseline_sweeps(end),'k')
    if ~isempty(config{3})
        vline(drug_sweeps(1))
        vline(drug_sweeps(end))
    end
    % Set labels and title
    xlabel('Sweep Number');
    ylabel('Frequency (Hz)');
    
    % Apply your preferred plotting style
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    
    % Optional: Add additional formatting if needed
    set(gca, 'FontSize', 12);
    
    sgtitle(['cell ' cell_id, ', ', group])
    print(fullfile(outpath,[cell_id, '_freqVsSweep.pdf']),'-dpdf');
end
clear baseline_freq_temp drug_freq_temp baseline_amp_temp drug_amp_temp 
clear allResults all_drug_amps_temp all_baseline_amps_temp cell_id filename
clear group overallFrequency sweepSummary reject_sweeps params baseline_start
clear drug_start duration_temp duration sweepNumbers


end
clear iSweep iFile drug_sweeps baseline_sweeps plottingEnabled d
%% Distributions of amplitudes
hypothetical_DART_amp = control_data.all_baseline_amps * .4;
%Plot normalized distributions for control and DART. 
figure;
histogram(control_data.all_baseline_amps,'BinWidth',4,'DisplayName', 'Control', 'FaceColor','k')
hold on
histogram(hypothetical_DART_amp,'BinWidth',4,'DisplayName', 'DART', 'FaceAlpha', 0.25,'FaceColor','b')
box off
set(gca,'TickDir','out')
xlabel('Amplitude (pA)')
vline(5)
print(fullfile(outpath,'hypothetical_amp_distribution.pdf'),'-dpdf');
%% Thresholding effect on frequency
control_all_frequency = length(control_data.all_baseline_amps)./control_data.total_baseline_durations;
hypothetical_DART_freq = length(find(hypothetical_DART_amp>5))./control_data.total_baseline_durations;
DART_all_frequency = length(DART_data.all_baseline_amps)./DART_data.total_baseline_durations;
table1 = table(control_all_frequency,hypothetical_DART_freq,DART_all_frequency,'VariableNames',{'Control total freq','Hypothetical thresholded DART freq','Observed total DART freq'})
% Option 2: Save as CSV file
writetable(table1, 'hypothetical_frequency_data.csv');
%% plotting frequency

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
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 25])
xlim([0.5 2.5])
ylabel('Hz')
subplot(1,2,2)
plot(combinedExpFreq','b','LineWidth',.25)
hold on
scatter([1,2],combinedExpFreq',20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",'b',"MarkerFaceAlpha",.25)
errorbar([1,2],meanExpFreq, seExperimentalFreq,'.-b','MarkerSize',15)
xticks([1,2])
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 25])
xlim([0.5 2.5])
x0=5;
y0=5;
width=2;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])


print(fullfile(outpath,'frequency.pdf'),'-dpdf');
%% plotting amplitude

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
subplot(1,2,1)
plot(combinedControlAmp','k','LineWidth',.25)
hold on
scatter([1,2],combinedControlAmp',20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",[0 0 0],"MarkerFaceAlpha",.25)
errorbar([1,2],meanControlAmp, seControlAmp,'.-k','MarkerSize',15)
xticks([1,2])
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 25])
xlim([0.5 2.5])
ylabel('pA')
subplot(1,2,2)
plot(combinedExpAmp','b','LineWidth',.25)
hold on
scatter([1,2],combinedExpAmp',20, ...
    "MarkerEdgeColor","none","MarkerFaceColor",'b',"MarkerFaceAlpha",.25)
errorbar([1,2],meanExpAmp, seExperimentalAmp,'.-b','MarkerSize',15)
xticks([1,2])
xticklabels({'Baseline','NBQX'})
set(gca,'TickDir','out')
box off
ylim([0 25])
xlim([0.5 2.5])
x0=5;
y0=5;
width=2;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])

print(fullfile(outpath,'amplitude.pdf'),'-dpdf');

%% Preparing data structures for ANOVAs

% Step 1: Prepare the data in the right format
% Make sure all data is in column vectors
control_baseline_amp = control_data.baseline_amplitude(:);
control_drug_amp = control_data.drug_amplitude(:);
dart_baseline_amp = DART_data.baseline_amplitude(:);
dart_drug_amp = DART_data.drug_amplitude(:);

% Add frequency data
control_baseline_freq = control_data.baseline_frequency(:);
control_drug_freq = control_data.drug_frequency(:);
dart_baseline_freq = DART_data.baseline_frequency(:);
dart_drug_freq = DART_data.drug_frequency(:);

% Step 2: Create a table in wide format 
n_control = length(control_baseline_amp);
n_dart = length(dart_baseline_amp);

% Create data table in wide format
cellID = [(1:n_control)'; (n_control+1:n_control+n_dart)'];
group = [repmat({'Control'}, n_control, 1); repmat({'DART'}, n_dart, 1)];
baseline_amp = [control_baseline_amp; dart_baseline_amp];
drug_amp = [control_drug_amp; dart_drug_amp];
baseline_freq = [control_baseline_freq; dart_baseline_freq];
drug_freq = [control_drug_freq; dart_drug_freq];

% Create the wide-format data table
data_wide = table(cellID, group, baseline_amp, drug_amp, baseline_freq, drug_freq);

% Convert categorical variables
data_wide.group = categorical(data_wide.group);

% Step 3: Create a long format table using stack
% First, stack the amplitude variables
data_long_amp = stack(data_wide, {'baseline_amp', 'drug_amp'}, ...
    'NewDataVariableName', 'amplitude', 'IndexVariableName', 'drug_condition');

% Stack the frequency variables
data_long_freq = stack(data_wide, {'baseline_freq', 'drug_freq'}, ...
    'NewDataVariableName', 'frequency', 'IndexVariableName', 'drug_condition');

% Combine amplitude and frequency
data_long_freq.amplitude = data_long_amp.amplitude;

% Keep only the necessary columns in the final long-format table
data_long = data_long_freq(:, {'cellID', 'group', 'drug_condition', 'amplitude', 'frequency'});

% Clean up variable names and values
data_long.drug_condition = categorical(data_long.drug_condition);

% Rename the categories for clarity
data_long.drug_condition = renamecats(data_long.drug_condition, ...
    {'baseline_freq', 'drug_freq'}, {'Baseline', 'Drug'});

% Display the first few rows to verify structure
disp('Long format data table (first few rows):');
disp(head(data_long));
%% Two-way ANOVAs for frequency and amplitude

lme_freq = fitlme(data_long, 'frequency~drug_condition*group+(1|cellID)');
anova(lme_freq)

lme_amp = fitlme(data_long, 'amplitude~drug_condition*group+(1|cellID)');
anova(lme_amp)
%% t-tests for DART and NBQX for frequency 
[h,p1,~,stats1] =ttest2(control_data.baseline_frequency,DART_data.baseline_frequency);
%two-sample t-test to test whether control cells have sig different frequency from DART in the baseline (ACSF)
[h,p2,~,stats2] =ttest2(control_data.drug_frequency,DART_data.drug_frequency);
%two-sample t-test to test whether control cells have sig different frequency from DART in the drug (NBQX)

[h,p3,~,stats3] =ttest(control_data.baseline_frequency,control_data.drug_frequency); 
%paired test whether control cells have sig different frequency in baseline vs NBQX
[h,p4,~,stats4] =ttest(DART_data.baseline_frequency,DART_data.drug_frequency); 
%paired test whether DART cells have sig different frequency in baseline vs NBQX

[p1*4, p2*4, p3*4, p4*4]

%% t-tests for DART and NBQX for amplitude 
[h,p1,~,stats1] =ttest2(control_data.baseline_amplitude,DART_data.baseline_amplitude);
%two-sample t-test to test whether control cells have sig different amplitude from DART in the baseline (ACSF)
[h,p2,~,stats2] =ttest2(control_data.drug_amplitude,DART_data.drug_amplitude);
%two-sample t-test to test whether control cells have sig different amplitude from DART in the drug (NBQX)

[h,p3,~,stats3] =ttest(control_data.baseline_amplitude,control_data.drug_amplitude); 
%paired test whether control cells have sig different amplitude in baseline vs NBQX
[h,p4,~,stats4] =ttest(DART_data.baseline_amplitude,DART_data.drug_amplitude); 
%paired test whether DART cells have sig different amplitude in baseline vs NBQX

[p1*4, p2*4, p3*4, p4*4]
%% example trace

CntrlExFile = 3
CntrlbaselineSweep=1;

config = file_config{CntrlExFile};
filename = config{1};
baseline_sweeps = config{2};
drug_sweeps = config{3};
group = config{4};
cell_id = char(filename(1:8));


load(filename); %load the current results file
abf_path = fullfile('\home\celine\Data\patch_data',[cell_id, '.abf']);
[data, si, h] = abfload(abf_path);
samplingInterval = si/1000; % Convert to ms
samplesPerMs = 1/samplingInterval;
timeVector = (0:size(data,1)-1)' * samplingInterval;
Fs = samplesPerMs*1000;

baselineEventsTimes_control = [];
baselineEventsTimes_control = [baselineEventsTimes_control, allResults(CntrlbaselineSweep).epscs.time];
eventTimes_control = baselineEventsTimes_control*samplesPerMs;

%collect the two example traces
tracesRaw = [data(:,1,CntrlbaselineSweep)];

if h.recChUnits{1,1}=='nA'
    tracesRaw = tracesRaw*1000;
end

%process the example traces 

cutoffFreq = 1300; % Cutoff frequency for lowpass filter (1.3 kHz)
notchFreq = 60;    % Frequency to remove (60 Hz)
order = 4;         % Filter order for lowpass filter

% Design notch filter to remove 60 Hz
notchQ = 35;       % Quality factor - higher Q = narrower notch
wo = notchFreq/(Fs/2);  % Normalized frequency
bw = wo/notchQ;         % Bandwidth
[b_notch, a_notch] = iirnotch(wo, bw);

% Design and apply lowpass filter
Wn_low = cutoffFreq/(Fs/2);
[b_low, a_low] = butter(order, Wn_low, 'low');

processedTraces=NaN(size(tracesRaw));
for iTrace = 1:1
    tempData = tracesRaw(:,iTrace);
    lowpassData = filtfilt(b_low, a_low, tempData);  
    % Apply notch filter to the already lowpass filtered data
    currentData = filtfilt(b_notch, a_notch, lowpassData);
    currentFiltered = medfilt1(currentData, 5);
    smoothWindow = ceil(1 / samplingInterval);
    currentSmoothed = smooth(currentFiltered, smoothWindow);
    currentCentered = currentSmoothed-mean(currentSmoothed);
    processedTraces(:,iTrace)=currentCentered;
end

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
samplesIn200ms = 200 * samplesPerMs;
x_start = 19000;  % Position near the end of the visible range
y_position = -35;
plot([x_start, x_start + samplesIn200ms], [y_position, y_position], 'k', 'LineWidth', 2)
text(x_start + samplesIn200ms/2, y_position - 3, '200 ms', 'HorizontalAlignment', 'center')

% Vertical scale bar (10 pA)
x_position = 19800;  % Position near the end of the visible range
y_start = -35;
plot([x_position, x_position], [y_start, y_start + 10], 'k', 'LineWidth', 2)
text(x_position + 200, y_start + 5, '10 pA', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')


DARTexFile = 11
DARTbaselineSweep=1;

config = file_config{DARTexFile};
filename = config{1};
baseline_sweeps = config{2};
drug_sweeps = config{3};
group = config{4};
cell_id = char(filename(1:8));


load(filename); %load the current results file
abf_path = fullfile('\home\celine\Data\patch_data',[cell_id, '.abf']);
[data, si, h] = abfload(abf_path);
samplingInterval = si/1000; % Convert to ms
samplesPerMs = 1/samplingInterval;
timeVector = (0:size(data,1)-1)' * samplingInterval;
Fs = samplesPerMs*1000;

baselineEventsTimes_DART = [];
baselineEventsTimes_DART = [baselineEventsTimes_DART, allResults(DARTbaselineSweep).epscs.time];
eventTimes_DART = baselineEventsTimes_DART*samplesPerMs;

%collect the two example traces
tracesRaw = [data(:,1,DARTbaselineSweep)];

if h.recChUnits{1,1}=='nA'
    tracesRaw = tracesRaw*1000;
end

%process the example traces 

cutoffFreq = 1300; % Cutoff frequency for lowpass filter (1.3 kHz)
notchFreq = 60;    % Frequency to remove (60 Hz)
order = 4;         % Filter order for lowpass filter

% Design notch filter to remove 60 Hz
notchQ = 35;       % Quality factor - higher Q = narrower notch
wo = notchFreq/(Fs/2);  % Normalized frequency
bw = wo/notchQ;         % Bandwidth
[b_notch, a_notch] = iirnotch(wo, bw);

% Design and apply lowpass filter
Wn_low = cutoffFreq/(Fs/2);
[b_low, a_low] = butter(order, Wn_low, 'low');

processedTraces=NaN(size(tracesRaw));
for iTrace = 1:1
    tempData = tracesRaw(:,iTrace);
    lowpassData = filtfilt(b_low, a_low, tempData);  
    % Apply notch filter to the already lowpass filtered data
    currentData = filtfilt(b_notch, a_notch, lowpassData);
    currentFiltered = medfilt1(currentData, 5);
    smoothWindow = ceil(1 / samplingInterval);
    currentSmoothed = smooth(currentFiltered, smoothWindow);
    currentCentered = currentSmoothed-mean(currentSmoothed);
    processedTraces(:,iTrace)=currentCentered;
end


subplot(2,1,2)
plot(processedTraces(:,1),'b')
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

print(fullfile(outpath,'ExampleTraces.pdf'),'-dpdf');