%% get the trigger from the raw data
Fs=double(Trigger.SamplingFreq);
dt = 1/Fs;

fpass  =    [80 400];  % pass freq
fstop  =    [70 500];  % stop freq outside pass
Rpass  =    0.5;       % Attenuation(dB) for passband
Astop  =    30;        % Attenuation(dB) for stopband
n      =    cheb2ord(fpass/Fs*2,fstop/Fs*2,Rpass,Astop);   % order of chey filter

[z,p,k] =   cheby2(n,Astop,fstop/Fs*2);   % zeros, poles, and gain
[s,g]   =   zp2sos(z,p,k);                  % create second order section
Hd      =   dfilt.df2sos(s,g);                 % dfilt object
photo_filter  =    filtfilthd(Hd,Trigger.data(Trigger.photoID,:));    % apply filter




photo_Subtract  = Trigger.data(Trigger.photoID,:)-mean(Trigger.data(Trigger.photoID,:));
photo_N         =   photo_Subtract./max(photo_Subtract);


photo_filter_N  =   photo_filter./max(photo_filter);
photo_threshold = photo_filter_N>=0.5;

%photo_threshold = photo_N  >=0.4;% deal with some odd conditions, date: 161208 ; 161215 ; 161214

labjack_N  =  Trigger.data(Trigger.labjackID,:)./max(Trigger.data(Trigger.labjackID,:));
labjack_threshold=labjack_N >=0.5;


%photo_threshold = photo_N  >=0.4;% deal with some odd conditions, date: 161208 ; 161215 ; 161214

labjack_N  =  Trigger.data(Trigger.labjackID,:)./max(Trigger.data(Trigger.labjackID,:));
labjack_threshold=labjack_N >=0.5; % threshold labjack signal from labjack noise? 

if input.doBlock2 ==1
if  isfield(input, 'tStimOffTimes')   
input.tStimOffTimes = cellfun(@(x) x(2:end),input.tStimOffTimes, 'UniformOutput', false);
end
% also get the led trigger 
Led_N  =  Trigger.data(Trigger.LedID,:)./max(Trigger.data(Trigger.LedID,:));
Led_threshold=Led_N >=0.3;
else
Led_threshold = [];    
end


%% clear up the photodiode trigger and sort the trial start
index_lab=find(labjack_threshold==1);
photo_clean=zeros(1,length(photo_threshold));
photo_clean(index_lab(1):(index_lab(end)+int64(Fs)))=photo_threshold(index_lab(1):(index_lab(end)+int64(Fs)));


index_photo=find(photo_clean==1);
diff_index=diff(index_photo);
a=find(diff_index>=(input.itiTimeMs./1000)*Fs)+1; 

TrialNum=numel(a)+1;
%%------manually put trial number if crashes--------
%  TrialNum = 66;
if TrialNum <Trialnum
   manual=1;
else
    manual=0;
end
TrialNum=min([TrialNum,Trialnum]);
% TrialNum = 306;
%% manually modify TrialNum if code breaks!
% need to modify the corresponding Index !
% Index_part = {}; % store as Index{Led,1}{Orien,1} orientation level start from smaller to bigger..without led to with led
% for i_led = 1:2
%     for orien_N = 1:numel(Orienlevel)
%         Index_part{i_led,1}{orien_N,1} = intersect(Index{i_led,1}{orien_N,1},1:TrialNum);
%     end
% end

trialstart=zeros(1,TrialNum);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(TrialNum-1))); % trial start timestamp align with photodiode signal
%index is on sampling rate space so divide by 30000 

%% get the timestamp of each trial onset
index_photo=find(photo_threshold==1);
diff_index=diff(index_photo);
a=find(diff_index>=(0.5)*Fs)+1; 
Trials = length(a);
trialstart=zeros(1,Trials);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(Trials-1))); % trial start timestamp align with photodiode signal

trialstart_s = trialstart./Fs;

%% get the led start and stop timestamps 

if input.doBlock2 ==1
idx_LED = find(Led_threshold ==1);
diff_LED = diff(idx_LED);

aa=find(diff_LED>=(input.itiTimeMs./1000)*Fs)+1; 

LED_time = NaN(length(aa)+1,2);
LED_time(1,1) = idx_LED(1);
LED_time(2:end,1) = idx_LED(aa); 
LED_time(1:(end-1),2) = idx_LED(aa-1);
LED_time(end,2) = idx_LED(end);
else
    LED_time = [];
end 
    


%% get the all the index about retinotopic position and orientation

%%get orientation degree and LED conditions
Contrast=cell2mat(input.tGratingContrast);
Contrast=Contrast(1:Trials);
unContrast=unique(Contrast);

Block2 =cell2mat(input.tBlock2TrialNumber);

isBlock2 = find(Block2 == 1) 

isControl = find(Block2 == 0)

ncon = length(unContrast);
% get orientation degree index
index_con_control = cell(1, ncon);
for i = 1: length(unContrast)
    index_con_control{i} = intersect(find(Contrast == unContrast(i)), find(Block2 == 0));
end

index_con_LED = cell(1, ncon);
for i = 1: length(unContrast)
    index_con_LED{i} = intersect(find(Contrast == unContrast(i)), find(Block2 == 1));
end

%% get the datacollapsed by contrast
Data_Collapsed_control = {}; % first position, then electrode

for unit_num = 1:length(spikes.spiketime_S)
for i_con = 1:ncon
           SpikeTimeS = spikes.spiketime_S{unit_num,1};
           VstartTimes = trialstart_s(index_con_control{1,i_con});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data_Collapsed_control{i_con,1}{unit_num,1}{j,1} = spike_aligned;
               else
               Data_Collapsed_control{i_con,1}{unit_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
end
    
Data_Collapsed_LED = {}; % first position, then electrode

for unit_num = 1:length(spikes.spiketime_S)
for i_con = 1:ncon
           SpikeTimeS = spikes.spiketime_S{unit_num,1};
           VstartTimes = trialstart_s(index_con_LED{1,i_con});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data_Collapsed_LED{i_con,1}{unit_num,1}{j,1} = spike_aligned;
               else
               Data_Collapsed_LED{i_con,1}{unit_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
    end

%% get the raw data
Data_control = {}; % data is structured with first separate the position, then orientation degree, then separate out the electrode
for unit_num = 1:length(spikes.spiketime_S)
for i_con = 1: ncon
           trialstart_index =[];
           VstartTimes =[];
           SpikeTimeS = [];
           SpikeTimeS = spikes.spiketime_S{unit_num,1};
           VstartTimes = trialstart_s(index_con_control{1,i_con});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data_control{i_con,1}{unit_num,1}{j,1} = spike_aligned;
               else
               Data_control{i_con,1}{unit_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
end
    
Data_LED = {}; % data is structured with first separate the position, then orientation degree, then separate out the electrode
for i_con = 1: ncon
       for unit_num = 1:length(spikes.spiketime_S)
           trialstart_index =[];
           VstartTimes =[];
           SpikeTimeS = [];
           SpikeTimeS = spikes.spiketime_S{unit_num,1};
           VstartTimes = trialstart_s(index_con_LED{1,i_con});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data_LED{i_con,1}{unit_num,1}{j,1} = spike_aligned;
               else
               Data_LED{i_con,1}{unit_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
    end


%% get the psth data for collapsed 
% difine the bin edges
binnum = 50;
edges = linspace(-0.2, 1.8, binnum);
binwidth = (edges(end)-edges(1)) ./ numel(edges);

for i_con = 1: ncon
       for unit_num = 1:length(spikes.spiketime_S)
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data_Collapsed_control{i_con,1}{unit_num,1}, repmat({edges}, numel(Data_Collapsed_control{i_con,1}{unit_num,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth_Collapsed_control{i_con,1}{unit_num,1 } = counts_per_bin;
       tuning_Collapsed_control{i_con,1}{unit_num,1} = sum(counts_stimOn,2); %sum of each column (unit)
       end
       psth_Collapsed_avg_control{i_con,1} = cellfun(@(x) mean(x,1), psth_Collapsed_control{i_con,1}, 'uniformoutput', false);
       psth_Collapsed_avg_control{i_con,1} = cell2mat(psth_Collapsed_avg_control{i_con,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_Collapsed_sem_control{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth_Collapsed_control{i_con,1}, 'uniformoutput', false);
       psth_Collapsed_sem_control{i_con,1} = cell2mat(psth_Collapsed_sem_control{i_con,1});
       tuning_Collapsed_avg_control{i_con,1} = cellfun(@(x) mean(x,1), tuning_Collapsed_control{i_con,1}, 'uniformoutput', false);
       tuning_Collapsed_avg_control{i_con,1} = cell2mat(tuning_Collapsed_avg_control{i_con,1});
       tuning_Collapsed_sem_control{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning_Collapsed_control{i_con,1}, 'uniformoutput', false);
       tuning_Collapsed_sem_control{i_con,1} = cell2mat(tuning_Collapsed_sem_control{i_con,1});
end

for i_con = 1: ncon
       for unit_num = 1:length(spikes.spiketime_S)
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data_Collapsed_LED{i_con,1}{unit_num,1}, repmat({edges}, numel(Data_Collapsed_LED{i_con,1}{unit_num,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth_Collapsed_LED{i_con,1}{unit_num,1 } = counts_per_bin;
       tuning_Collapsed_LED{i_con,1}{unit_num,1} = sum(counts_stimOn,2);
       end
       psth_Collapsed_avg_LED{i_con,1} = cellfun(@(x) mean(x,1), psth_Collapsed_LED{i_con,1}, 'uniformoutput', false);
       psth_Collapsed_avg_LED{i_con,1} = cell2mat(psth_Collapsed_avg_LED{i_con,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_Collapsed_sem_LED{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth_Collapsed_LED{i_con,1}, 'uniformoutput', false);
       psth_Collapsed_sem_LED{i_con,1} = cell2mat(psth_Collapsed_sem_LED{i_con,1});
       tuning_Collapsed_avg_LED{i_con,1} = cellfun(@(x) mean(x,1), tuning_Collapsed_LED{i_con,1}, 'uniformoutput', false);
       tuning_Collapsed_avg_LED{i_con,1} = cell2mat(tuning_Collapsed_avg_LED{i_con,1});
       tuning_Collapsed_sem_LED{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning_Collapsed_LED{i_con,1}, 'uniformoutput', false);
       tuning_Collapsed_sem_LED{i_con,1} = cell2mat(tuning_Collapsed_sem_LED{i_con,1});
end

% get the psth data


    for i_con = 1: ncon
       for unit_num = 1:length(spikes.spiketime_S)   
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data_control{i_con,1}{unit_num,1}, repmat({edges}, numel(Data_control{i_con,1}{unit_num,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth_control{i_con,1}{unit_num,1 } = counts_per_bin;
       tuning_control{i_con,1}{unit_num,1} = sum(counts_stimOn,2);
       end
       psth_avg_control{i_con,1}= cellfun(@(x) mean(x,1), psth_control{i_con,1}, 'uniformoutput', false);
       psth_avg_control{i_con,1} = cell2mat(psth_avg_control{i_con,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_sem_control{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth_control{i_con,1}, 'uniformoutput', false);
       psth_sem_control{i_con,1} = cell2mat(psth_sem_control{i_con,1});
       tuning_avg_control{i_con,1}= cellfun(@(x) mean(x,1), tuning_control{i_con,1}, 'uniformoutput', false);
       tuning_avg_control{i_con,1} = cell2mat(tuning_avg_control{i_con,1});
       tuning_sem_control{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning_control{i_con,1}, 'uniformoutput', false);
       tuning_sem_control{i_con,1} = cell2mat(tuning_sem_control{i_con,1});
    end
    
    
    for i_con = 1: ncon
       for unit_num = 1:length(spikes.spiketime_S)   
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data_LED{i_con,1}{unit_num,1}, repmat({edges}, numel(Data_LED{i_con,1}{unit_num,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth_LED{i_con,1}{unit_num,1 } = counts_per_bin;
       tuning_LED{i_con,1}{unit_num,1} = sum(counts_stimOn,2);
       end
       psth_avg_LED{i_con,1}= cellfun(@(x) mean(x,1), psth_LED{i_con,1}, 'uniformoutput', false);
       psth_avg_LED{i_con,1} = cell2mat(psth_avg_LED{i_con,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_sem_LED{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth_LED{i_con,1}, 'uniformoutput', false);
       psth_sem_LED{i_con,1} = cell2mat(psth_sem_LED{i_con,1});
       tuning_avg_LED{i_con,1}= cellfun(@(x) mean(x,1), tuning_LED{i_con,1}, 'uniformoutput', false);
       tuning_avg_LED{i_con,1} = cell2mat(tuning_avg_LED{i_con,1});
       tuning_sem_LED{i_con,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning_LED{i_con,1}, 'uniformoutput', false);
       tuning_sem_LED{i_con,1} = cell2mat(tuning_sem_LED{i_con,1});
    end


% get the collapsed retinotopic data
 tuningCollapsedplot={};
 tuningCollapsedplot_sem={};
 for i_con = 1:ncon
     for unit = 1:length(spikes.spiketime_S)   
        
             tuningCollapsedplot{i_con,1}(unit,1) = tuning_Collapsed_avg{i_con,1}(unit);
             tuningCollapsedplot_sem{i_con,1}(unit,1)= tuning_Collapsed_sem{i_con,1}(unit);
            
         end
     end

% need to retinotopic and orientation tuning curve data
 tuningplot={};
 tuningplot_sem={};
 for i_con = 1:ncon
     for unit = 1:length(spikes.spiketime_S)   
             tuningplot{i_con,1}{unit,1} = tuning_avg{i_con,1}(unit);
             tuningplot_sem{i_con,1}{unit,1}= tuning_sem{i_con,1}(unit);
            
         end
     end

%%
 singleshank=1;
if singleshank==1

 shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
 shanks = [1,2];
 electrodes=16;
else
  shankchannel=[7 14 4 6 8 5 3 13;11 1 12 10 15 2 16 9;32 22 30 21 31 18 19 17;27 26 24 28 25 20 23 29];
  shanks = [1:4];
  electrodes=8;
end
%% plot retinotopic figure for collapsed data
% plot retinotopic for each electrode site, aligned by eah shank
% row: high elevation->low ; column: media azimuth to lateral
retino_organized = {};
for i_elec = 1:32
    for i_pos = 1: npos
        retino_organized{i_elec,1}(i_pos,1)=tuningCollapsedplot{i_pos,1}(i_elec);
    end
end

retino_collapsed={};
best_response = [];
for i_elec=1:32
    best_response(i_elec,1)=max(retino_organized{i_elec});
    retino_collapsed{i_elec,1}(1,1:3)=retino_organized{i_elec}(4:6,1);
    retino_collapsed{i_elec,1}(2,1:3)=retino_organized{i_elec}(1:3,1);
end
%% plot retinotopic for each electrode site, aligned by eah shank
azimuth={num2str(pos(1,2));num2str(pos(2,2));num2str(pos(3,2))};
elevation = {num2str(pos(4,1));num2str(pos(1,1))};
for shank=shanks
    figure
    for j=1:electrodes
    elec=shankchannel(shank,j);
    if best_response(elec,1)
    subplot(1,electrodes,j)
    col=[0 1];
    imagesc(retino_collapsed{elec}./best_response(elec,1), col) % normalized /best_response(elec,1)
%     colorbar
    set(gca,'XTick',1:1:3,'YTick',1:1:2,'XGrid','off','XTickLabel',azimuth,'YTickLabel',elevation);
    title(['S:',num2str(shank),'  E:',num2str(elec)]);
    end
    end
   
end
%% find LED time offset

trialstart_block2 = trialstart(find(cell2mat(input.tBlock2TrialNumber)))
LED_ontime = trialstart_block2 - LED_time(:,1)
LED_ontime_S = LED_ontime./30000
avg_LEDontime_S = mean(ans)


%% plot the traces for each electrode
figure; hold on
plot(edges,psth_avg_control{1,1},'Color', [.2 .2 .2])
ylim([0 2])
xlim([-.2 1.8])
axis('square')
figure; hold on
plot(edges,psth_avg_control{2,1},'Color', [.4 .4 .4])
ylim([0 2])
xlim([-.2 1.8])
axis('square')
figure; 
plot(edges,psth_avg_control{3,1},'Color', [.6 .6 .6])
ylim([0 2])
xlim([-.2 1.8])
axis('square')

figure; hold on
plot(edges,psth_Collapsed_avg_LED{1,1},'Color', [.2 .2 .2])
ylim([0 10])
xlim([-.2 1.8])
axis('square')
figure; hold on
plot(edges,psth_Collapsed_avg_LED{2,1},'Color', [.4 .4 .4])
ylim([0 5])
xlim([-.2 1.8])
axis('square')
figure; hold on 
plot(edges,psth_Collapsed_avg_LED{3,1},'Color', [.6 .6 .6])
plot(-.0304,1, '.' ) 
plot(.4696,1, '.' )
ylim([0 12])
xlim([-.2 .5])


ylim([0 12])
xlim([-.2 .5])
axis('square')

for i = 1:length(spikes.spiketime_S)  
figure; hold on
title('Grating Contrast = .8') 
h = refline(1000000, 0) 
h.Color = 'k' 
h.LineStyle = '--'
plot(edges,psth_avg_control{1,1}(i,:),'Color', [.2 .2 .2])
fill([edges, fliplr(edges)], [(psth_avg_control{1,1}(i,:) - psth_sem_control{1,1}(i,:)), fliplr(psth_avg_control{1,1}(i,:) + psth_sem_control{1,1}(i,:))], [.2 .2 .2], 'EdgeColor', [.2 .2 .2])
alpha .2 
ylim([0 5])
xlim([-.2 .5])
axis('square')
hold on
plot(edges,psth_avg_LED{1,1}(i,:),'Color', [.2 .2 1])
fill([edges, fliplr(edges)], [(psth_avg_LED{1,1}(i,:) - psth_sem_LED{1,1}(i,:)), fliplr(psth_avg_LED{1,1}(i,:) + psth_sem_LED{1,1}(i,:))], [.2 .2 1], 'EdgeColor', [.2 .2 1])
alpha .2 
ylim([0 5])
xlim([-.2 .5])
axis('square')
figure; hold on
title('Grating Contrast = .5') 
h = refline(1000000, 0) 
h.Color = 'k' 
h.LineStyle = '--'
plot(edges,psth_avg_control{2,1}(i,:),'Color', [.4 .4 .4])
fill([edges, fliplr(edges)], [(psth_sem_control{2,1}), fliplr(psth_sem_control{2,1})],[.4 .4 .4], 'EdgeColor', 'none')
fill([edges, fliplr(edges)], [(psth_avg_control{2,1}(i,:) - psth_sem_control{2,1}(i,:)), fliplr(psth_avg_control{2,1}(i,:) + psth_sem_control{2,1}(i,:))], [.4 .4 .4], 'EdgeColor', [.4 .4 .4])
alpha .2
ylim([0 5])
xlim([-.2 .5])
axis('square')
hold on
plot(edges,psth_avg_LED{2,1}(i,:),'Color', [.4 .4 1])
fill([edges, fliplr(edges)], [(psth_sem_LED{2,1}), fliplr(psth_sem_LED{2,1})],[.4 .4 .4], 'EdgeColor', 'none')
fill([edges, fliplr(edges)], [(psth_avg_LED{2,1}(i,:) - psth_sem_LED{2,1}(i,:)), fliplr(psth_avg_LED{2,1}(i,:) + psth_sem_LED{2,1}(i,:))], [.4 .4 1], 'EdgeColor', [.4 .4 1])
alpha .2
ylim([0 5])
xlim([-.2 .5])
axis('square')
figure; hold on
title('Grating Contrast = .2') 
h = refline(1000000, 0) 
h.Color = 'k' 
h.LineStyle = '--'
plot(edges,psth_avg_control{3,1}(i,:),'Color', [.6 .6 .6])
fill([edges, fliplr(edges)], [(psth_avg_control{3,1}(i,:) - psth_sem_control{3,1}(i,:)), fliplr(psth_avg_control{3,1}(i,:) + psth_sem_control{3,1}(i,:))], [.6 .6 .6], 'EdgeColor', [.6 .6 .6])
alpha .2
ylim([0 5])
xlim([-.2 .5])
axis('square')
hold on
plot(edges,psth_avg_LED{3,1}(i,:),'Color', [.6 .6 1])
fill([edges, fliplr(edges)], [(psth_avg_LED{3,1}(i,:) - psth_sem_LED{3,1}(i,:)), fliplr(psth_avg_LED{3,1}(i,:) + psth_sem_LED{3,1}(i,:))], [.6 .6 1], 'EdgeColor', [.6 .6 1])
alpha .2
ylim([0 10])
xlim([-.2 .5])
axis('square')
end



% (# of spikes in time window for all trials)/ (# of trials * time window
% duration) 


