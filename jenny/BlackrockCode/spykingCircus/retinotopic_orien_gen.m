

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


%% get the timestamp of each trial onset
index_photo=find(photo_threshold==1);
diff_index=diff(index_photo);
a=find(diff_index>=(0.5)*Fs)+1; 
Trials = length(a);
trialstart=zeros(1,Trials);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(Trials-1))); % trial start timestamp align with photodiode signal

trialstart_s = trialstart./Fs;

%% get the all the index about retinotopic position and orientation

%%get orientation degree and LED conditions
Orientation=cell2mat(input.tGratingDirectionDeg);
Orientation=Orientation(1:Trials);
Degree=unique(Orientation);
%get retinotopic position
Azimuth=cell2mat(input.tGratingAzimuthDeg); 
Azimuth = Azimuth(1:Trials);
Elevation = cellfun(@(x) double(x),input.tGratingElevationDeg);

Elevation = Elevation(1:Trials);


A_pos = unique(Azimuth);
E_pos = unique(Elevation);

ndeg = length(Degree);
nA = length(A_pos);
nE = length(E_pos);
npos = nA*nE;
% get orientation degree index
index_deg = cell(1, ndeg);
for i = 1: length(Degree)
    index_deg{i} = find(Orientation == Degree(i));
end

index_pos = cell(npos,1);
pos = [];% for saving: first column is elevation, second column is azimuth
j=1;
for i_E = 1:nE
    for i_A = 1:nA
    index_pos{j} = intersect(find(Azimuth==A_pos(i_A)), find(Elevation==E_pos(i_E)));
    pos(j,1)=E_pos(i_E); % first column elevation
    pos(j,2)=A_pos(i_A); % second column azimuth
    j=j+1;
    end
end

% 6 position x 4 orientation conditions %
index_all = {}; % in this cell, start with 6 position in sequence, then inside with each orientation

for i_pos = 1:npos
    for i_deg = 1: ndeg
        index_all{i_pos,1}{i_deg,1} = intersect(index_pos{i_pos,1}, index_deg{i_deg});
    end
end
%% get the datacollapsed by orientation
Data_Collapsed = {}; % first position, then electrode
for i_pos = 1: npos
       for Elec_num = 1:32
           
           VstartTimes =[];
           SpikeTimeS = [];
           Elec = num2str(Elec_num);
           Unit = 0; % only for unsorted recording data, *****if sorted, code needs modification************
           [allTimestamps allSnippets allIndices] = findEventTimes(NEV, str2num(Elec), Unit);
           SpikeTimeS = allTimestamps/Fs;
          
           VstartTimes = trialstart_s(index_pos{i_pos,1});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data_Collapsed{i_pos,1}{Elec_num,1}{j,1} = spike_aligned;
               else
               Data_Collapsed{i_pos,1}{Elec_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
    end

%% get the raw data
Data = {}; % data is structured with first separate the position, then orientation degree, then separate out the electrode
for i_pos = 1: npos
    for i_deg = 1: ndeg  
       for Elec_num = 1:32
           trialstart_index =[];
           VstartTimes =[];
           SpikeTimeS = [];
           Elec = num2str(Elec_num);
           Unit = 0; % only for unsorted recording data, *****if sorted, code needs modification************
           [allTimestamps allSnippets allIndices] = findEventTimes(NEV, str2num(Elec), Unit);
           SpikeTimeS = allTimestamps/Fs;
          
           VstartTimes = trialstart_s(index_all{i_pos,1}{i_deg,1});
           trial_num = length(VstartTimes);
           for j = 1: trial_num % for each start timestamp
               spike =[];
               spike_aligned = [];
               spike = SpikeTimeS - VstartTimes(j);% align with visual stim onset take account for delay in screen
               spike_aligned = spike(find(spike>=-0.2 & spike<=1.8));% 1s on 1s off
               if ~isempty(spike_aligned)
               Data{i_pos,1}{i_deg,1}{Elec_num,1}{j,1} = spike_aligned;
               else
               Data{i_pos,1}{i_deg,1}{Elec_num,1}{j,1} = NaN;   
               end
           end
           
       end 
       
    end
end
%% get the psth data for collapsed 
% difine the bin edges
binnum = 50;
edges = linspace(-0.2, 1.8, binnum);
binwidth = (edges(end)-edges(1)) ./ numel(edges);

for i_pos = 1: npos
       for i_etrode = 1:32   
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data_Collapsed{i_pos,1}{i_etrode,1}, repmat({edges}, numel(Data_Collapsed{i_pos,1}{i_etrode,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth_Collapsed{i_pos,1}{i_etrode,1 } = counts_per_bin;
       tuning_Collapsed{i_pos,1}{i_etrode,1} = sum(counts_stimOn,2);
       end
       psth_Collapsed_avg{i_pos,1} = cellfun(@(x) mean(x,1), psth_Collapsed{i_pos,1}, 'uniformoutput', false);
       psth_Collapsed_avg{i_pos,1} = cell2mat(psth_Collapsed_avg{i_pos,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_Collapsed_sem{i_pos,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth_Collapsed{i_pos,1}, 'uniformoutput', false);
       psth_Collapsed_sem{i_pos,1} = cell2mat(psth_Collapsed_sem{i_pos,1});
       tuning_Collapsed_avg{i_pos,1} = cellfun(@(x) mean(x,1), tuning_Collapsed{i_pos,1}, 'uniformoutput', false);
       tuning_Collapsed_avg{i_pos,1} = cell2mat(tuning_Collapsed_avg{i_pos,1});
       tuning_Collapsed_sem{i_pos,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning_Collapsed{i_pos,1}, 'uniformoutput', false);
       tuning_Collapsed_sem{i_pos,1} = cell2mat(tuning_Collapsed_sem{i_pos,1});
end


% get the psth data



for i_pos = 1: npos
    for i_deg = 1: ndeg
       for i_etrode = 1:32   
         % use cellfun to get N spikes per bin
       counts_per_bin=[];
       counts_stimOn =[];
       counts_per_bin = cellfun(@(x,y) histc(x,y), Data{i_pos,1}{i_deg,1}{i_etrode,1}, repmat({edges}, numel(Data{i_pos,1}{i_deg,1}{i_etrode,1}), 1), 'uniformoutput', false); 
       counts_per_bin = cell2mat(counts_per_bin); % convert to a matrix for easy operations later.
       counts_stimOn = counts_per_bin(:,6:((binnum/2)+6));
       psth{i_pos,1}{i_deg,1}{i_etrode,1 } = counts_per_bin;
       tuning{i_pos,1}{i_deg,1}{i_etrode,1} = sum(counts_stimOn,2);
       end
       psth_avg{i_pos,1}{i_deg,1} = cellfun(@(x) mean(x,1), psth{i_pos,1}{i_deg,1}, 'uniformoutput', false);
       psth_avg{i_pos,1}{i_deg,1} = cell2mat(psth_avg{i_pos,1}{i_deg,1}); % each row is a different etrode, and the data are the average across trials (for that etrode)
       psth_sem{i_pos,1}{i_deg,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), psth{i_pos,1}{i_deg,1}, 'uniformoutput', false);
       psth_sem{i_pos,1}{i_deg,1} = cell2mat(psth_sem{i_pos,1}{i_deg,1});
       tuning_avg{i_pos,1}{i_deg,1} = cellfun(@(x) mean(x,1), tuning{i_pos,1}{i_deg,1}, 'uniformoutput', false);
       tuning_avg{i_pos,1}{i_deg,1} = cell2mat(tuning_avg{i_pos,1}{i_deg,1});
       tuning_sem{i_pos,1}{i_deg,1} = cellfun(@(x) std(x,[],1)./sqrt(size(x,1)), tuning{i_pos,1}{i_deg,1}, 'uniformoutput', false);
       tuning_sem{i_pos,1}{i_deg,1} = cell2mat(tuning_sem{i_pos,1}{i_deg,1});
    end
end
% get the collapsed retinotopic data
 tuningCollapsedplot={};
 tuningCollapsedplot_sem={};
 for i_pos = 1:npos
     for elec = 1:32
        
             tuningCollapsedplot{i_pos,1}(elec,1) = tuning_Collapsed_avg{i_pos,1}(elec);
             tuningCollapsedplot_sem{i_pos,1}(elec,1)= tuning_Collapsed_sem{i_pos,1}(elec);
            
         end
     end

% need to retinotopic and orientation tuning curve data
 tuningplot={};
 tuningplot_sem={};
 for i_pos = 1:npos
     for elec = 1:32
         for i_deg = 1:ndeg
             tuningplot{i_pos,1}{elec,1}(i_deg,1) = tuning_avg{i_pos,1}{i_deg,1}(elec);
             tuningplot_sem{i_pos,1}{elec,1}(i_deg,1)= tuning_sem{i_pos,1}{i_deg,1}(elec);
            
         end
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
%% plot the traces for each electrode
i_deg = 1;
for i = 1:32
    figure
    for i_pos = 1:npos
    subplot(2,3,i_pos)
    plot(edges,psth_avg{i_pos,1}{i_deg,1}(i,:))
    xlabel(num2str(pos(i_pos,:)))
    xlim([-0.2 1.8])
    end 
end 
