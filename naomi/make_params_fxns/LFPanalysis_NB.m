% function layerBound = LFPanalysis(filename,trigFile)

filename=[rawData_path,layersFile,'.ns4']
trigFile=params.info.triggersFile
%
NS4 = openNSx('read',filename);
%%
% select etype

param.filename = NS4.MetaTags.Filename;
param.eLabels = vertcat(NS4.ElectrodesInfo.Label);
tmp_filename = strsplit(param.filename,'.');
temp_date = strsplit(NS4.MetaTags.FilePath,'\');
%%

temp_table = readtable('Z:\home\Naomi\Electrophysiology\New Analysis\probeLogNB.xlsx','Sheet',1);
% temp_table = readtable('Z:\All_Staff\home\jen\Notebook files\probeLog.xlsx','Sheet',1);
date_idx = contains(temp_table.Dates,temp_date{6});
if ~isempty(date_idx)
sheet_idx = temp_table.comboID{date_idx};
else
    disp([temp_table.comboID temp_table.Probe temp_table.HS]);
    waitflag = input('did not find log date. select recording config (A-E):','s');
end
%%
coordLog = readtable('Z:\home\Naomi\Electrophysiology\New Analysis\probeLogNB.xlsx','Sheet',sheet_idx);
% coordLog = readtable('Z:\All_Staff\home\jen\Notebook files\probeLog.xlsx','Sheet',sheet_idx);
shankchannel = coordLog.chanID;
elec_pos = coordLog.yCoord';
elec_col = coordLog.xCoord;

param.eSizeCol = numel(unique(elec_col));
param.eSizeRow = numel(shankchannel)/param.eSizeCol;

shankchannel = reshape(shankchannel,param.eSizeRow,param.eSizeCol)';
elec_pos = reshape(elec_pos,param.eSizeRow,param.eSizeCol)';
%%

temp_filename = [temp_date{6},'_',tmp_filename{1},'_triggers.mat'];
if isempty(trigFile)
% addpath(genpath(['Z:\All_staff\home\jen\Analysis\Spyking_Circus\',temp_date{end},'\']));
addpath(genpath(['Z:\home\Naomi\Electrophysiology\BlackRock\',temp_date{6},'\']));

load(temp_filename);
else
    if exist(trigFile)
        load(trigFile);
        
       tmpName = strsplit(trigFile,'_');
        saveName = tmpName{end-1};
        
        try 
            temp = strsplit(saveName,'-');
        catch
            temp = {saveName};
        end
        file_idx = find(cell2mat(cellfun(@(x) strcmp(x(end),filename(end-4)),temp,'un',0)));
        try
        Trigger.data = Trigger.data{file_idx};
        catch
        Trigger.data = Trigger.data;
        end
        
    else
        disp('cannot find triggers. exiting');
    end
end
    
triggers = Trigger.data(Trigger.photoID,:);
param.sampleInterval = 1/Trigger.SamplingFreq;

% clear Trigger 
%%
if iscell(NS4.Data)
    LFP.rawTrace = cell2mat(NS4.Data);
else
    LFP.rawTrace = NS4.Data;
end

param.sampleRate = NS4.MetaTags.SamplingFreq;
param.threshold = .4;
param.bandPassFreq = [0.5 250];
param.nElectrodes = param.eSizeCol*param.eSizeRow;
param.afterStimS = 0.5;
time = 0:(1/param.sampleRate):param.afterStimS;

%% photodiode stim on

photo_dio = triggers;

Fs = 30000; % this is from the NS6 file
fpass  =    [80 400];  % pass freq
fstop  =    [70 500];  % stop freq outside pass
Rpass  =    0.5;       % Attenuation(dB) for passband
Astop  =    30;        % Attenuation(dB) for stopband
n      =    cheb2ord(fpass/Fs*2,fstop/Fs*2,Rpass,Astop);   % order of chey filter

[z,p,k] =   cheby2(n,Astop,fstop/Fs*2);   % zeros, poles, and gain
[s,g]   =   zp2sos(z,p,k);                  % create second order section
Hd      =   dfilt.df2sos(s,g);                 % dfilt object


photo_dio =    filtfilthd(Hd,photo_dio);    % apply filter

photo_dio = bsxfun(@rdivide,photo_dio,(max(photo_dio,[],2)));
is_dio_on = photo_dio > param.threshold * (max(photo_dio(:)));

temp_stimOn = find(is_dio_on);
temp_stimOn2 = diff([-10000 temp_stimOn])>Fs;
stimOn = temp_stimOn(temp_stimOn2);

clear waitflag
waitflag = input('1 or 2 stim');

if waitflag == 2 % case two stim
    stimOn = floor(stimOn/3);
    try
    param.stimOn = [stimOn(1:2:end)' stimOn(2:2:end)'];
    catch
            param.stimOn = [stimOn(1:2:end-1)' stimOn(2:2:end)'];
    end
    dio_idx = 2;
elseif waitflag == 1  % one stim 
    stimOn = floor(stimOn/3);
    param.stimOn = stimOn';
    dio_idx = 1;
end


%% filter for LFP and get stim on
[b,a] = butter(1,param.bandPassFreq/(param.sampleRate));

for column = 1:param.eSizeCol
    for row = 1:param.eSizeRow
        elec = shankchannel(column,row);
% elec = row;
        LFP_smooth = filtfilt(b,a,double(LFP.rawTrace(elec,:)));
        temp_struct = cleanData(LFP_smooth,1,param.stimOn(:,dio_idx),'pclamp',false,'sampleFreq',param.sampleRate,'afterStimS',param.afterStimS,'rebaseline',false);
        LFPbyStim{row,column} = temp_struct.pts;      
    end
end

%% make some corrections for dead channels and bad trials

LFP.mean = cellfun(@(x) nanmean(x,2)',LFPbyStim,'un',0);
for column = 1:param.eSizeCol
   LFP.meanByCol{column} = cell2mat(LFP.mean(:,column));
end

%
for column = 1:param.eSizeCol
figure
if param.eSizeRow < 10
for row = 1:param.eSizeRow
    depth = elec_pos(column,:);
    ax(row) = subplot(param.eSizeRow,1,row);
    hold on;
    plot(time,LFP.meanByCol{column}(row,:))
    title([num2str(depth(row)),'(',num2str(row),')']); 
    ylim([-800 800]);xlim([0 param.afterStimS]);
end
else 
for row = 1:param.eSizeRow
    depth = elec_pos(column,:);
    ax(row) = subplot(param.eSizeRow/4,4,row);
    hold on;
    plot(time,LFP.meanByCol{column}(row,:))
    title([num2str(depth(row)),'(',num2str(row),')']); 
    ylim([-800 800]);xlim([0 param.afterStimS]);
end
end
xlabel('time (s)')
linkaxes(ax);


    bad_channels_col = input('which channels are bad on this col?');
    for badchannel_i = bad_channels_col
        if badchannel_i == 1
            LFP.mean{1,1} = LFP.mean{2,1};
        elseif badchannel_i == param.eSizeRow
            LFP.mean{badchannel_i,1} = LFP.mean{badchannel_i-1,1};
        else
            temp_diff = setdiff(1:param.eSizeRow,bad_channels_col);
            temp_replace = temp_diff-badchannel_i;
            ch1_max = temp_diff(find(temp_replace>0,1,'first'));
            ch1_min = temp_diff(find(temp_replace<0,1,'last'));
            LFP.mean{badchannel_i,1} = mean([LFP.mean{ch1_max,1}; LFP.mean{ch1_min,1}]);
        end
    end
if param.eSizeCol == 1
space_i = 2;
else
    space_i = 1;
end
end
%%
for column = 1:param.eSizeCol
   LFP.meanByCol{column} = cell2mat(LFP.mean(:,column));
end

%%
cmap_depth = brewermap(param.eSizeRow,'RdGy');
figure; hold on;
for row = 1:param.eSizeRow
    plot(1000*time,LFP.meanByCol{1}(row,:),'Color',cmap_depth(row,:))
end
fix_axes(gcf,20,'time','uV');
legend(param.eLabels(shankchannel(1:param.eSizeRow),:))
% xlim([0 500])


%%

for column = 1:param.eSizeCol
    LFP.devDev=-(diff(LFP.meanByCol{column}(1:space_i:end,:),2,1));
    depth = elec_pos(column,1:space_i:end);
    figure
    pcolor(time,-depth(2:end-1),LFP.devDev);
    xlim([0 0.5])
    shading('interp');
    colormap('jet');
    caxis([-80 80]);

end


%%

% click layer depth distinguishers
depthL23 = -ginput(2);
depthL4 = -ginput(1);

% find contact site closest to depth clicked
% tempMin = abs(elec_pos-floor(depthL23(1,2)));
% [~,idx1] = min(tempMin(:));
% [pos1,pos2] = ind2sub(size(tempMin),idx1);
% layerBound.layer23_depth(1) = elec_pos(pos1,pos2);
layerBound.layer23_depth(1) = depthL23(1,2);

% tempMin = abs(elec_pos-floor(depthL23(2,2)));
% [~,idx2] = min(tempMin(:));
% [pos1,pos2] = ind2sub(size(tempMin),idx2);
% layerBound.layer23_depth(2) = elec_pos(pos1,pos2);
layerBound.layer23_depth(2) = depthL23(2,2);

layerBound.layer4_depth(1) = layerBound.layer23_depth(2);

% tempMin = abs(elec_pos-floor(depthL4(1,2)));
% [~,idx3] = min(tempMin(:));
% [pos1,pos2] = ind2sub(size(tempMin),idx3);
% layerBound.layer4_depth(2) = elec_pos(pos1,pos2);
layerBound.layer4_depth(2) = depthL4(1,2);

%%


temp_filename = strsplit(temp_filename,'_');
% if isempty(trigFile)
    saveName = temp_filename{2};
% save(['Z:\All_staff\home\jen\Analysis\Spyking_Circus\',temp_date{end},'\',...
    save(['Z:\home\Naomi\Electrophysiology\Kilosort_Analysis\',temp_date{6},'\',temp_filename{1},'_',saveName,'_layerBound'],'layerBound')
% else
%     % save(['Z:\All_staff\home\jen\Analysis\Spyking_Circus\',temp_date{end},'\',...
%      save(['Z:\home\Naomi\Electrophysiology\BlackRock\',temp_date{6},'\',temp_filename{1},'_',saveName,'_layerBound'],'layerBound')
% end