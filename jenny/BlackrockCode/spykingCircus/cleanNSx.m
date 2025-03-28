clear


date = input('Date','s');
data_number = input('Data number','s');

doAll = input('merge all ns6 files?','s');
if strcmp(doAll,'y')
    tempPath = ['Z:\home\Naomi\From_Jenny\Practice_Analysis',date,'\'];
    % tempPath = ['Z:\All_staff\home\jen\Raw Data\Multiunit\',date,'\'];
    tempPath = getFolder;
    ns6_files = dir(fullfile(tempPath,'*.ns6'));
    NS62 = cellfun(@(x) openNSx([tempPath,'\', x]),{ns6_files.name},'un',0);
    count = numel(NS62);
else
    
NS62{1} = openNSx;

waitFlag = input('merge?','s');
count = 1;
while waitFlag == 'y'
    NS62{count+1} = openNSx;
    waitFlag = input('merge?','s');
    count = count+1;
end
end
%% check if the channel ID is correct for poly shanks 

camID=34;
photoID=33;
LedID=35;
data = [];

for i = 1:count
    if iscell(NS62{i}.Data)
        NS62{i}.Data = cell2mat(NS62{i}.Data);
    end
    data = [data NS62{i}.Data(:,:)];
        
    Trigger.data{i}(1,:)= double(NS62{i}.Data(photoID,:));% 1 is the photodiode trigger
    Trigger.data{i}(2,:)= double(NS62{i}.Data(camID,:));% 2 is the camera frames trigger
    Trigger.data{i}(3,:)= double(NS62{i}.Data(LedID,:));% 3 is the LED trigger
end

if count == 1
    Trigger.data = Trigger.data{1};
end
Trigger.photoID = 1;
Trigger.camID = 2;
Trigger.LedID = 3;
Trigger.SamplingFreq =  NS62{1}.MetaTags.SamplingFreq ; 
Trigger.Timestamp = NS62{1}.MetaTags.Timestamp;

% data = data(:,1:105000000);
data = data(1:end-3,:);
data = data(:);

% savepath = ['Z:\All_staff\home\jen\Analysis\Spyking_Circus\',date,'\',data_number,'\data\'];
savepath = ['Z:\home\Naomi\From_Jenny\Practice_Analysis\SPCircus\',date,'\',data_number,'\data\'];
if ~exist(savepath)
    % mkdir('Z:\All_staff\home\jen\Analysis\Spyking_Circus\',[date,'\',data_number,'\data\']);
    mkdir('Z:\home\Naomi\From_Jenny\Practice_Analysis\SPCircus\',[date,'\',data_number,'\data\']);
end
save_names = strsplit(savepath,'\');
save([savepath,save_names{7},'_',save_names{8},'_triggers'],'Trigger','-v7.3')
newFilename = [savepath,'data.dat'];

% Opening the output file for saving
FIDw = fopen(newFilename, 'w+', 'ieee-le');

% Writing data into file
disp('Writing the converted data into the new .dat file...');
fwrite(FIDw, data, 'int16');
fclose(FIDw);