%% reduced the total size of NS6, keep only trigger dat
%% and also need to compensate for timestamp difference!
function trigger_data(filename,NS6)
% deal with if it is a cell
if iscell(NS6.Data)
    NS6.Data = cell2mat(NS6.Data);
end 
labjackID=find(NS6.MetaTags.ChannelID==136);
photoID=find(NS6.MetaTags.ChannelID==135);
LedID=find(NS6.MetaTags.ChannelID==137);

Trigger.data(1,:)=double(NS6.Data(photoID,:));% 1 is the photodiode trigger
Trigger.data(2,:)=double(NS6.Data(labjackID,:));% 2 is the labjack trigger
Trigger.data(3,:)=double(NS6.Data(LedID,:));% 2 is the labjack trigger
Trigger.photoID = 1;
Trigger.labjackID = 2;
Trigger.LedID = 3;
Trigger.SamplingFreq =  NS6.MetaTags.SamplingFreq ; 
Trigger.Timestamp = NS6.MetaTags.Timestamp;
save(filename,'Trigger','-v7.3')
end

