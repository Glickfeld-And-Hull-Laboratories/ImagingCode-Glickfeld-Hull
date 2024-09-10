function [Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,LED_time]=PhotoISI(Trigger,input)
global Fs stimOnTimeS Offs


target_compensate = 0;


% define the compensate adjustment for photodiode signal based on elevation
if unique(double(cell2mat(input.tGratingHeightDeg)))==30
    compensate =8-unique(double(cell2mat(input.tGratingElevationDeg)))*0.4;
else
    compensate =0;
end

Cycle = double(unique(cell2mat(input.tCyclesOn))); % number of cycles
Cycle_trl = double(cell2mat(input.tCyclesOn)); % how many cycles each trial
Trialnum = input.trialSinceReset;

stimOnTimeS = double(input.stimOnTimeMs)./1000;

% try to align the onset trigger times trial to trial basis?
if isfield(input,'tStimOffTimes')
    Offs = unique(double(cell2mat(input.tStimOffTimes))); % how many off conditions
elseif isfield(input,'tStimOffTimeMs')
    
    Offs = unique(double(cell2mat(input.tStimOffTimeMs)));
else
    Offs = input.stimOffTimeMs;
    
end



% get the all the normalized and filtered trigger
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
photo_threshold = photo_filter_N>=0.5; % threshold photodiode signal from photodiode noise?

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
    



%% get the timestamp of each cycle start
% if block2 is on, then tStimOffTimes has extra 250ms in the
% beginning, get rid of this
if input.doBlock2 ==1
    Led = cell(TrialNum,1);
   
end

photo_cycle_start = cell(TrialNum,1);
photo_cycle_stop = cell(TrialNum,1);
labjack = cell(TrialNum,1);
photo = cell(TrialNum,1);
Time_start = cell(TrialNum,1);
Time_stop = cell(TrialNum,1);
Time_start_new = cell(TrialNum,1);
Time_stop_new = cell(TrialNum,1);

photo_duration = cell(TrialNum,1);
photo_inter = cell(TrialNum,1);
photo_duration_MS = cell(TrialNum,1);
photo_inter_MS =cell(TrialNum,1);
for i=1:1
    Time_start{i,1}=zeros(1,(Cycle_trl(i)+1));
    Time_stop{i,1}=zeros(1,(Cycle_trl(i)+1));
    Time_start_new{i,1} = zeros(1,(Cycle_trl(i)+1));
    Time_stop_new{i,1} = zeros(1,(Cycle_trl(i)+1));
    photo_cycle_start{i,1}=zeros(1,(Cycle_trl(i)+1)); % the last one always is
    photo_cycle_stop{i,1}=zeros(1,(Cycle_trl(i)+1));
    if isfield(input, 'tStimOffTimes')
    offs_trl = double(input.tStimOffTimes{i})./1000; % off times for each trial, the second last one is the one before target on
    duration=int64(Fs.*(sum(offs_trl)+stimOnTimeS*(Cycle_trl(i)+1))); % for the duration of each trial
    elseif isfield(input,'tStimOffTimeMs')
    offs_trl = double(input.tStimOffTimeMs{i})./1000; % off times for each trial, the second last one is the one before target on
    duration=int64(Fs.*(offs_trl+stimOnTimeS)*(Cycle_trl(i)+1)); % for the duration of each trial     
    else
    offs_trl = double(input.stimOffTimeMs)./1000; % off times for each trial, the second last one is the one before target on
    duration=int64(Fs.*(offs_trl+stimOnTimeS)*(Cycle_trl(i)+1)); 
    
    end 
    try
        labjack{i,:}=labjack_threshold(trialstart(i)-3000:trialstart(i)+duration); % 3000 is the baseline...
        photo{i,:}=photo_clean(trialstart(i)-3000:trialstart(i)+duration); %use arryfun later
        
    catch
        
        warning('trial cycles not complete');
        
        labjack{i,:}=zeros(1,3000+duration); % 3000 is the baseline...
        photo{i,:}=zeros(1,3000+duration);
    end
    if input.doBlock2 ==1
     try
        Led{i,:}=Led_threshold(trialstart(i)-3000:trialstart(i)+duration); % 3000 is the baseline...
        
    catch
        warning('trial cycles not complete');
        Led{i,:}=zeros(1,3000+duration); % 3000 is the baseline...
      
    end    
        
    end
    
    index_cycle=[];
    diff_cycle=[];
    c=[];
    index_cycle=find(photo{i,:}==1);
    diff_cycle=diff(index_cycle);
    if length(Offs)>3 % deal with pair pulse condition
    
    c=intersect(find(diff_cycle>=6000),find(diff_cycle<=4*Fs+6000))+1;
        
    else
    c=intersect(find(diff_cycle>=6000),find(diff_cycle<=52000))+1;
    end
    c_new = c(1:Cycle_trl(i));
    photo_cycle_start{i,1}(1)= index_cycle(1);
    photo_cycle_start{i,1}(2:(Cycle_trl(i)+1))=index_cycle(c_new);
    photo_cycle_stop{i,1}(1:Cycle_trl(i))=index_cycle(c_new-1);
    photo_cycle_stop{i,1}((Cycle_trl(i)+1))=index_cycle(end);
    photo_duration{i,1}=photo_cycle_stop{i,1}-photo_cycle_start{i,1};
    photo_inter{i,1}=photo_cycle_start{i,1}(2:(Cycle_trl(i)+1))-photo_cycle_stop{i,1}(1:Cycle_trl(i));
    photo_duration_MS{i,1}=photo_duration{i,1}.*1000/Fs;
    photo_inter_MS{i,1}=photo_inter{i,1}.*1000./Fs;
    
    trial_duration=0;
    
    for j=1:(Cycle_trl(i)+1)
        Time_start{i,1}(j)=trialstart(i)+trial_duration;
        Time_start_new{i,1}(j) =Time_start{i,1}(j) +(compensate.*Fs./1000); % compensate the difference between photodiode and visual stimus signal
        
        
        
        if j<=Cycle_trl(i)
            trial_duration=trial_duration + photo_duration{i,1}(j)+photo_inter{i,1}(j);
        end
        
        
        Time_stop{i,1}(j)=Time_start{i,1}(j)+photo_duration{i,1}(j);
        Time_stop_new{i,1}(j) = Time_stop{i,1}(j)+(compensate.*Fs./1000);
    end
    
    
    
    
    
    
end
%% get the LED time stamp relative to the visual onset 
if input.doBlock2 ==1 && input.block2DoTrialLaser == 1
block2 = cell2mat(input.tBlock2TrialNumber); % one is led trials

temp = Led(block2==1,:);
temp1 = photo(block2==1,:);

idx = cellfun(@(x) find(x==1),temp,'UniformOutput',false);
idx1 = cellfun(@(x) find(x==1),temp1,'UniformOutput',false);
time = cell2mat(cellfun(@(x,y) (y(1,1)-x(1,1))./Fs,idx,idx1,'UniformOutput',false));
% calculat led duration
led_dura = cell2mat(cellfun(@(x) sum(x==1)./Fs,temp,'UniformOutput',false));
LED_visual(1) = -mean(time);
LED_visual(2)= LED_visual(1)  + mean(led_dura); % duration is 200ms 
LED_visual = LED_visual - compensate./1000;
else
    LED_visual = NaN;
end





%% check if it is correct sort out non led condition
trial=25;
 x = (1:length(labjack{trial,:}))./Fs ; 
figure

subplot(2,1,1)
plot(x,labjack{trial,:},'r')

subplot(2,1,2)
plot(x,photo{trial,:},'g')
if input.doBlock2==1
    hold on
plot(x,Led{trial,:},'b')
end
supertitle('representative trial: red: labjack trigger; green: photodiode signal')

%% check if it is correct

figure
vline(Time_start_new{trial,1},'k')
hold on
vline(Time_stop_new{trial,1},'r')




end