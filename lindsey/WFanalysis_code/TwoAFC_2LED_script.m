mouse = 'i560';
date = '170810';
time = strvcat('1643');
nrun = 2;
mkdir(['Z:\home\lindsey\Analysis\Widefield_imaging\' mouse '\' date '_' mouse])
offset = [];
data = [];
temp = [];
for irun = 1:nrun
    if size(time,1)==nrun
        fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    else
        fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
    end
    load(fName);    
    
    pName = fullfile('Z:\home\lindsey\Data\Widefield_images', [date '_' mouse], [date '_' mouse '_' num2str(irun)], [date '_' mouse '_' num2str(irun) '_MMStack.ome.tif']);
    data_temp = readtiff(pName,'single');
    sz(irun) = size(data_temp,3);
    
    if isfield(input, 'cStimOn') 
        if size(time,1)>1
            if irun>1
                ntrials = size(input.trialOutcomeCell,2);
                for itrial = 1:ntrials
                    temp(irun).cItiStart{itrial} = temp(irun).cItiStart{itrial}+offset;
                    temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset;
                    temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                    temp(irun).cStimOff{itrial} = temp(irun).cStimOff{itrial}+offset;
                    temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset;
                end
            end
        end
    end
    offset = offset+sz(irun);
    data = cat(3,data,data_temp);
end
clear data_temp

figure; 
for i = 1:6
    subplot(2,3,i); 
    imagesc(data(:,:,i));
    clim([0 216])
end

cItiStart = cell2mat(input.cItiStart);
trial_rem = max(find(cItiStart<sz(1)),[],2);

data_470 = data(:,:,1:2:end);
data_530 = data(:,:,2:2:end);
siz = size(data_470);
data_470_long = reshape(data_470, [siz(1)*siz(2) siz(3)]);
data_530_long = reshape(data_530, [siz(1)*siz(2) size(data_530,3)]);
clear data_470 data_530
data_470_interp = interp1(1:2:(siz(3)*2), double(data_470_long)',1:(siz(3)*2));
data_470_interp = reshape(data_470_interp',[siz(1) siz(2) siz(3)*2]);
data_530_interp = interp1(2:2:(size(data_530_long,2)*2), double(data_530_long)',1:(siz(3)*2));
data_530_interp = reshape(data_530_interp',[siz(1) siz(2) siz(3)*2]);
clear data_470_long data_530_long


data_sub = data_470_interp - (1.5.*data_530_interp);

cStimOn = cell2mat(input.cStimOn);
nTrials = size(cStimOn,2);
SIx = strcmp(input.trialOutcomeCell, 'success');
IIx = strcmp(input.trialOutcomeCell, 'incorrect');
tLeftTrial = cell2mat(input.tLeftTrial);
tGratingContrast = cell2mat(input.tGratingContrast);
cons = unique(tGratingContrast);
ncon = length(cons);
data_trial = nan(siz(1),siz(2),30,nTrials);

for itrial = 1:nTrials
    if trial_rem ~= itrial
        if size(data_sub,3) > cStimOn(itrial)+20 
            data_trial(:,:,:,itrial) = data_sub(:,:,cStimOn(itrial)-9:cStimOn(itrial)+20);
        end
    end
end

data_f = mean(data_trial(:,:,1:10,:),3);
data_trial_dfof = bsxfun(@rdivide, bsxfun(@minus,data_trial,data_f),data_f);

figure;
start = 1;
for ii = 1:2
    for i = 1:ncon
    indR = intersect(intersect(find(SIx), find(tGratingContrast == cons(i))),find(tLeftTrial==ii-1));
    data_trial_R1 = mean(nanmean(data_trial_dfof(:,:,13:15,indR),4),3);
    subplot(2,3,start)
    imagesc(data_trial_R1)
    title([num2str(cons(i)) '; ' num2str(length(indR)) ' trials'])    
    clim([0 .03])
    start = start+ 1;
    colormap gray
    end
end
suptitle([mouse ' ' date ' Top: Contra; Bottom: Ipsi'])
print(['Z:\home\lindsey\Analysis\Widefield_imaging\' mouse '\' date '_' mouse '\' date '_' mouse '_respByConByPos.pdf'], '-dpdf', '-bestfit')

figure
start = 1;
for i = 1:ncon
    indR = intersect(intersect(find(SIx), find(tGratingContrast == cons(i))),find(tLeftTrial==0));
    data_trial_R1 = mean(nanmean(data_trial_dfof(:,:,13:15,indR),4),3);
    subplot(2,3,start)
    imagesc(data_trial_R1)
    title([num2str(cons(i)) '; ' num2str(length(indR)) ' trials'])
    clim([0 .03])
    start = start+ 1;
    colormap gray
end
for i = 1:ncon
    indR = intersect(intersect(find(IIx), find(tGratingContrast == cons(i))),find(tLeftTrial==0));
    data_trial_R1 = mean(nanmean(data_trial_dfof(:,:,13:15,indR),4),3);
    subplot(2,3,start)
    imagesc(data_trial_R1)
    title([num2str(cons(i)) '; ' num2str(length(indR)) ' trials'])
    clim([0 .03])
    start = start+ 1;
    colormap gray
end
suptitle([mouse ' ' date ' Contra Stim; Top: Correct; Bottom: Incorrect'])
print(['Z:\home\lindsey\Analysis\Widefield_imaging\' mouse '\' date '_' mouse '\' date '_' mouse '_contraRespByConByOutcome.pdf'], '-dpdf', '-bestfit')

figure
start = 1;
for i = 1:ncon
    indR = intersect(intersect(find(SIx), find(tGratingContrast == cons(i))),find(tLeftTrial==1));
    data_trial_R1 = mean(nanmean(data_trial_dfof(:,:,13:15,indR),4),3);
    subplot(2,3,start)
    imagesc(data_trial_R1)
    title([num2str(cons(i)) '; ' num2str(length(indR)) ' trials'])
    clim([0 .03])
    start = start+ 1;
    colormap gray
end
for i = 1:ncon
    indR = intersect(intersect(find(IIx), find(tGratingContrast == cons(i))),find(tLeftTrial==1));
    data_trial_R1 = mean(nanmean(data_trial_dfof(:,:,13:15,indR),4),3);
    subplot(2,3,start)
    imagesc(data_trial_R1)
    title([num2str(cons(i)) '; ' num2str(length(indR)) ' trials'])
    clim([0 .03])
    start = start+ 1;
    colormap gray
end
suptitle([mouse ' ' date ' Ipsi Stim; Top: Correct; Bottom: Incorrect'])
print(['Z:\home\lindsey\Analysis\Widefield_imaging\' mouse '\' date '_' mouse '\' date '_' mouse '_IpsiRespByConByOutcome.pdf'], '-dpdf', '-bestfit')
