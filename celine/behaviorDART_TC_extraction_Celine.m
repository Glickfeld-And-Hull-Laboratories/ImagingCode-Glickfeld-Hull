%% get path names

close all;clear all;clc;
ds = 'DART_behavior_ExptList';
experimentFolder = 'SST_behavior';
iexp = 11; 
eval(ds)

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).behavior_runs;
time = expt(iexp).behavior_time;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    database = fullfile('/All_Staff/home/ACh/Data/2p_data');
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis/',experimentFolder);
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/data-');
elseif string(hostname) == 'NB-NUKE'
    isilonName = 'Z:/All_Staff';
    base = fullfile(isilonName,'/home/ACh/Analysis/2p_analysis',experimentFolder);
    database = fullfile(isilonName,'/home/ACh/Data/2p_data/');
else
    isilonName = 'duhs-user-nc1.dhe.duke.edu/';
    base = fullfile('/home/ACh/Analysis/2p_analysis',experimentFolder);
    database = fullfile('/home/ACh/Data/2p_data/');
   
   beh_prefix = strcat('Z:\Behavior\Data\data-');
end

fprintf(['2P imaging DART behavior analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

fnout = fullfile(base,mouse,date,cell2mat(ImgFolder));
mkdir(fnout)
%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(database,mouse,date,ImgFolder{irun});
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

%% Choose register interval
step = 5000;
nep = floor(size(data,3)./step);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*step):500+((i-1)*step)),3)); title([num2str(1+((i-1)*step)) '-' num2str(500+((i-1)*step))]); end

data_avg = mean(data(:,:,55000:55500),3);
%% Register data
[database mouse '\' date '\' ImgFolder{irun}]
if exist(fullfile(fnout,'\regOuts&Img.mat'))
    load(fullfile(fnout,'regOuts&Img.mat'))
    save(fullfile(fnout,'input.mat'), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(fnout))
    save(fullfile(fnout,'regOuts&Img.mat'), 'out', 'data_avg')
    save(fullfile(fnout,'input.mat'), 'input')
end
clear data out


%% find activated cells

cTarget = celleqel2mat_padded(input.cTargetOn);
cStart = celleqel2mat_padded(input.cLeverDown);
nTrials = length(cTarget);
sz = size(data_reg);
data_f_base = zeros(sz(1),sz(2),nTrials);
data_base = nan(sz(1),sz(2),nTrials);
data_f_targ = zeros(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials    
    if ~isnan(cStart(itrial))
        if cStart(itrial)+19 < sz(3) % making a buffer of enough frames prior to the end of the trial. why 19?
            data_f_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
            data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+5:cStart(itrial)+10),3);
        end
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+19 < sz(3)
            data_f_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)-20:cTarget(itrial)-1),3);
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+5:cTarget(itrial)+10),3);
        end
    end
    
end
data_base_dfof = (data_base-data_f_base)./data_f_base;
data_targ_dfof = (data_targ-data_f_targ)./data_f_targ;
trialCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(trialCon);
nCon = length(cons);
data_dfof_con = zeros(sz(1),sz(2),nCon);
[n n2] = subplotn(nCon);
for iCon = 1:nCon
    ind_con = find(trialCon == cons(iCon));
    data_dfof_con(:,:,iCon) = nanmean(data_targ_dfof(:,:,ind_con),3);
    subplot(n,n2,iCon)
    imagesc(data_dfof_con(:,:,iCon))
    title(cons(iCon))
end


data_dfof = cat(3, reshape(data_dfof_con,[sz(1), sz(2), nCon]), data_dfof_con);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)
data_dfof = cat(3, data_dfof, data_dfof_max);

down = 100;
data_reg_down  = stackGroupProject(data_reg,down);
pixelCorr = getPixelCorrelationImage(data_reg_down);
figure; 
imagesc(pixelCorr)


save(fullfile(fnout, 'cellSelect.mat'), 'data_dfof')
%% load red data
data_reg_avg =  mean(data_reg(:,:,:),3);
%this is where we use the 1040, 1000-frame run
if exist(fullfile(fnout,'redImage.mat'))
    load(fullfile(fnout,'redImage'))
elseif ~isempty(expt(iexp).redChannelRun) %if there IS a red channel run, find and load it
    redRun = expt(iexp).redChannelRun;
    imgMatFile = [redRun '_000_000.mat'];
    cd(fullfile(database,mouse,date,redRun))
    load(imgMatFile);

    fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
    if size(data_temp,1) == 2
        data_rg = squeeze(data_temp(1,:,:,:));
        data_rr = squeeze(data_temp(2,:,:,:));
    else
        data_rr = squeeze(data_temp(1,:,:,:));
    end
    clear data_temp

    if exist('redChImg')
        [out, data_rr_reg] = stackRegister(stackGroupProject(data_rr,100), redChImg);
        redChImg = mean(data_rr_reg,3);
    elseif info.config.pmt0_gain>0.5 %if there is a green channel in this run, it gets registered to the registration image from green channel from the 920 run
        %data_rr = padarray(data_rr,9,0,'pre');
        redAvg = mean(data_rr,3);%take the average of the red at 1040
        [out, data_rr_reg] = stackRegister(data_rr,redAvg); %register red at 1040 to average of itself
        [~, data_rg_reg] = stackRegister_MA(data_rg,[],[],out); %register green at 1040 using the shifts from the registration of the red at 1040
        redChImgTemp = mean(data_rr_reg,3);  %take an average of the registered red at 1040
        rg_avg = mean(data_rg_reg,3); %take an average of the registered green at 1040
        [out2, ~] = stackRegister(rg_avg,data_reg_avg); %register that green average to the registered data from my experimental run
        [~,redChImg]=stackRegister_MA(redChImgTemp,[],[],out2); %apply those shifts to the mean registered red at 1040
        
%         [out, data_rg_reg] = stackRegister(data_rg,data_avg); %register the green channel from the 1040 run to the green channel from the 920 run
%         [~, data_rr_reg]=stackRegister_MA(data_rr,[],[],out); %use those shifts to register the red 1040 run
%         redChImg = mean(data_rr_reg,3);
        
        
    else %if there is no green channel in this run
        redAvg = mean(data_rr,3);
        [out, data_rr_reg] = stackRegister(data_rr,redAvg);
        redChImgTemp = mean(data_rr_reg,3);
        [~,redChImg] = stackRegister(redChImgTemp,data_reg_avg);
    end
    
    figure; colormap gray; imagesc(redChImg);  movegui('center');title('registration image for red channel');
   
    rgb = zeros(sz(1),sz(2),3);
    rgb(:,:,1) = redChImg./max(redChImg(:));
    rgb(:,:,2) = data_reg_avg./max(data_reg_avg(:));
    figure; image(rgb);  movegui('center')
    title('Green-920 + Red-1040')
    print(fullfile(fnout,'red_green_FOV.pdf'),'-dpdf','-bestfit')
    save(fullfile(fnout,'redImage'),'redChImg')

elseif ~exist('redChImg')
    redChImg = zeros(size(regImg));
end

clear data_rr data_rg data_rg_reg data_rr_reg

%% cell segmentation 
close all
data_dfof = cat(3, data_dfof, data_reg_avg,data_reg_avg,data_reg_avg, pixelCorr);
redForSegmenting = cat(3, redChImg,redChImg,redChImg); %make a dataframe that repeats the red channel image twice
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
%find and label the red cells - this is the first segmentation figure that
%comes up
if ~isempty(expt(iexp).redChannelRun)
    for iStim=1:size(redForSegmenting,3)
        mask_data_temp=redForSegmenting(:,:,iStim);
        mask_data_temp(find(mask_exp >= 1)) = 0;
        bwout = imCellEditInteractiveLG(mask_data_temp);
        mask_all = mask_all+bwout; 
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        close all
    end
end

mask_cell_red = bwlabel(mask_all);
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)

nCells = max(mask_cell(:));
%mask_label = ones(1,nCells); %this is for the EMX and similar lines only
mask_label = zeros(1,nCells);
for i = 1:nCells
    if mask_cell_red(find(mask_cell == i, 1))
        mask_label(1,i) = 1;
    end
end

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(fnout,'mask_cell.mat'), 'data_dfof', 'mask_cell', 'mask_np','mask_label','data_reg_avg')

%clear data_Avg data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
%clear data_reg data_reg_down

save(fullfile(fnout,'TCs.mat'), 'data_tc', 'np_tc', 'npSub_tc')



%% finding trial types
%find trials within nframes
lastTrial = find(cell2mat(input.cTrialEnd)<min(nframes),1,'last');

%these are indices of trials
hit =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'success'));
FA =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'failure'));
miss =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'ignore'));
block2trials = find(cell2mat(input.tBlock2TrialNumber(1:lastTrial)));
block1trials = find(~cell2mat(input.tBlock2TrialNumber(1:lastTrial)));

hitRate = nan(nCon,2);

%get the mean neural response to hits and misses at each contrast for each
%block. Currently looking 1 second before through 1 second after the stim
%turns on
trialAvrg_tcs = nan(60,nCells,nCon,4); %empty matrix, 60 frames by nCells by nCon by hit B1, miss B1, hit B2, miss B2
for iCon = 1:nCon
    conInds=find(trialCon==cons(iCon)); %find trials for this contrast
    block1Inds = intersect(conInds,block1trials); %find trials for block 1 at this contrast
    hitInds = intersect(hit, block1Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss,block1Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds)); %find when the target comes on for the hit trials for this contrast, in frames
    missTargetTimes = cell2mat(input.cTargetOn(missInds)); %find when the target comes on for the miss trials for this contrast, in frames
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitInds) %loop through the hit trials for this contrast
        hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(hitTargetTimes(iHit)-15:hitTargetTimes(iHit),:),1);
        %to get dfof, take the np subtracted timecourse for t=-1 before stim through 1
        %second, and subtract the mean of the -0.5 through 0 (using the
        %last half second before stim onset as the baseline F
    end
    trialAvrg_tcs(:,:,iCon,1)=mean(hitsTemp,3);
    for iMiss = 1:length(missInds) %loop through the miss trials for this contrast
        missTemp(:,:,iMiss)=npSub_tc(missTargetTimes(iMiss)-29:missTargetTimes(iMiss)+30,:)-mean(npSub_tc(missTargetTimes(iMiss)-15:missTargetTimes(iMiss),:),1);
    end
    trialAvrg_tcs(:,:,iCon,2)=mean(missTemp,3);
    hitRate(iCon, 1) = length(hitInds)/(length(hitInds)+length(missInds));

    %now do block 2
    block2Inds = intersect(conInds,block2trials); %find trials for block 1 at this contrast
    hitInds = intersect(hit, block2Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss,block2Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds));
    missTargetTimes = cell2mat(input.cTargetOn(missInds));
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitInds)
        hitsTemp(:,:,iHit)=npSub_tc(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc(hitTargetTimes(iHit)-15:hitTargetTimes(iHit),:),1);
    end
    trialAvrg_tcs(:,:,iCon,3)=mean(hitsTemp,3);
    for iMiss = 1:length(missInds)
        missTemp(:,:,iMiss)=npSub_tc(missTargetTimes(iMiss)-29:missTargetTimes(iMiss)+30,:)-mean(npSub_tc(missTargetTimes(iMiss)-15:missTargetTimes(iMiss),:),1);
    end
    trialAvrg_tcs(:,:,iCon,4)=mean(missTemp,3);
    hitRate(iCon, 2) = length(hitInds)/(length(hitInds)+length(missInds));

end
save(fullfile(fnout,'trialAvrg_tcs.mat'),'trialAvrg_tcs')

%% plot mean responses
pyr = ~mask_label;
sst = mask_label;

t=(-29:30)/30;
positions=[1,3,5,7,9,11];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs(:,pyr,iCon,1),2);
    stdHit=std(trialAvrg_tcs(:,pyr,iCon,1),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,pyr,iCon,2),2);
    stdMiss=std(trialAvrg_tcs(:,pyr,iCon,2),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(['Block1, contrast ' num2str(cons(iCon))])
    box off
    ylim([-30 40])
    txt = {num2str(hitRate(iCon,1))};
    text(-0.75,30,txt)
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanHit=mean(trialAvrg_tcs(:,pyr,iCon,3),2);
    stdHit=std(trialAvrg_tcs(:,pyr,iCon,3),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,pyr,iCon,4),2);
    stdMiss=std(trialAvrg_tcs(:,pyr,iCon,4),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(['Block2, contrast ' num2str(cons(iCon))])
    box off
    ylim([-30 40])
    txt = {num2str(hitRate(iCon,2))};
    text(-0.75,30,txt)
end
sgtitle('Pyr cells')

x0=5;
y0=0;
width=6;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'Pyr_Active_HitvsMiss.pdf'),'-dpdf', '-bestfit')

%% plot mean responses

sst = logical(mask_label);

t=(-29:30)/30;
positions=[1,3,5,7,9,11];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs(:,sst,iCon,1),2);
    stdHit=std(trialAvrg_tcs(:,sst,iCon,1),[],2);
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,sst,iCon,2),2);
    stdMiss=std(trialAvrg_tcs(:,sst,iCon,2),[],2);
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(['Block1, contrast ' num2str(cons(iCon))])
    box off
    ylim([-30 40])
    txt = {num2str(hitRate(iCon,1))};
    text(-0.75,30,txt)
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanHit=mean(trialAvrg_tcs(:,sst,iCon,3),2);
    stdHit=std(trialAvrg_tcs(:,sst,iCon,3),[],2);
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,sst,iCon,4),2);
    stdMiss=std(trialAvrg_tcs(:,sst,iCon,4),[],2);
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(['Block2, contrast ' num2str(cons(iCon))])
    box off
    ylim([-30 40])
    txt = {num2str(hitRate(iCon,2))};
    text(-0.75,30,txt)
end
sgtitle('SST cells')

x0=5;
y0=0;
width=6;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
print(fullfile(fnout,'SST_Active_HitvsMiss.pdf'),'-dpdf', '-bestfit')

%% read in the passive dataset
ImgFolderPassive = expt(iexp).passive_runs;
time = expt(iexp).passive_time;
nrun = length(ImgFolderPassive);
run_str_passive = catRunName(cell2mat(ImgFolderPassive), nrun);

tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(database,mouse,date,ImgFolderPassive{irun});
    cd(CD);
    imgMatFile = [ImgFolderPassive{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc
%% register the passive data to the active data
fnout_passive = fullfile(base,mouse,date,cell2mat(ImgFolderPassive));
[database mouse '\' date '\' ImgFolderPassive{irun}]
if exist(fullfile(fnout_passive,'\regOuts&Img.mat'))
    load(fullfile(fnout_passive,'regOuts&Img.mat'))
    save(fullfile(fnout_passive,'input.mat'), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(fnout_passive))
    save(fullfile(fnout_passive,'regOuts&Img.mat'), 'out', 'data_avg')
    save(fullfile(fnout_passive,'input.mat'), 'input')
end
clear data out

figure;
subplot(2,1,1)
imagesc(data_avg);
title('Active run');
subplot(2,1,2)
imagesc(mean(data_reg(:,:,:),3));
title('Passive run')
%% apply the masks from the active run
% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc_passive = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down_passive  = stackGroupProject(data_reg,down);
data_tc_down_passive = stackGetTimeCourses(data_reg_down_passive, mask_cell);
np_tc_passive = zeros(sz(3),nCells);
np_tc_down_passive = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc_passive(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down_passive(:,i) = stackGetTimeCourses(data_reg_down_passive,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down_passive-tcRemoveDC(np_tc_down_passive*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w_passive = 0.01*ind;
npSub_tc_passive = data_tc_passive-bsxfun(@times,tcRemoveDC(np_tc_passive),np_w_passive);
%clear data_reg data_reg_down

save(fullfile(fnout,'TCs_passive.mat'), 'data_tc_passive', 'np_tc_passive', 'npSub_tc_passive')

clear data_tc_passive data_tc_down_passive np_tc_passive np_tc_down_passive mask_np 

%% finding trial types
%find trials within nframes
trialCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(trialCon);
nCon = length(cons);
lastTrial = find(cell2mat(input.cTrialEnd)<min(nframes),1,'last');

%these are indices of trials
hit_passive =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'success'));
FA_passive =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'failure'));
miss_passive =  find(strcmp(input.trialOutcomeCell(1:lastTrial), 'ignore'));
block2trials_passive = find(cell2mat(input.tBlock2TrialNumber(1:lastTrial)));
block1trials_passive = find(~cell2mat(input.tBlock2TrialNumber(1:lastTrial)));

hitRate_passive = nan(nCon,2);

%get the mean neural response to hits and misses at each contrast for each
%block. Currently looking 1 second before through 1 second after the stim
%turns on
trialAvrg_tcs_passive = nan(60,nCells,nCon,4); %empty matrix, 60 frames by nCells by nCon by hit B1, miss B1, hit B2, miss B2
for iCon = 1:nCon
    conInds=find(trialCon==cons(iCon)); %find trials for this contrast
    block1Inds = intersect(conInds,block1trials_passive); %find trials for block 1 at this contrast
    hitInds = intersect(hit_passive, block1Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss_passive,block1Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds));
    hitTargetTimes(hitTargetTimes>(min(nframes)-30))=[];
    missTargetTimes = cell2mat(input.cTargetOn(missInds));
    missTargetTimes(missTargetTimes>(min(nframes)-30))=[];
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitTargetTimes)
        hitsTemp(:,:,iHit)=npSub_tc_passive(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc_passive(hitTargetTimes(iHit)-29:hitTargetTimes(iHit),:),1);
    end
    trialAvrg_tcs_passive(:,:,iCon,1)=mean(hitsTemp,3);
    for iMiss = 1:length(missTargetTimes)
        missTemp(:,:,iMiss)=npSub_tc_passive(missTargetTimes(iMiss)-29:missTargetTimes(iMiss)+30,:)-mean(npSub_tc_passive(missTargetTimes(iMiss)-29:missTargetTimes(iMiss),:),1);
    end
    trialAvrg_tcs_passive(:,:,iCon,2)=mean(missTemp,3);
    hitRate_passive(iCon, 1) = length(hitInds)/(length(hitInds)+length(missInds));

    %now do block 2
    block2Inds = intersect(conInds,block2trials_passive); %find trials for block 1 at this contrast
    hitInds = intersect(hit_passive, block2Inds); %hit trials for block 1 at this contrast
    missInds = intersect(miss_passive,block2Inds); %miss trials for block1 at this contrast
    hitTargetTimes = cell2mat(input.cTargetOn(hitInds));
    hitTargetTimes(hitTargetTimes>(min(nframes)-30))=[];
    missTargetTimes = cell2mat(input.cTargetOn(missInds));
    missTargetTimes(missTargetTimes>(min(nframes)-30))=[];
    hitsTemp = nan(60,nCells,length(hitInds));
    missTemp = nan(60,nCells,length(missInds));
    for iHit = 1:length(hitTargetTimes)
        hitsTemp(:,:,iHit)=npSub_tc_passive(hitTargetTimes(iHit)-29:hitTargetTimes(iHit)+30,:)-mean(npSub_tc_passive(hitTargetTimes(iHit)-29:hitTargetTimes(iHit),:),1);
    end
    trialAvrg_tcs_passive(:,:,iCon,3)=mean(hitsTemp,3);
    for iMiss = 1:length(missTargetTimes)
        missTemp(:,:,iMiss)=npSub_tc_passive(missTargetTimes(iMiss)-29:missTargetTimes(iMiss)+30,:)-mean(npSub_tc_passive(missTargetTimes(iMiss)-29:missTargetTimes(iMiss),:),1);
    end
    trialAvrg_tcs_passive(:,:,iCon,4)=mean(missTemp,3);
    hitRate_passive(iCon, 2) = length(hitInds)/(length(hitInds)+length(missInds));

end
save(fullfile(fnout,'trialAvrg_tcs_passive.mat'),'trialAvrg_tcs_passive')

%% plot mean responses
pyr = ~mask_label;

t=(-29:30)/30;
positions=[1,3,5,7];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs_passive(:,pyr,iCon,1),2);
    stdHit=std(trialAvrg_tcs_passive(:,pyr,iCon,1),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs_passive(:,pyr,iCon,2),2);
    stdMiss=std(trialAvrg_tcs_passive(:,pyr,iCon,2),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
    % txt = {num2str(hitRate_passive(iCon,1))};
    % text(-0.75,75,txt)
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanHit=mean(trialAvrg_tcs_passive(:,pyr,iCon,3),2);
    stdHit=std(trialAvrg_tcs_passive(:,pyr,iCon,3),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs_passive(:,pyr,iCon,4),2);
    stdMiss=std(trialAvrg_tcs_passive(:,pyr,iCon,4),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
    % txt = {num2str(hitRate_passive(iCon,2))};
    % text(-0.75,75,txt)
end

%% plot responses comparing block 1 and 2
positions=[1,3,5,7];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs_passive(:,pyr,iCon,1),2);
    stdHit=std(trialAvrg_tcs_passive(:,pyr,iCon,1),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit,'b')
    hold on
    meanHit=mean(trialAvrg_tcs_passive(:,pyr,iCon,3),2);
    stdHit=std(trialAvrg_tcs_passive(:,pyr,iCon,3),[],2);
    seHit=stdHit./sqrt(sum(pyr));
    shadedErrorBar(t,meanHit,seHit,'g')

    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])

    txt = {num2str(hitRate_passive(iCon,1))};
    text(-0.75,85,txt,'Color','b')

    txt = {num2str(hitRate_passive(iCon,2))};
    text(-0.75,55,txt,'Color','g')
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanMiss=mean(trialAvrg_tcs_passive(:,pyr,iCon,2),2);
    stdMiss=std(trialAvrg_tcs_passive(:,pyr,iCon,2),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'b')
    hold on
    meanMiss=mean(trialAvrg_tcs_passive(:,pyr,iCon,4),2);
    stdMiss=std(trialAvrg_tcs_passive(:,pyr,iCon,4),[],2);
    seMiss=stdMiss./sqrt(sum(pyr));
    shadedErrorBar(t,meanMiss,seMiss,'g')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
end

print(fullfile(fnout,'Passive_B1vsB2.pdf'),'-dpdf', '-bestfit')
%% plot mean responses SST active
t=(-29:30)/30;
positions=[1,3,5,7,9,11];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs(:,logical(sst),iCon,1),2,'omitmissing');
    stdHit=std(trialAvrg_tcs(:,logical(sst),iCon,1),[],2,'omitmissing');
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,logical(sst),iCon,2),2,'omitmissing');
    stdMiss=std(trialAvrg_tcs(:,logical(sst),iCon,2),[],2,'omitmissing');
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
    % txt = {num2str(hitRate(iCon,1))};
    % text(-0.75,75,txt)
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanHit=mean(trialAvrg_tcs(:,logical(sst),iCon,3),2,'omitmissing');
    stdHit=std(trialAvrg_tcs(:,logical(sst),iCon,3),[],2,'omitmissing');
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs(:,logical(sst),iCon,4),2,'omitmissing');
    stdMiss=std(trialAvrg_tcs(:,logical(sst),iCon,4),[],2,'omitmissing');
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
    % txt = {num2str(hitRate(iCon,2))};
    % text(-0.75,75,txt)
end

%% plot mean responses SST passive


t=(-29:30)/30;
positions=[1,3,5,7];

figure;
for iCon = 1:nCon
    subplot(nCon,2,positions(iCon))
    meanHit=mean(trialAvrg_tcs_passive(:,logical(sst),iCon,1),2,'omitmissing');
    stdHit=std(trialAvrg_tcs_passive(:,logical(sst),iCon,1),[],2,'omitmissing');
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs_passive(:,logical(sst),iCon,2),2,'omitmissing');
    stdMiss=std(trialAvrg_tcs_passive(:,logical(sst),iCon,2),[],2,'omitmissing');
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    box off
    ylim([-20 200])
    % txt = {num2str(hitRate_passive(iCon,1))};
    % text(-0.75,75,txt)
end

for iCon = 1:nCon
    subplot(nCon,2,iCon*2)
    meanHit=mean(trialAvrg_tcs_passive(:,logical(sst),iCon,3),2,'omitmissing');
    stdHit=std(trialAvrg_tcs_passive(:,logical(sst),iCon,3),[],2,'omitmissing');
    seHit=stdHit./sqrt(sum(sst));
    shadedErrorBar(t,meanHit,seHit)
    hold on
    meanMiss=mean(trialAvrg_tcs_passive(:,logical(sst),iCon,4),2,'omitmissing');
    stdMiss=std(trialAvrg_tcs_passive(:,logical(sst),iCon,4),[],2,'omitmissing');
    seMiss=stdMiss./sqrt(sum(sst));
    shadedErrorBar(t,meanMiss,seMiss,'r')
    title(num2str(cons(iCon)))
    % box off
    ylim([-20 200])
    % txt = {num2str(hitRate_passive(iCon,2))};
    % text(-0.75,75,txt)
end