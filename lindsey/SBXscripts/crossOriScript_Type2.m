clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirType2_ExptList';
iexp = 20; 
eval(ds)

frame_rate = params.frameRate;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(ImgFolder, nrun);

base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];


fprintf(['2P imaging cross ori analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load data
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
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
    
    if isfield(input, 'cStimOneOn') 
        if irun>1
            ntrials = size(input.cStimOneOn,2);
            for itrial = 1:ntrials
                temp(irun).cStimOneOn{itrial} = temp(irun).cStimOneOn{itrial}+offset;
                temp(irun).cStimOneOff{itrial} = temp(irun).cStimOneOff{itrial}+offset;
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateStructuresLG(temp);
clear data_temp
clear temp
toc

%% Register data

if exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg] = stackRegister_MA(data,[],[],out);
    nep = 9;
else
    totframes = size(data,3);
    nep = 9;
    nframes = 500; %nframes to average for target
    nskip = floor(totframes./nep); %nframes to skip for each average
    [n, n2] = subplotn(nep); 
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:nep 
        subplot(n,n2,i); 
        imagesc(mean(data(:,:,1+((i-1)*nskip):nframes+((i-1)*nskip)),3)); 
        title([num2str(1+((i-1)*nskip)) '-' num2str(nframes+((i-1)*nskip))]); 
    end
    
    %    b. GUI to select target image- choose one that is sharp and close to center of stack
    f=gcf;
    w = waitforbuttonpress; %click on subplot
    if w == 0
        axesClicked = gca;
        allAxes = flipud(findobj(f.Children,'Type','axes'));
        numClicked = find(axesClicked==allAxes);
        close all
    end
    fprintf(['Selected subplot ' num2str(numClicked)])
    %    c. Create target image
    data_avg = mean(data(:,:,1+((numClicked-1)*nskip):nframes+((numClicked-1)*nskip)),3); %average 500 frames to make target
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end

data_reg_avg = mean(data_reg,3);
% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*nskip):500+((i-1)*nskip)),3)); title([num2str(1+((i-1)*nskip)) '-' num2str(500+((i-1)*nskip))]); end
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
movegui('center')
figure; imagesq(mean(data_reg(:,:,1:nskip),3)); truesize;
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells
clear data out
[cStimOn cStimOff] = photoFrameFinder_Sanworks(info.frame);
nTrials = length(cStimOn);
if length(cStimOn)>length(cStimOff)
    nTrials = length(cStimOff);
    cStimOn = cStimOn(1:nTrials);
end
if size(input.cStimOneOn,2)>length(cStimOn)
   input = trialChopper(input,1:length(cStimOn));
end
sz = size(data_reg);

data_resp = nan(sz(1),sz(2),nTrials);
data_f = nan(sz(1),sz(2),nTrials);

for itrial = 1:nTrials
    if cStimOn(itrial) + 20 < sz(3)
        data_resp(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+20),3);
    end
    data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-15:cStimOn(itrial)),3);
end

data_resp_dfof = (data_resp-data_f)./data_f;

tPlaid = celleqel2mat_padded(input.tMaskOneGratingContrast);
stimDir_all = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
stimDirs = unique(stimDir_all);
nStimDir = length(stimDirs);
stimTF_all = celleqel2mat_padded(input.tStimOneGratingTemporalFreqCPS);
maskTF_all = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
TFs = unique(stimTF_all);
nTF = length(TFs);

if nTF>1
    nStim = nStimDir*nTF;
else
    nStim  = nStimDir*2;
end
data_dfof = nan(sz(1),sz(2),nStim);
ind_grating = find(tPlaid==0);
start = 1;
for id = 1:nStimDir
    ind_dir = find(stimDir_all == stimDirs(id));
    for it = 1:nTF
        ind_tf = find(stimTF_all == TFs(it));
        ind_use = intersect(ind_tf,intersect(ind_dir,ind_grating));
        data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_use),3);
        start = start+1;
        if nTF == 1
            ind_use = intersect(find(tPlaid),ind_dir);
            data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_use),3);
            start = start+1;
        end
    end
end

data_dfof(:,:,isnan(mean(mean(data_dfof,1),2))) = [];
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
data_dfof = cat(3,data_dfof_max,cat(3,data_dfof_max,data_dfof));
figure;
movegui('center')
imagesc(data_dfof_max)
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_maxdfof.pdf']),'-dpdf')

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'tPlaid','stimDir_all', 'stimDirs', 'nStimDir', 'stimTF_all','maskTF_all', 'TFs', 'nTF', 'frame_rate', 'nTrials')

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)   
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end

mask_cell= bwlabel(mask_all);
if doRedChannel
    red_cells = unique(mask_cell(find(red_mask)));
else
    red_cells = [];
end
figure; movegui('center')
imagesc(mask_cell)    

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg data_resp data_resp_dfof bwout

 %% neuropil subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_cells')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 


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
clear data_reg data_reg_down

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np

if strcmp(expt(iexp).img_strct,'dendrites')
    figure; image(imShade(data_reg_avg,mask_cell))
    truesize
    print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellMaskFOV.pdf']),'-dpdf','-bestfit')
end
 