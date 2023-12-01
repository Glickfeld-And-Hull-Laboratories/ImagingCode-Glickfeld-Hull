clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriSingleStimRandDirAdapt_ExptList';
iexp = 5; 
rc = behavConstsAV;
eval(ds)

frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging retinotopy analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
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
    
    if isfield(input, 'cAdaptStart') 
        if irun>1
            ntrials = size(input.cAdaptStart,2);
            for itrial = 1:ntrials
                temp(irun).cAdaptStart{itrial} = temp(irun).cAdaptStart{itrial}+offset;
                temp(irun).cTestOn{itrial} = temp(irun).cTestOn{itrial}+offset;
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

% Choose register interval
step = 10000;
nep = floor(size(data,3)./step);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*step):500+((i-1)*step)),3)); title([num2str(1+((i-1)*step)) '-' num2str(500+((i-1)*step))]); colormap gray; clim([0 3000]); end

data_avg = mean(data(:,:,60001:60500),3);
%% Register data

[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')

figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*step):500+((i-1)*step)),3)); title([num2str(1+((i-1)*step)) '-' num2str(500+((i-1)*step))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% if red channel data
if doRedChannel
    CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolderRed(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolderRed(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    nframes = info.config.frames;
    fprintf(['Reading red data- ' num2str(nframes) ' frames \r\n'])
    data_red = sbxread(imgMatFile(1,1:11),0,nframes);
    data_red_g = squeeze(data_red(1,:,:,:));
    data_red_r = squeeze(data_red(2,:,:,:));

    [rg_out rg_reg] = stackRegister(data_red_g,data_avg);
    [rr_out rr_reg] = stackRegister_MA(data_red_r,[],[],rg_out);

    data_red_avg = mean(rr_reg,3);
    figure; imagesc(data_red_avg);
end

%% find activated cells
clear data out
cAdapt= celleqel2mat_padded(input.cAdaptStart);
cTest = celleqel2mat_padded(input.cTestOn);
nTrials = length(cTest);
sz = size(data_reg);

data_adapt = nan(sz(1),sz(2),nTrials);
data_test = nan(sz(1),sz(2),nTrials);
data_f = nan(sz(1),sz(2),nTrials);

ind_contadapt = find(celleqel2mat_padded(input.tDoAsynchronousAdapt)==0 &celleqel2mat_padded(input.tDoContingentAdapt)==1 ...
    & celleqel2mat_padded(input.tDoSingleStimAdapt)==0);
ind_asynadapt = find(celleqel2mat_padded(input.tDoAsynchronousAdapt)==1 &celleqel2mat_padded(input.tDoContingentAdapt)==0 ...
    & celleqel2mat_padded(input.tDoSingleStimAdapt)==0);
ind_singadapt = find(celleqel2mat_padded(input.tDoAsynchronousAdapt)==0 &celleqel2mat_padded(input.tDoContingentAdapt)==0 ...
    & celleqel2mat_padded(input.tDoSingleStimAdapt)==1);
ind_noadapt = find(celleqel2mat_padded(input.tDoAsynchronousAdapt)==0 &celleqel2mat_padded(input.tDoContingentAdapt)==0 ...
    & celleqel2mat_padded(input.tDoSingleStimAdapt)==0);
ind_noadapt_start = intersect(ind_noadapt, find(celleqel2mat_padded(input.tAdaptTimeMs) ==40000));

for itrial = 1:nTrials
    if cAdapt(itrial) + 30 < sz(3)
        data_adapt(:,:,itrial) = mean(data_reg(:,:,cAdapt(itrial)+5:cAdapt(itrial)+30),3);
    end
    if cTest(itrial) + 20 < sz(3)
        data_test(:,:,itrial) = mean(data_reg(:,:,cTest(itrial)+5:cTest(itrial)+20),3);
    end
    data_f(:,:,itrial) = mean(data_reg(:,:,cTest(itrial)-15:cTest(itrial)),3);
end

data_adapt_dfof = (data_adapt-data_f)./data_f;
data_test_dfof = (data_test-data_f)./data_f;

data_adapt_avg = zeros(sz(1),sz(2),3);
n = unique(input.trialsPerAdaptBlock);
for i = 1:3
    x = [];
    for ii = 1:20
        x = [x; (1+((i-1).*n):n+((i-1).*n))'];
    end
    ind = intersect(1:nTrials, x);
    data_adapt_avg(:,:,i) = nanmean(data_adapt_dfof(:,:,ind),3);
end

testCon = celleqel2mat_padded(input.tTestStimGratingContrast);
maskCon = celleqel2mat_padded(input.tTestMaskGratingContrast);
testCons = unique(testCon);
maskCons = unique(maskCon);
nMaskCons = length(maskCons);
nTestCons = length(testCons);

testDir = celleqel2mat_padded(input.tTestStimGratingDirectionDeg);
maskDir = celleqel2mat_padded(input.tTestMaskGratingDirectionDeg);
testDirs = unique(testDir);
maskDirs = unique(maskDir);
nTestDirs = length(testDirs);           
nMaskDirs = length(maskDirs);

if ~input.doRandDir
    nTest = unique(input.testGratingContrastStepN).*unique(input.maskGratingContrastStepN);
else
    nTest = unique(input.testGratingDirectionStepN).*2;
end

data_test_avg = zeros(sz(1),sz(2),nTest);

% s = celleqel2mat_padded(input.tStimulusNumber);
% t = rem(s,8);
% maskPhase = 0 + (45.*t);
% maskPhases = unique(maskPhase);
% nMaskPhase = length(maskPhases);
maskPhase = celleqel2mat_padded(input.tTestMaskGratingPhaseDeg);
maskPhases = unique(maskPhase);
nMaskPhase = length(maskPhases);

ind_stimAlone = intersect(find(testCon),find(maskCon==0));
ind_maskAlone = intersect(find(testCon==0),find(maskCon));
ind_plaid = intersect(find(testCon),find(maskCon));

start = 1;
if ~input.doRandDir
    for it = 1:nTestCons
        ind_test = find(testCon == testCons(it));
        for im = 1:nMaskCons
            if im>1 & it>1
                for ip = 1:nMaskPhase
                    ind_mask = intersect(find(maskPhase == maskPhases(ip)),find(maskCon == maskCons(im)));
                    %ind = intersect(ind_noadapt, intersect(ind_test,ind_mask));
                    data_test_temp(:,:,1) = nanmean(data_adapt_dfof(:,:,intersect(ind_noadapt, intersect(ind_test,ind_mask))),3);
                    data_test_temp(:,:,2) = nanmean(data_adapt_dfof(:,:,intersect(ind_contadapt, intersect(ind_test,ind_mask))),3);
                    data_test_temp(:,:,3) = nanmean(data_adapt_dfof(:,:,intersect(ind_asynadapt, intersect(ind_test,ind_mask))),3);
                    data_test_temp(:,:,4) = nanmean(data_adapt_dfof(:,:,intersect(ind_singadapt, intersect(ind_test,ind_mask))),3);
                    data_test_avg(:,:,start) = max(data_test_temp,[],3);
                    start = start+1;
                end
            else
                ind_mask = find(maskCon == maskCons(im));
                %ind = intersect(ind_noadapt, intersect(ind_test,ind_mask));
                data_test_temp(:,:,1) = nanmean(data_adapt_dfof(:,:,intersect(ind_noadapt, intersect(ind_test,ind_mask))),3);
                data_test_temp(:,:,2) = nanmean(data_adapt_dfof(:,:,intersect(ind_contadapt, intersect(ind_test,ind_mask))),3);
                data_test_temp(:,:,3) = nanmean(data_adapt_dfof(:,:,intersect(ind_asynadapt, intersect(ind_test,ind_mask))),3);
                data_test_temp(:,:,4) = nanmean(data_adapt_dfof(:,:,intersect(ind_singadapt, intersect(ind_test,ind_mask))),3);
                data_test_avg(:,:,start) = max(data_test_temp,[],3);
                start = start+1;
            end
        end
    end
elseif input.doRandDir
    for iDir = 1:nTestDirs
        ind_stimdir = intersect(ind_noadapt,find(testDir == testDirs(iDir)));
        ind_maskdir = intersect(ind_noadapt,find(maskDir == maskDirs(iDir)));
        ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
        ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
        data_test_avg(:,:,start) = nanmean(data_adapt_dfof(:,:,ind_diralone),3);
        data_test_avg(:,:,start+1) = nanmean(data_adapt_dfof(:,:,ind_dirplaid),3);
        start = start+2;
    end
end

data_dfof = cat(3,data_adapt_avg,data_test_avg);
%data_dfof = data_test_avg;
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)
data_dfof = cat(3, data_dfof, data_dfof_max);
if doRedChannel
    data_dfof = cat(3,data_dfof,data_red_avg);
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cTest', 'testDir', 'maskDir', 'testDirs', 'maskDirs', 'nTestDirs', 'nMaskDirs', 'maskCon', 'testCon', 'testCons', 'maskCons', 'nMaskCons', 'nTestCons', 'maskPhase', 'maskPhases', 'nMaskPhase','frame_rate', 'nTestCons', 'nTrials', 'ind_noadapt','ind_asynadapt','ind_contadapt','ind_singadapt','data_dfof')

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)   
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    if doRedChannel & iStim==1
        red_mask = bwout;
    end
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
figure; imagesc(mask_cell)

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg
%% neuropil subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_cells')

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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%% Test stim analysis
prewin_frames = 15;
postwin_frames = 45;
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_test = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

for itrial = 1:nTrials
    if cTest(itrial) + postwin_frames < sz(3)
        data_test(:,:,itrial) = npSub_tc(cTest(itrial)-prewin_frames:cTest(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_test(1:prewin_frames,:,:),1);
data_dfof_tc = (data_test-data_f)./data_f;

nMaskCons = length(maskCons);
nTestCons = length(testCons);
noadapt_resp_square = zeros(nCells,nMaskCons,nTestCons);
noadapt_resp_cell = cell(nMaskCons,nTestCons);
trialsperstim = zeros(nMaskCons,nTestCons);
h_resp =zeros(nCells,nMaskCons,nTestCons);
p_resp =zeros(nCells,nMaskCons,nTestCons);
for im = 1:nMaskCons
    ind_mask = find(maskCon == maskCons(im));
    for it = 1:nTestCons
        ind_test = find(testCon == testCons(it));
        ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
        trialsperstim(im,it) = length(ind);
        noadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
        noadapt_resp_square(:,im,it) = mean(noadapt_resp_cell{im,it},2);
        [h_resp(:,im,it) p_resp(:,im,it)] = ttest2(noadapt_resp_cell{im,it},noadapt_resp_cell{1,1,1},'dim',2,'tail','right','alpha',0.05./(nMaskCons*nTestCons));
    end
end

% figure;
% [n,n2] = subplotn(nCells);
% for iCell = 1:nCells
%     subplot(n,n2,iCell)
%     indT = intersect(find(testCon == testCons(5)),find(maskCon == maskCons(1)));
%     indM = intersect(find(testCon == testCons(1)),find(maskCon == maskCons(5)));
%     indP = intersect(find(testCon == testCons(5)),find(maskCon == maskCons(5)));
%     plot(tt,mean(data_dfof_tc(:,iCell,indT),3))
%     hold on
%     plot(tt,mean(data_dfof_tc(:,iCell,indM),3))
%     plot(tt,mean(data_dfof_tc(:,iCell,indP),3))
%     %title(num2str(chop(squeeze(max(max(noadapt_resp_square(iCell,:,:),3),2)),2)))
% end

test_resp = noadapt_resp_cell{1,nTestCons};
mask_resp = noadapt_resp_cell{nMaskCons,1};
plaid_resp = noadapt_resp_cell{nMaskCons,nTestCons};
null_resp = noadapt_resp_cell{1,1};

[h_resptest p_resptest] =  ttest2(test_resp,null_resp,'dim',2,'tail','right');
[h_respmask p_respmask] = ttest2(mask_resp,null_resp,'dim',2,'tail','right');
[h_respplaid p_respplaid] = ttest2(plaid_resp,null_resp,'dim',2,'tail','right');
resptest_ind = find(h_resptest);
respmask_ind = find(h_respmask);
respplaid_ind = find(h_respplaid);

resptest_ind = setdiff(resptest_ind,red_cells);
respmask_ind = setdiff(respmask_ind,red_cells);
respplaid_ind = setdiff(respplaid_ind,red_cells);

resp_ind = unique([resptest_ind; respmask_ind; respplaid_ind]);

[h_preftest p_preftest] = ttest2(test_resp,mask_resp,'dim',2,'tail','right');
[h_prefmask p_prefmask] = ttest2(test_resp,mask_resp,'dim',2,'tail','left');
[h_prefplaid1 p_prefplaid1] = ttest2(plaid_resp, test_resp ,'dim',2,'tail','right');
[h_prefplaid2 p_prefplaid2] = ttest2(plaid_resp, mask_resp ,'dim',2,'tail','right');
preftest_ind = intersect(resp_ind,find(h_preftest));
prefmask_ind = intersect(resp_ind,find(h_prefmask));
prefplaid_ind = intersect(resp_ind,intersect(find(h_prefplaid1),find(h_prefplaid2)));
prefplaid_nonlinear = find(mean(mask_resp,2)+mean(test_resp,2)<mean(plaid_resp,2));
prefplaidonly_ind =  intersect(prefplaid_ind,prefplaid_nonlinear);
preftestonly_ind = setdiff(preftest_ind,prefplaidonly_ind);
prefmaskonly_ind = setdiff(prefmask_ind,prefplaidonly_ind);

preftestonly_ind = setdiff(preftestonly_ind,red_cells);
prefmaskonly_ind = setdiff(prefmaskonly_ind,red_cells);
prefplaidonly_ind = setdiff(prefplaidonly_ind,red_cells);

for i = 1:length(resp_ind)
    ii = resp_ind(i);
    figure;
    start = 1;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(data_dfof_tc(:,ii,ind),3))
            start = start+1;
        end
    end
    suptitle(['Cell ' num2str(ii)])
end

figure
[n,n2] = subplotn(length(resp_ind));
for i = 1:length(resp_ind)
    ii = resp_ind(i);
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
            temp(im,it) = mean(nanmean(data_dfof_tc(prewin_frames+6:prewin_frames+21,ii,ind),3),1);
        end
    end
    subplot(n,n2,i)
    imagesc(temp)
    title(['Cell ' num2str(ii)])
end

figure;
start = 1;
prefTest_noadapt_resp = zeros(nMaskCons,nTestCons,2);
prefMask_noadapt_resp = zeros(nMaskCons,nTestCons,2);
for im = 1:nMaskCons
    ind_mask = find(maskCon == maskCons(im));
    for it = 1:nTestCons
        ind_test = find(testCon == testCons(it));
        ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
        subplot(nMaskCons,nTestCons,start)
        plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind),2),3))
        hold on
        plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmask_ind,ind),2),3))
        plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind),2),3))
        prefTest_noadapt_resp(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(preftest_ind,:),1),2);
        prefTest_noadapt_resp(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
        prefMask_noadapt_resp(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefmask_ind,:),1),2);
        prefMask_noadapt_resp(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefmask_ind,:),2),[],1)./sqrt(length(prefmask_ind));
        title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
        ylim([-0.1 1])
        xlabel('Time (s)')
        start = start+1;
    end
end
suptitle(['No adapt- Blue Test Pref (n = ' num2str(length(preftest_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmask_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaid_ind)) ')'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noAdapt_TCs_groups.pdf']),'-dpdf','-bestfit');

if doRedChannel
    figure;
    start = 1;
    prefTest_noadapt_resp = zeros(nMaskCons,nTestCons,2);
    prefMask_noadapt_resp = zeros(nMaskCons,nTestCons,2);
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,setdiff(preftest_ind,red_cells),ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,setdiff(prefmask_ind,red_cells),ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,setdiff(prefplaid_ind,red_cells),ind),2),3))
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 .4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['No adapt- No red cells- Blue Test Pref (n = ' num2str(length(setdiff(preftest_ind,red_cells))) ') ; Red Mask Pref (n = ' num2str(length(setdiff(prefmask_ind,red_cells))) '); Yellow Plaid Pref (n = ' num2str(length(setdiff(prefplaid_ind,red_cells))) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noAdapt_TCs_groups_noRed.pdf']),'-dpdf','-bestfit');

    figure;
    start = 1;
    prefTest_noadapt_resp = zeros(nMaskCons,nTestCons,2);
    prefMask_noadapt_resp = zeros(nMaskCons,nTestCons,2);
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,intersect(preftest_ind,red_cells),ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,intersect(prefmask_ind,red_cells),ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,intersect(prefplaid_ind,red_cells),ind),2),3))
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 .4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['No adapt- Only red cells- Blue Test Pref (n = ' num2str(length(intersect(preftest_ind,red_cells))) ') ; Red Mask Pref (n = ' num2str(length(intersect(prefmask_ind,red_cells))) '); Yellow Plaid Pref (n = ' num2str(length(intersect(prefplaid_ind,red_cells))) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noAdapt_TCs_groups_redOnly.pdf']),'-dpdf','-bestfit');
end

%%
close all
start = 1;
contadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTest_contadapt_resp = nan(nMaskCons,nTestCons,2);
prefMask_contadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_contadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_contadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmask_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind),2),3))
            contadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTest_contadapt_resp(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(preftest_ind,:),1),2);
            prefTest_contadapt_resp(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
            prefMask_contadapt_resp(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefmask_ind,:),1),2);
            prefMask_contadapt_resp(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefmask_ind,:),2),[],1)./sqrt(length(prefmask_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Contingent adapt-  Blue Test Pref (n = ' num2str(length(preftest_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmask_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaid_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Contingent_TCs.pdf']),'-dpdf','-bestfit');
end


start = 1;
asynadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTest_asynadapt_resp = nan(nMaskCons,nTestCons,2);
prefMask_asynadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_asynadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_asynadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmask_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind),2),3))
            asynadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTest_asynadapt_resp(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(preftest_ind,:),1),2);
            prefTest_asynadapt_resp(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
            prefMask_asynadapt_resp(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefmask_ind,:),1),2);
            prefMask_asynadapt_resp(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefmask_ind,:),2),[],1)./sqrt(length(prefmask_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Asynchronous adapt-  Blue Test Pref (n = ' num2str(length(preftest_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmask_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaid_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Asynchronous_TCs.pdf']),'-dpdf','-bestfit');
end


start = 1;
singadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTest_singadapt_resp = nan(nMaskCons,nTestCons,2);
prefMask_singadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_singadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_singadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftest_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmask_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaid_ind,ind),2),3))
            singadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTest_singadapt_resp(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(preftest_ind,:),1),2);
            prefTest_singadapt_resp(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(preftest_ind,:),2),[],1)./sqrt(length(preftest_ind));
            prefMask_singadapt_resp(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(prefmask_ind,:),1),2);
            prefMask_singadapt_resp(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(prefmask_ind,:),2),[],1)./sqrt(length(prefmask_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Single adapt-  Blue Test Pref (n = ' num2str(length(preftest_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmask_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaid_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Single_TCs.pdf']),'-dpdf','-bestfit');
end

figure;
start = 1;
prefTestOnly_noadapt_resp = zeros(nMaskCons,nTestCons,2);
prefMaskOnly_noadapt_resp = zeros(nMaskCons,nTestCons,2);
prefPlaidOnly_noadapt_resp = zeros(nMaskCons,nTestCons,2);
for im = 1:nMaskCons
    ind_mask = find(maskCon == maskCons(im));
    for it = 1:nTestCons
        ind_test = find(testCon == testCons(it));
        ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
        subplot(nMaskCons,nTestCons,start)
        plot(tt,nanmean(nanmean(data_dfof_tc(:,preftestonly_ind,ind),2),3))
        hold on
        plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmaskonly_ind,ind),2),3))
        plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaidonly_ind,ind),2),3))
        prefTestOnly_noadapt_resp(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(preftestonly_ind,:),1),2);
        prefTestOnly_noadapt_resp(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(preftestonly_ind,:),2),[],1)./sqrt(length(preftestonly_ind));
        prefMaskOnly_noadapt_resp(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefmaskonly_ind,:),1),2);
        prefMaskOnly_noadapt_resp(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefmaskonly_ind,:),2),[],1)./sqrt(length(prefmaskonly_ind));
        prefPlaidOnly_noadapt_resp(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefplaidonly_ind,:),1),2);
        prefPlaidOnly_noadapt_resp(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefplaidonly_ind,:),2),[],1)./sqrt(length(prefplaidonly_ind));
        title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
        ylim([-0.1 0.4])
        xlabel('Time (s)')
        start = start+1;
    end
end
suptitle(['No adapt- Blue Test Pref (n = ' num2str(length(preftestonly_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmaskonly_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaidonly_ind)) ')'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noAdapt_TCs_onlygroups.pdf']),'-dpdf','-bestfit');

%%

start = 1;
contadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTestOnly_contadapt_resp = nan(nMaskCons,nTestCons,2);
prefMaskOnly_contadapt_resp = nan(nMaskCons,nTestCons,2);
prefPlaidOnly_contadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_contadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_contadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftestonly_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmaskonly_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaidonly_ind,ind),2),3))
            contadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTestOnly_contadapt_resp(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(preftestonly_ind,:),1),2);
            prefTestOnly_contadapt_resp(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(preftestonly_ind,:),2),[],1)./sqrt(length(preftestonly_ind));
            prefMaskOnly_contadapt_resp(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefmaskonly_ind,:),1),2);
            prefMaskOnly_contadapt_resp(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefmaskonly_ind,:),2),[],1)./sqrt(length(prefmaskonly_ind));
            prefPlaidOnly_contadapt_resp(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefplaidonly_ind,:),1),2);
            prefPlaidOnly_contadapt_resp(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefplaidonly_ind,:),2),[],1)./sqrt(length(prefplaidonly_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Contingent adapt-  Blue Test Pref (n = ' num2str(length(preftestonly_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmaskonly_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaidonly_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Contingent_TCs_only.pdf']),'-dpdf','-bestfit');
end


start = 1;
asynadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTestOnly_asynadapt_resp = nan(nMaskCons,nTestCons,2);
prefMaskOnly_asynadapt_resp = nan(nMaskCons,nTestCons,2);
prefPlaidOnly_asynadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_asynadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_asynadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftestonly_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmaskonly_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaidonly_ind,ind),2),3))
            asynadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTestOnly_asynadapt_resp(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(preftestonly_ind,:),1),2);
            prefTestOnly_asynadapt_resp(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(preftestonly_ind,:),2),[],1)./sqrt(length(preftestonly_ind));
            prefMaskOnly_asynadapt_resp(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefmaskonly_ind,:),1),2);
            prefMaskOnly_asynadapt_resp(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefmaskonly_ind,:),2),[],1)./sqrt(length(prefmaskonly_ind));
            prefPlaidOnly_asynadapt_resp(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefplaidonly_ind,:),1),2);
            prefPlaidOnly_asynadapt_resp(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefplaidonly_ind,:),2),[],1)./sqrt(length(prefplaidonly_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Asynchronous adapt-  Blue Test Pref (n = ' num2str(length(preftestonly_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmaskonly_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaidonly_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Asynchronous_TCs_only.pdf']),'-dpdf','-bestfit');
end

start = 1;
singadapt_resp_cell = cell(nMaskCons,nTestCons);
prefTestOnly_singadapt_resp = nan(nMaskCons,nTestCons,2);
prefMaskOnly_singadapt_resp = nan(nMaskCons,nTestCons,2);
prefPlaidOnly_singadapt_resp = nan(nMaskCons,nTestCons,2);
if length(ind_singadapt)>0
    figure;
    for im = 1:nMaskCons
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTestCons
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_singadapt,intersect(ind_test,ind_mask));
            subplot(nMaskCons,nTestCons,start)
            plot(tt,nanmean(nanmean(data_dfof_tc(:,preftestonly_ind,ind),2),3))
            hold on
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefmaskonly_ind,ind),2),3))
            plot(tt,nanmean(nanmean(data_dfof_tc(:,prefplaidonly_ind,ind),2),3))
            singadapt_resp_cell{im,it} = squeeze(mean(data_dfof_tc(prewin_frames+6:prewin_frames+21,:,ind),1));
            prefTestOnly_singadapt_resp(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(preftestonly_ind,:),1),2);
            prefTestOnly_singadapt_resp(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(preftestonly_ind,:),2),[],1)./sqrt(length(preftestonly_ind));
            prefMaskOnly_singadapt_resp(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(prefmaskonly_ind,:),1),2);
            prefMaskOnly_singadapt_resp(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(prefmaskonly_ind,:),2),[],1)./sqrt(length(prefmaskonly_ind));
            prefPlaidOnly_singadapt_resp(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(prefplaidonly_ind,:),1),2);
            prefPlaidOnly_singadapt_resp(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(prefplaidonly_ind,:),2),[],1)./sqrt(length(prefplaidonly_ind));
            title(['T: ' num2str(testCons(it)) '; M: ' num2str(maskCons(im))])
            ylim([-0.1 0.4])
            xlabel('Time (s)')
            start = start+1;
        end
    end
    suptitle(['Single adapt-  Blue Test Pref (n = ' num2str(length(preftestonly_ind)) ') ; Red Mask Pref (n = ' num2str(length(prefmaskonly_ind)) '); Yellow Plaid Pref (n = ' num2str(length(prefplaidonly_ind)) ')'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Single_TCs_only.pdf']),'-dpdf','-bestfit');
end


save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allCellResp.mat']),'noadapt_resp_cell','asynadapt_resp_cell','singadapt_resp_cell', 'contadapt_resp_cell','preftestonly_ind','prefmaskonly_ind','prefplaidonly_ind','resptest_ind','respmask_ind','respplaid_ind');


%%
figure;
subplot(2,2,1)
for im = 1:nMaskCons
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for im = 1:nMaskCons
    errorbar(testCons, prefTestOnly_contadapt_resp(im,:,1),prefTestOnly_contadapt_resp(im,:,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for im = 1:nMaskCons
    errorbar(testCons, prefTestOnly_asynadapt_resp(im,:,1),prefTestOnly_asynadapt_resp(im,:,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for im = 1:nMaskCons
    errorbar(testCons, prefTestOnly_singadapt_resp(im,:,1),prefTestOnly_singadapt_resp(im,:,2),'-o')
    hold on
end
title('Single')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(maskCons'))
suptitle(['Test preferring cells- n = ' num2str(length(preftestonly_ind))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_TestPrefOnly.pdf']),'-dpdf','-bestfit');


figure;
subplot(2,2,1)
for it = 1:nTestCons
    errorbar(maskCons, prefMaskOnly_noadapt_resp(:,it,1),prefMaskOnly_noadapt_resp(:,it,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for it = 1:nTestCons
    errorbar(maskCons, prefMaskOnly_contadapt_resp(:,it,1),prefMaskOnly_contadapt_resp(:,it,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for it = 1:nTestCons
    errorbar(maskCons, prefMaskOnly_asynadapt_resp(:,it,1),prefMaskOnly_asynadapt_resp(:,it,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for it = 1:nTestCons
    errorbar(maskCons, prefMaskOnly_singadapt_resp(:,it,1),prefMaskOnly_singadapt_resp(:,it,2),'-o')
    hold on
end
title('Single')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(testCons'))
suptitle(['Mask preferring cells- n = ' num2str(length(prefmaskonly_ind))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_MaskPrefOnly.pdf']),'-dpdf','-bestfit');

figure;
subplot(2,2,1)
for it = 1:nTestCons
    errorbar(maskCons, prefPlaidOnly_noadapt_resp(:,it,1),prefPlaidOnly_noadapt_resp(:,it,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for it = 1:nTestCons
    errorbar(maskCons, prefPlaidOnly_contadapt_resp(:,it,1),prefPlaidOnly_contadapt_resp(:,it,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for it = 1:nTestCons
    errorbar(maskCons, prefPlaidOnly_asynadapt_resp(:,it,1),prefPlaidOnly_asynadapt_resp(:,it,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for it = 1:nTestCons
    errorbar(maskCons, prefPlaidOnly_singadapt_resp(:,it,1),prefPlaidOnly_singadapt_resp(:,it,2),'-o')
    hold on
end
title('Single')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(testCons'))
suptitle(['Plaid preferring cells- n = ' num2str(length(prefplaidonly_ind))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_PlaidPrefOnly.pdf']),'-dpdf','-bestfit');

%%
figure;
subplot(1,2,1)
for im = 1:nMaskCons
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('Test Preferring')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.7])

subplot(1,2,2)
for it = 1:nTestCons
    errorbar(maskCons, prefMaskOnly_noadapt_resp(:,it,1),prefMaskOnly_noadapt_resp(:,it,2),'-o')
    hold on
end
title('Mask Preferring')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 .7])

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_NoAdapt.pdf']),'-dpdf','-bestfit');

%%

% T_noadapt = [nanmean(noadapt_resp_cell{1,4}(preftest_ind,:),2); nanmean(noadapt_resp_cell{2,1}(prefmask_ind,:),2)];
% TM1_noadapt = [nanmean(noadapt_resp_cell{2,4}(preftest_ind,:),2); nanmean(noadapt_resp_cell{2,4}(prefmask_ind,:),2)];
% TM2_noadapt = [nanmean(noadapt_resp_cell{3,4}(preftest_ind,:),2); nanmean(noadapt_resp_cell{2,5}(prefmask_ind,:),2)];
% T_noadapt(find(T_noadapt<0)) = 0;
% TM1_noadapt(find(TM1_noadapt<0)) = 0;
% TM2_noadapt(find(TM2_noadapt<0)) = 0;
% MI1_noadapt = (T_noadapt-TM1_noadapt)./(T_noadapt+TM1_noadapt);
% MI2_noadapt = (T_noadapt-TM2_noadapt)./(T_noadapt+TM2_noadapt);
% 
% noTM_asynadapt = [nanmean(asynadapt_resp_cell{1,1}(preftest_ind,:),2); nanmean(asynadapt_resp_cell{1,1}(prefmask_ind,:),2)];
% T_asynadapt = [nanmean(asynadapt_resp_cell{1,4}(preftest_ind,:),2); nanmean(asynadapt_resp_cell{2,1}(prefmask_ind,:),2)]-noTM_asynadapt;
% TM1_asynadapt = [nanmean(asynadapt_resp_cell{2,4}(preftest_ind,:),2); nanmean(asynadapt_resp_cell{2,4}(prefmask_ind,:),2)]-noTM_asynadapt;
% TM2_asynadapt = [nanmean(asynadapt_resp_cell{3,4}(preftest_ind,:),2); nanmean(asynadapt_resp_cell{2,5}(prefmask_ind,:),2)]-noTM_asynadapt;
% T_asynadapt(find(T_asynadapt<0)) = 0;
% TM1_asynadapt(find(TM1_asynadapt<0)) = 0;
% TM2_asynadapt(find(TM2_asynadapt<0)) = 0;
% MI1_asynadapt = (T_asynadapt-TM1_asynadapt)./(T_asynadapt+TM1_asynadapt);
% MI2_asynadapt = (T_asynadapt-TM2_asynadapt)./(T_asynadapt+TM2_asynadapt);
% 
% noTM_contadapt = [nanmean(contadapt_resp_cell{1,1}(preftest_ind,:),2); nanmean(contadapt_resp_cell{1,1}(prefmask_ind,:),2)];
% T_contadapt = [nanmean(contadapt_resp_cell{1,4}(preftest_ind,:),2); nanmean(contadapt_resp_cell{2,1}(prefmask_ind,:),2)]-noTM_contadapt;
% TM1_contadapt = [nanmean(contadapt_resp_cell{2,4}(preftest_ind,:),2); nanmean(contadapt_resp_cell{2,4}(prefmask_ind,:),2)]-noTM_contadapt;
% TM2_contadapt = [nanmean(contadapt_resp_cell{3,4}(preftest_ind,:),2); nanmean(contadapt_resp_cell{2,5}(prefmask_ind,:),2)]-noTM_contadapt;
% T_contadapt(find(T_contadapt<0)) = 0;
% TM1_contadapt(find(TM1_contadapt<0)) = 0;
% TM2_contadapt(find(TM2_contadapt<0)) = 0;
% MI1_contadapt = (T_contadapt-TM1_contadapt)./(T_contadapt+TM1_contadapt);
% MI2_contadapt = (T_contadapt-TM2_contadapt)./(T_contadapt+TM2_contadapt);
% 
% figure; 
% subplot(1,2,1)
% errorbar([0.25 0.5], [nanmean(MI1_noadapt) nanmean(MI2_noadapt)], [nanstd(MI1_noadapt)./sqrt(size(MI1_noadapt,1)) nanstd(MI2_noadapt)./sqrt(size(MI2_noadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt) nanmean(MI2_asynadapt)], [nanstd(MI1_asynadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt) nanmean(MI2_contadapt)], [nanstd(MI1_contadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% legend('No adapt','Asynchronous','Contingent');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index')
% 
% subplot(1,2,2)
% plot(nan,nan)
% hold on 
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt-MI1_noadapt) nanmean(MI2_asynadapt-MI2_noadapt)], [nanstd(MI1_asynadapt-MI1_noadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt-MI2_noadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt-MI1_noadapt) nanmean(MI2_contadapt-MI2_noadapt)], [nanstd(MI1_contadapt-MI1_noadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt-MI2_noadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% legend('','Asynchronous','Contingent');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index (Post-Pre)')
% suptitle(['All cells- n = ' num2str(size(MI1_contadapt,1)) ' cells'])
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ModulationSummary.pdf']),'-dpdf','-bestfit');
% 
% T_noadapt = nanmean(noadapt_resp_cell{1,4}(preftestonly_ind,:),2);
% TM1_noadapt = nanmean(noadapt_resp_cell{2,4}(preftestonly_ind,:),2);
% TM2_noadapt = nanmean(noadapt_resp_cell{3,4}(preftestonly_ind,:),2);
% T_noadapt(find(T_noadapt<0)) = 0;
% TM1_noadapt(find(TM1_noadapt<0)) = 0;
% TM2_noadapt(find(TM2_noadapt<0)) = 0;
% MI1_noadapt = (T_noadapt-TM1_noadapt)./(T_noadapt+TM1_noadapt);
% MI2_noadapt = (T_noadapt-TM2_noadapt)./(T_noadapt+TM2_noadapt);
% 
% noTM_asynadapt = nanmean(asynadapt_resp_cell{1,1}(preftestonly_ind,:),2);
% T_asynadapt = nanmean(asynadapt_resp_cell{1,4}(preftestonly_ind,:),2)-noTM_asynadapt;
% TM1_asynadapt = nanmean(asynadapt_resp_cell{2,4}(preftestonly_ind,:),2)-noTM_asynadapt;
% TM2_asynadapt = nanmean(asynadapt_resp_cell{3,4}(preftestonly_ind,:),2)-noTM_asynadapt;
% T_asynadapt(find(T_asynadapt<0)) = 0;
% TM1_asynadapt(find(TM1_asynadapt<0)) = 0;
% TM2_asynadapt(find(TM2_asynadapt<0)) = 0;
% MI1_asynadapt = (T_asynadapt-TM1_asynadapt)./(T_asynadapt+TM1_asynadapt);
% MI2_asynadapt = (T_asynadapt-TM2_asynadapt)./(T_asynadapt+TM2_asynadapt);
% 
% noTM_contadapt = nanmean(contadapt_resp_cell{1,1}(preftestonly_ind,:),2);
% T_contadapt = nanmean(contadapt_resp_cell{1,4}(preftestonly_ind,:),2)-noTM_contadapt;
% TM1_contadapt = nanmean(contadapt_resp_cell{2,4}(preftestonly_ind,:),2)-noTM_contadapt;
% TM2_contadapt = nanmean(contadapt_resp_cell{3,4}(preftestonly_ind,:),2)-noTM_contadapt;
% T_contadapt(find(T_contadapt<0)) = 0;
% TM1_contadapt(find(TM1_contadapt<0)) = 0;
% TM2_contadapt(find(TM2_contadapt<0)) = 0;
% MI1_contadapt = (T_contadapt-TM1_contadapt)./(T_contadapt+TM1_contadapt);
% MI2_contadapt = (T_contadapt-TM2_contadapt)./(T_contadapt+TM2_contadapt);
% 
% figure; 
% subplot(1,2,1)
% errorbar([0.25 0.5], [nanmean(MI1_noadapt) nanmean(MI2_noadapt)], [nanstd(MI1_noadapt)./sqrt(size(MI1_noadapt,1)) nanstd(MI2_noadapt)./sqrt(size(MI2_noadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt) nanmean(MI2_asynadapt)], [nanstd(MI1_asynadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt) nanmean(MI2_contadapt)], [nanstd(MI1_contadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% legend('No adapt','Asynchronous','Contingent','location','southeast');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index')
% 
% subplot(1,2,2)
% plot(nan,nan)
% hold on 
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt-MI1_noadapt) nanmean(MI2_asynadapt-MI2_noadapt)], [nanstd(MI1_asynadapt-MI1_noadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt-MI2_noadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt-MI1_noadapt) nanmean(MI2_contadapt-MI2_noadapt)], [nanstd(MI1_contadapt-MI1_noadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt-MI2_noadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index (Post-Pre)')
% suptitle(['Test Pref only cells- n = ' num2str(size(MI1_contadapt,1)) ' cells'])
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ModulationSummary_TestPrefOnly.pdf']),'-dpdf','-bestfit');


[h_no_testresp p_no_testresp] = ttest2(noadapt_resp_cell{1,nTestCons},noadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_no_maskresp p_no_maskresp] = ttest2(noadapt_resp_cell{nMaskCons,1},noadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_no_plaidresp p_no_plaidresp] = ttest2(noadapt_resp_cell{nMaskCons,nTestCons},noadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_asyn_testresp p_asyn_testresp] = ttest2(asynadapt_resp_cell{1,nTestCons},asynadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_asyn_maskresp p_asyn_maskresp] = ttest2(asynadapt_resp_cell{nMaskCons,1},asynadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_asyn_plaidresp p_asyn_plaidresp] = ttest2(asynadapt_resp_cell{nMaskCons,nTestCons},asynadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_sing_testresp p_sing_testresp] = ttest2(singadapt_resp_cell{1,nTestCons},singadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_sing_maskresp p_sing_maskresp] = ttest2(singadapt_resp_cell{nMaskCons,1},singadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_sing_plaidresp p_sing_plaidresp] = ttest2(singadapt_resp_cell{nMaskCons,nTestCons},singadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_cont_testresp p_cont_testresp] = ttest2(contadapt_resp_cell{1,nTestCons},contadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_cont_maskresp p_cont_maskresp] = ttest2(contadapt_resp_cell{nMaskCons,1},contadapt_resp_cell{1,1},'dim',2,'tail','right');
[h_cont_plaidresp p_cont_plaidresp] = ttest2(contadapt_resp_cell{nMaskCons,nTestCons},contadapt_resp_cell{1,1},'dim',2,'tail','right');

preftest_noadaptresp_ind = intersect(preftest_ind, find(h_no_testresp));
prefmask_noadaptresp_ind = intersect(prefmask_ind, find(h_no_maskresp));
prefplaid_noadaptresp_ind = intersect(prefplaid_ind, find(h_no_plaidresp));
preftest_asynresp_ind = intersect(preftest_ind, find(h_asyn_testresp));
prefmask_asynresp_ind = intersect(prefmask_ind, find(h_asyn_maskresp));
prefplaid_asynresp_ind = intersect(prefplaid_ind, find(h_asyn_plaidresp));
preftest_singresp_ind = intersect(preftest_ind, find(h_sing_testresp));
prefmask_singresp_ind = intersect(prefmask_ind, find(h_sing_maskresp));
prefplaid_singresp_ind = intersect(prefplaid_ind, find(h_sing_plaidresp));
preftest_contresp_ind = intersect(preftest_ind, find(h_cont_testresp));
prefmask_contresp_ind = intersect(prefmask_ind, find(h_cont_maskresp));
prefplaid_contresp_ind = intersect(prefplaid_ind, find(h_cont_plaidresp));

prefTest_noadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefMask_noadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefPlaid_noadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefTest_asynadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefMask_asynadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefPlaid_asynadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefTest_singadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefMask_singadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefPlaid_singadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefTest_contadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefMask_contadapt_respcells = zeros(nMaskCons,nTestCons,2);
prefPlaid_contadapt_respcells = zeros(nMaskCons,nTestCons,2);

for im = 1:nMaskCons
    for it = 1:nTestCons
        prefTest_noadapt_respcells(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind,:),1),2);
        prefTest_noadapt_respcells(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind,:),2),[],1)./sqrt(length(preftest_noadaptresp_ind));
        prefMask_noadapt_respcells(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind,:),1),2);
        prefMask_noadapt_respcells(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind,:),2),[],1)./sqrt(length(prefmask_noadaptresp_ind));
        prefPlaid_noadapt_respcells(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefplaid_noadaptresp_ind,:),1),2);
        prefPlaid_noadapt_respcells(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefplaid_noadaptresp_ind,:),2),[],1)./sqrt(length(prefplaid_noadaptresp_ind));
        prefTest_asynadapt_respcells(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind,:),1),2);
        prefTest_asynadapt_respcells(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind,:),2),[],1)./sqrt(length(preftest_asynresp_ind));
        prefMask_asynadapt_respcells(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind,:),1),2);
        prefMask_asynadapt_respcells(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind,:),2),[],1)./sqrt(length(prefmask_asynresp_ind));
        prefPlaid_asynadapt_respcells(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefplaid_asynresp_ind,:),1),2);
        prefPlaid_asynadapt_respcells(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefplaid_asynresp_ind,:),2),[],1)./sqrt(length(prefplaid_asynresp_ind));
        prefTest_singadapt_respcells(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(preftest_singresp_ind,:),1),2);
        prefTest_singadapt_respcells(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(preftest_singresp_ind,:),2),[],1)./sqrt(length(preftest_singresp_ind));
        prefMask_singadapt_respcells(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(prefmask_singresp_ind,:),1),2);
        prefMask_singadapt_respcells(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(prefmask_singresp_ind,:),2),[],1)./sqrt(length(prefmask_singresp_ind));
        prefPlaid_singadapt_respcells(im,it,1) = nanmean(nanmean(singadapt_resp_cell{im,it}(prefplaid_singresp_ind,:),1),2);
        prefPlaid_singadapt_respcells(im,it,2) = nanstd(nanmean(singadapt_resp_cell{im,it}(prefplaid_singresp_ind,:),2),[],1)./sqrt(length(prefplaid_singresp_ind));
        prefTest_contadapt_respcells(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind,:),1),2);
        prefTest_contadapt_respcells(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind,:),2),[],1)./sqrt(length(preftest_contresp_ind));
        prefMask_contadapt_respcells(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind,:),1),2);
        prefMask_contadapt_respcells(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind,:),2),[],1)./sqrt(length(prefmask_contresp_ind));
        prefPlaid_contadapt_respcells(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefplaid_contresp_ind,:),1),2);
        prefPlaid_contadapt_respcells(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefplaid_contresp_ind,:),2),[],1)./sqrt(length(prefplaid_contresp_ind));
    end
end

figure;
subplot(2,2,1)
for im = 1:nMaskCons
    errorbar(testCons, prefTest_noadapt_respcells(im,:,1),prefTest_noadapt_respcells(im,:,2),'-o')
    hold on
end
title(['No adapt- n = ' num2str(length(preftest_noadaptresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for im = 1:nMaskCons
    errorbar(testCons, prefTest_contadapt_respcells(im,:,1),prefTest_contadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Contingent- n = ' num2str(length(preftest_contresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for im = 1:nMaskCons
    errorbar(testCons, prefTest_asynadapt_respcells(im,:,1),prefTest_asynadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Asynchronous- n = ' num2str(length(preftest_asynresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for im = 1:nMaskCons
    errorbar(testCons, prefTest_singadapt_respcells(im,:,1),prefTest_singadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Single- n = ' num2str(length(preftest_singresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(maskCons'))
suptitle(['Test preferring (n = ' num2str(length(preftest_ind)) ') and responsive cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_TestPref_RespCells.pdf']),'-dpdf','-bestfit');


figure;
subplot(2,2,1)
for it = 1:nTestCons
    errorbar(maskCons, prefMask_noadapt_respcells(:,it,1),prefMask_noadapt_respcells(:,it,2),'-o')
    hold on
end
title(['No adapt- n = ' num2str(length(prefmask_noadaptresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for it = 1:nTestCons
    errorbar(maskCons, prefMask_contadapt_respcells(:,it,1),prefMask_contadapt_respcells(:,it,2),'-o')
    hold on
end
title(['Contingent- n = ' num2str(length(prefmask_contresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for it = 1:nTestCons
    errorbar(maskCons, prefMask_asynadapt_respcells(:,it,1),prefMask_asynadapt_respcells(:,it,2),'-o')
    hold on
end
title(['Asynchronous- n = ' num2str(length(prefmask_asynresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for it = 1:nTestCons
    errorbar(maskCons, prefMask_singadapt_respcells(:,it,1),prefMask_singadapt_respcells(:,it,2),'-o')
    hold on
end
title(['Single- n = ' num2str(length(prefmask_singresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(testCons'))
suptitle(['Mask preferring (n = ' num2str(length(prefmask_ind)) ') and responsive cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_MaskPref_RespCells.pdf']),'-dpdf','-bestfit');

figure;
subplot(2,2,1)
for im = 1:nMaskCons
    errorbar(testCons, prefPlaid_noadapt_respcells(im,:,1),prefPlaid_noadapt_respcells(im,:,2),'-o')
    hold on
end
title(['No adapt- n = ' num2str(length(prefplaid_noadaptresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,2)
for im = 1:nMaskCons
    errorbar(testCons, prefPlaid_contadapt_respcells(im,:,1),prefPlaid_contadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Contingent- n = ' num2str(length(prefplaid_contresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,3)
for im = 1:nMaskCons
    errorbar(testCons, prefPlaid_asynadapt_respcells(im,:,1),prefPlaid_asynadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Asynchronous- n = ' num2str(length(prefplaid_asynresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
subplot(2,2,4)
for im = 1:nMaskCons
    errorbar(testCons, prefPlaid_singadapt_respcells(im,:,1),prefPlaid_singadapt_respcells(im,:,2),'-o')
    hold on
end
title(['Single- n = ' num2str(length(prefplaid_singresp_ind))])
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 1])
legend(num2str(maskCons'))
suptitle(['Plaid preferring (n = ' num2str(length(prefplaid_ind)) ') and responsive cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_PlaidPref_RespCells.pdf']),'-dpdf','-bestfit');


if doRedChannel
    preftest_noadaptresp_ind_red = intersect(intersect(preftest_ind, find(h_no_testresp)),red_cells);
    prefmask_noadaptresp_ind_red = intersect(intersect(prefmask_ind, find(h_no_maskresp)),red_cells);
    preftest_asynresp_ind_red = intersect(intersect(preftest_ind, find(h_asyn_testresp)),red_cells);
    prefmask_asynresp_ind_red = intersect(intersect(prefmask_ind, find(h_asyn_maskresp)),red_cells);
    preftest_contresp_ind_red = intersect(intersect(preftest_ind, find(h_cont_testresp)),red_cells);
    prefmask_contresp_ind_red = intersect(intersect(prefmask_ind, find(h_cont_maskresp)),red_cells);
    
    preftest_noadaptresp_ind_notred = setdiff(intersect(preftest_ind, find(h_no_testresp)),red_cells);
    prefmask_noadaptresp_ind_notred = setdiff(intersect(prefmask_ind, find(h_no_maskresp)),red_cells);
    preftest_asynresp_ind_notred = setdiff(intersect(preftest_ind, find(h_asyn_testresp)),red_cells);
    prefmask_asynresp_ind_notred = setdiff(intersect(prefmask_ind, find(h_asyn_maskresp)),red_cells);
    preftest_contresp_ind_notred = setdiff(intersect(preftest_ind, find(h_cont_testresp)),red_cells);
    prefmask_contresp_ind_notred = setdiff(intersect(prefmask_ind, find(h_cont_maskresp)),red_cells);

    prefTest_noadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    prefMask_noadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    prefTest_asynadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    prefMask_asynadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    prefTest_contadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    prefMask_contadapt_respcells_red = zeros(nMaskCons,nTestCons,2);
    
    prefTest_noadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    prefMask_noadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    prefTest_asynadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    prefMask_asynadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    prefTest_contadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    prefMask_contadapt_respcells_notred = zeros(nMaskCons,nTestCons,2);
    for im = 1:nMaskCons
        for it = 1:nTestCons
            prefTest_noadapt_respcells_red(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind_red,:),1),2);
            prefTest_noadapt_respcells_red(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind_red,:),2),[],1)./sqrt(length(preftest_noadaptresp_ind_red));
            prefMask_noadapt_respcells_red(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind_red,:),1),2);
            prefMask_noadapt_respcells_red(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind_red,:),2),[],1)./sqrt(length(prefmask_noadaptresp_ind_red));
            prefTest_asynadapt_respcells_red(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind_red,:),1),2);
            prefTest_asynadapt_respcells_red(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind_red,:),2),[],1)./sqrt(length(preftest_asynresp_ind_red));
            prefMask_asynadapt_respcells_red(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind_red,:),1),2);
            prefMask_asynadapt_respcells_red(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind_red,:),2),[],1)./sqrt(length(prefmask_asynresp_ind_red));
            prefTest_contadapt_respcells_red(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind_red,:),1),2);
            prefTest_contadapt_respcells_red(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind_red,:),2),[],1)./sqrt(length(preftest_contresp_ind_red));
            prefMask_contadapt_respcells_red(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind_red,:),1),2);
            prefMask_contadapt_respcells_red(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind_red,:),2),[],1)./sqrt(length(prefmask_contresp_ind_red));
            
            prefTest_noadapt_respcells_notred(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind_notred,:),1),2);
            prefTest_noadapt_respcells_notred(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(preftest_noadaptresp_ind_notred,:),2),[],1)./sqrt(length(preftest_noadaptresp_ind_notred));
            prefMask_noadapt_respcells_notred(im,it,1) = nanmean(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind_notred,:),1),2);
            prefMask_noadapt_respcells_notred(im,it,2) = nanstd(nanmean(noadapt_resp_cell{im,it}(prefmask_noadaptresp_ind_notred,:),2),[],1)./sqrt(length(prefmask_noadaptresp_ind_notred));
            prefTest_asynadapt_respcells_notred(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind_notred,:),1),2);
            prefTest_asynadapt_respcells_notred(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(preftest_asynresp_ind_notred,:),2),[],1)./sqrt(length(preftest_asynresp_ind_notred));
            prefMask_asynadapt_respcells_notred(im,it,1) = nanmean(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind_notred,:),1),2);
            prefMask_asynadapt_respcells_notred(im,it,2) = nanstd(nanmean(asynadapt_resp_cell{im,it}(prefmask_asynresp_ind_notred,:),2),[],1)./sqrt(length(prefmask_asynresp_ind_notred));
            prefTest_contadapt_respcells_notred(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind_notred,:),1),2);
            prefTest_contadapt_respcells_notred(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(preftest_contresp_ind_notred,:),2),[],1)./sqrt(length(preftest_contresp_ind_notred));
            prefMask_contadapt_respcells_notred(im,it,1) = nanmean(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind_notred,:),1),2);
            prefMask_contadapt_respcells_notred(im,it,2) = nanstd(nanmean(contadapt_resp_cell{im,it}(prefmask_contresp_ind_notred,:),2),[],1)./sqrt(length(prefmask_contresp_ind_notred));
        end
    end

    figure;
    subplot(1,3,1)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_noadapt_respcells_red(im,:,1),prefTest_noadapt_respcells_red(im,:,2),'-o')
        hold on
    end
    title(['No adapt- n = ' num2str(length(preftest_noadaptresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,2)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_contadapt_respcells_red(im,:,1),prefTest_contadapt_respcells_red(im,:,2),'-o')
        hold on
    end
    title(['Contingent- n = ' num2str(length(preftest_contresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,3)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_asynadapt_respcells_red(im,:,1),prefTest_asynadapt_respcells_red(im,:,2),'-o')
        hold on
    end
    title(['Asynchronous- n = ' num2str(length(preftest_asynresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    legend(num2str(maskCons'))
    suptitle(['Test preferring (n = ' num2str(length(intersect(red_cells,preftest_ind))) ') and responsive red cells'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_TestPref_RespCells_Red.pdf']),'-dpdf','-bestfit');


    figure;
    subplot(1,3,1)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_noadapt_respcells_red(:,it,1),prefMask_noadapt_respcells_red(:,it,2),'-o')
        hold on
    end
    title(['No adapt- n = ' num2str(length(prefmask_noadaptresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,2)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_contadapt_respcells_red(:,it,1),prefMask_contadapt_respcells_red(:,it,2),'-o')
        hold on
    end
    title(['Contingent- n = ' num2str(length(prefmask_contresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,3)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_asynadapt_respcells_red(:,it,1),prefMask_asynadapt_respcells_red(:,it,2),'-o')
        hold on
    end
    title(['Asynchronous- n = ' num2str(length(prefmask_asynresp_ind_red))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    legend(num2str(testCons'))
    suptitle(['Mask preferring (n = ' num2str(length(intersect(red_cells,prefmask_ind))) ') and responsive red cells'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_MaskPref_RespCells_Red.pdf']),'-dpdf','-bestfit');

    figure;
    subplot(1,3,1)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_noadapt_respcells_notred(im,:,1),prefTest_noadapt_respcells_notred(im,:,2),'-o')
        hold on
    end
    title(['No adapt- n = ' num2str(length(preftest_noadaptresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,2)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_contadapt_respcells_notred(im,:,1),prefTest_contadapt_respcells_notred(im,:,2),'-o')
        hold on
    end
    title(['Contingent- n = ' num2str(length(preftest_contresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,3)
    for im = 1:nMaskCons
        errorbar(testCons, prefTest_asynadapt_respcells_notred(im,:,1),prefTest_asynadapt_respcells_notred(im,:,2),'-o')
        hold on
    end
    title(['Asynchronous- n = ' num2str(length(preftest_asynresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    legend(num2str(maskCons'))
    suptitle(['Test preferring (n = ' num2str(length(setdiff(preftest_ind,red_cells))) ') and responsive not red cells'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_TestPref_RespCells_NotRed.pdf']),'-dpdf','-bestfit');


    figure;
    subplot(1,3,1)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_noadapt_respcells_notred(:,it,1),prefMask_noadapt_respcells_notred(:,it,2),'-o')
        hold on
    end
    title(['No adapt- n = ' num2str(length(prefmask_noadaptresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,2)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_contadapt_respcells_notred(:,it,1),prefMask_contadapt_respcells_notred(:,it,2),'-o')
        hold on
    end
    title(['Contingent- n = ' num2str(length(prefmask_contresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    subplot(1,3,3)
    for it = 1:nTestCons
        errorbar(maskCons, prefMask_asynadapt_respcells_notred(:,it,1),prefMask_asynadapt_respcells_notred(:,it,2),'-o')
        hold on
    end
    title(['Asynchronous- n = ' num2str(length(prefmask_asynresp_ind_notred))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
    legend(num2str(testCons'))
    suptitle(['Mask preferring (n = ' num2str(length(setdiff(prefmask_ind,red_cells))) ') and responsive not red cells'])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allContrastResponse_MaskPref_RespCells_NotRed.pdf']),'-dpdf','-bestfit');

end
% T_noadapt = [nanmean(noadapt_resp_cell{1,4}(preftest_noadaptresp_ind,:),2); nanmean(noadapt_resp_cell{2,1}(prefmask_noadaptresp_ind,:),2)];
% TM1_noadapt = [nanmean(noadapt_resp_cell{2,4}(preftest_noadaptresp_ind,:),2); nanmean(noadapt_resp_cell{2,4}(prefmask_noadaptresp_ind,:),2)];
% TM2_noadapt = [nanmean(noadapt_resp_cell{3,4}(preftest_noadaptresp_ind,:),2); nanmean(noadapt_resp_cell{2,5}(prefmask_noadaptresp_ind,:),2)];
% T_noadapt(find(T_noadapt<0)) = 0;
% TM1_noadapt(find(TM1_noadapt<0)) = 0;
% TM2_noadapt(find(TM2_noadapt<0)) = 0;
% MI1_noadapt = (T_noadapt-TM1_noadapt)./(T_noadapt+TM1_noadapt);
% MI2_noadapt = (T_noadapt-TM2_noadapt)./(T_noadapt+TM2_noadapt);
% 
% noTM_asynadapt = [nanmean(asynadapt_resp_cell{1,1}(preftest_asynresp_ind,:),2); nanmean(asynadapt_resp_cell{1,1}(prefmask_asynresp_ind,:),2)];
% T_asynadapt = [nanmean(asynadapt_resp_cell{1,4}(preftest_asynresp_ind,:),2); nanmean(asynadapt_resp_cell{2,1}(prefmask_asynresp_ind,:),2)]-noTM_asynadapt;
% TM1_asynadapt = [nanmean(asynadapt_resp_cell{2,4}(preftest_asynresp_ind,:),2); nanmean(asynadapt_resp_cell{2,4}(prefmask_asynresp_ind,:),2)]-noTM_asynadapt;
% TM2_asynadapt = [nanmean(asynadapt_resp_cell{3,4}(preftest_asynresp_ind,:),2); nanmean(asynadapt_resp_cell{2,5}(prefmask_asynresp_ind,:),2)]-noTM_asynadapt;
% T_asynadapt(find(T_asynadapt<0)) = 0;
% TM1_asynadapt(find(TM1_asynadapt<0)) = 0;
% TM2_asynadapt(find(TM2_asynadapt<0)) = 0;
% MI1_asynadapt = (T_asynadapt-TM1_asynadapt)./(T_asynadapt+TM1_asynadapt);
% MI2_asynadapt = (T_asynadapt-TM2_asynadapt)./(T_asynadapt+TM2_asynadapt);
% 
% noTM_contadapt = [nanmean(contadapt_resp_cell{1,1}(preftest_contresp_ind,:),2); nanmean(contadapt_resp_cell{1,1}(prefmask_contresp_ind,:),2)];
% T_contadapt = [nanmean(contadapt_resp_cell{1,4}(preftest_contresp_ind,:),2); nanmean(contadapt_resp_cell{2,1}(prefmask_contresp_ind,:),2)]-noTM_contadapt;
% TM1_contadapt = [nanmean(contadapt_resp_cell{2,4}(preftest_contresp_ind,:),2); nanmean(contadapt_resp_cell{2,4}(prefmask_contresp_ind,:),2)]-noTM_contadapt;
% TM2_contadapt = [nanmean(contadapt_resp_cell{3,4}(preftest_contresp_ind,:),2); nanmean(contadapt_resp_cell{2,5}(prefmask_contresp_ind,:),2)]-noTM_contadapt;
% T_contadapt(find(T_contadapt<0)) = 0;
% TM1_contadapt(find(TM1_contadapt<0)) = 0;
% TM2_contadapt(find(TM2_contadapt<0)) = 0;
% MI1_contadapt = (T_contadapt-TM1_contadapt)./(T_contadapt+TM1_contadapt);
% MI2_contadapt = (T_contadapt-TM2_contadapt)./(T_contadapt+TM2_contadapt);
% 
% T_noadapt_asyn = [nanmean(noadapt_resp_cell{1,4}(preftest_asynresp_ind,:),2); nanmean(noadapt_resp_cell{2,1}(prefmask_asynresp_ind,:),2)];
% TM1_noadapt_asyn = [nanmean(noadapt_resp_cell{2,4}(preftest_asynresp_ind,:),2); nanmean(noadapt_resp_cell{2,4}(prefmask_asynresp_ind,:),2)];
% TM2_noadapt_asyn = [nanmean(noadapt_resp_cell{3,4}(preftest_asynresp_ind,:),2); nanmean(noadapt_resp_cell{2,5}(prefmask_asynresp_ind,:),2)];
% T_noadapt_asyn(find(T_noadapt_asyn<0)) = 0;
% TM1_noadapt_asyn(find(TM1_noadapt_asyn<0)) = 0;
% TM2_noadapt_asyn(find(TM2_noadapt_asyn<0)) = 0;
% MI1_noadapt_asyn = (T_noadapt_asyn-TM1_noadapt_asyn)./(T_noadapt_asyn+TM1_noadapt_asyn);
% MI2_noadapt_asyn = (T_noadapt_asyn-TM2_noadapt_asyn)./(T_noadapt_asyn+TM2_noadapt_asyn);
% 
% T_noadapt_cont = [nanmean(noadapt_resp_cell{1,4}(preftest_contresp_ind,:),2); nanmean(noadapt_resp_cell{2,1}(prefmask_contresp_ind,:),2)];
% TM1_noadapt_cont = [nanmean(noadapt_resp_cell{2,4}(preftest_contresp_ind,:),2); nanmean(noadapt_resp_cell{2,4}(prefmask_contresp_ind,:),2)];
% TM2_noadapt_cont = [nanmean(noadapt_resp_cell{3,4}(preftest_contresp_ind,:),2); nanmean(noadapt_resp_cell{2,5}(prefmask_contresp_ind,:),2)];
% T_noadapt_cont(find(T_noadapt_cont<0)) = 0;
% TM1_noadapt_cont(find(TM1_noadapt_cont<0)) = 0;
% TM2_noadapt_cont(find(TM2_noadapt_cont<0)) = 0;
% MI1_noadapt_cont = (T_noadapt_cont-TM1_noadapt_cont)./(T_noadapt_cont+TM1_noadapt_cont);
% MI2_noadapt_cont = (T_noadapt_cont-TM2_noadapt_cont)./(T_noadapt_cont+TM2_noadapt_cont);
% 
% figure; 
% subplot(1,2,1)
% errorbar([0.25 0.5], [nanmean(MI1_noadapt) nanmean(MI2_noadapt)], [nanstd(MI1_noadapt)./sqrt(size(MI1_noadapt,1)) nanstd(MI2_noadapt)./sqrt(size(MI2_noadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt) nanmean(MI2_asynadapt)], [nanstd(MI1_asynadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt) nanmean(MI2_contadapt)], [nanstd(MI1_contadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% legend(['No adapt- n = ' num2str(size(MI1_noadapt,1))],['Asynchronous- n = ' num2str(size(MI1_asynadapt,1))],['Contingent - n = ' num2str(size(MI1_contadapt,1))]);
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index')
% 
% subplot(1,2,2)
% plot(nan,nan)
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt-MI1_noadapt_asyn) nanmean(MI2_asynadapt-MI2_noadapt_asyn)], [nanstd(MI1_asynadapt-MI1_noadapt_asyn)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt-MI2_noadapt_asyn)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt-MI1_noadapt_cont) nanmean(MI2_contadapt-MI2_noadapt_cont)], [nanstd(MI1_contadapt-MI1_noadapt_cont)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt-MI2_noadapt_cont)./sqrt(size(MI2_contadapt,1))],'-o');
% legend('','Asynchronous','Contingent');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index (Post-Pre)')
% suptitle('Only cells responsive at max preferred contrast')
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ModulationSummary_RespCells.pdf']),'-dpdf','-bestfit');
% 
% T_noadapt = [nanmean(noadapt_resp_cell{1,4}(preftest_noadaptresp_ind,:),2)];
% TM1_noadapt = [nanmean(noadapt_resp_cell{2,4}(preftest_noadaptresp_ind,:),2)];
% TM2_noadapt = [nanmean(noadapt_resp_cell{3,4}(preftest_noadaptresp_ind,:),2)];
% T_noadapt(find(T_noadapt<0)) = 0;
% TM1_noadapt(find(TM1_noadapt<0)) = 0;
% TM2_noadapt(find(TM2_noadapt<0)) = 0;
% MI1_noadapt = (T_noadapt-TM1_noadapt)./(T_noadapt+TM1_noadapt);
% MI2_noadapt = (T_noadapt-TM2_noadapt)./(T_noadapt+TM2_noadapt);
% 
% noTM_asynadapt = [nanmean(asynadapt_resp_cell{1,1}(preftest_asynresp_ind,:),2)];
% T_asynadapt = [nanmean(asynadapt_resp_cell{1,4}(preftest_asynresp_ind,:),2)]-noTM_asynadapt;
% TM1_asynadapt = [nanmean(asynadapt_resp_cell{2,4}(preftest_asynresp_ind,:),2)]-noTM_asynadapt;
% TM2_asynadapt = [nanmean(asynadapt_resp_cell{3,4}(preftest_asynresp_ind,:),2)]-noTM_asynadapt;
% T_asynadapt(find(T_asynadapt<0)) = 0;
% TM1_asynadapt(find(TM1_asynadapt<0)) = 0;
% TM2_asynadapt(find(TM2_asynadapt<0)) = 0;
% MI1_asynadapt = (T_asynadapt-TM1_asynadapt)./(T_asynadapt+TM1_asynadapt);
% MI2_asynadapt = (T_asynadapt-TM2_asynadapt)./(T_asynadapt+TM2_asynadapt);
% 
% noTM_contadapt = [nanmean(contadapt_resp_cell{1,1}(preftest_contresp_ind,:),2)];
% T_contadapt = [nanmean(contadapt_resp_cell{1,4}(preftest_contresp_ind,:),2)]-noTM_contadapt;
% TM1_contadapt = [nanmean(contadapt_resp_cell{2,4}(preftest_contresp_ind,:),2)]-noTM_contadapt;
% TM2_contadapt = [nanmean(contadapt_resp_cell{3,4}(preftest_contresp_ind,:),2)]-noTM_contadapt;
% T_contadapt(find(T_contadapt<0)) = 0;
% TM1_contadapt(find(TM1_contadapt<0)) = 0;
% TM2_contadapt(find(TM2_contadapt<0)) = 0;
% MI1_contadapt = (T_contadapt-TM1_contadapt)./(T_contadapt+TM1_contadapt);
% MI2_contadapt = (T_contadapt-TM2_contadapt)./(T_contadapt+TM2_contadapt);
% 
% T_noadapt_asyn = [nanmean(noadapt_resp_cell{1,4}(preftest_asynresp_ind,:),2)];
% TM1_noadapt_asyn = [nanmean(noadapt_resp_cell{2,4}(preftest_asynresp_ind,:),2)];
% TM2_noadapt_asyn = [nanmean(noadapt_resp_cell{3,4}(preftest_asynresp_ind,:),2)];
% T_noadapt_asyn(find(T_noadapt_asyn<0)) = 0;
% TM1_noadapt_asyn(find(TM1_noadapt_asyn<0)) = 0;
% TM2_noadapt_asyn(find(TM2_noadapt_asyn<0)) = 0;
% MI1_noadapt_asyn = (T_noadapt_asyn-TM1_noadapt_asyn)./(T_noadapt_asyn+TM1_noadapt_asyn);
% MI2_noadapt_asyn = (T_noadapt_asyn-TM2_noadapt_asyn)./(T_noadapt_asyn+TM2_noadapt_asyn);
% 
% T_noadapt_cont = [nanmean(noadapt_resp_cell{1,4}(preftest_contresp_ind,:),2)];
% TM1_noadapt_cont = [nanmean(noadapt_resp_cell{2,4}(preftest_contresp_ind,:),2)];
% TM2_noadapt_cont = [nanmean(noadapt_resp_cell{3,4}(preftest_contresp_ind,:),2)];
% T_noadapt_cont(find(T_noadapt_cont<0)) = 0;
% TM1_noadapt_cont(find(TM1_noadapt_cont<0)) = 0;
% TM2_noadapt_cont(find(TM2_noadapt_cont<0)) = 0;
% MI1_noadapt_cont = (T_noadapt_cont-TM1_noadapt_cont)./(T_noadapt_cont+TM1_noadapt_cont);
% MI2_noadapt_cont = (T_noadapt_cont-TM2_noadapt_cont)./(T_noadapt_cont+TM2_noadapt_cont);
% 
% 
% figure; 
% subplot(1,2,1)
% errorbar([0.25 0.5], [nanmean(MI1_noadapt) nanmean(MI2_noadapt)], [nanstd(MI1_noadapt)./sqrt(size(MI1_noadapt,1)) nanstd(MI2_noadapt)./sqrt(size(MI2_noadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt) nanmean(MI2_asynadapt)], [nanstd(MI1_asynadapt)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt) nanmean(MI2_contadapt)], [nanstd(MI1_contadapt)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt)./sqrt(size(MI2_contadapt,1))],'-o');
% legend(['No adapt- n = ' num2str(size(MI1_noadapt,1))],['Asynchronous- n = ' num2str(size(MI1_asynadapt,1))],['Contingent - n = ' num2str(size(MI1_contadapt,1))],'location','southwest');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index')
% 
% subplot(1,2,2)
% plot(nan,nan)
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_asynadapt-MI1_noadapt_asyn) nanmean(MI2_asynadapt-MI2_noadapt_asyn)], [nanstd(MI1_asynadapt-MI1_noadapt_asyn)./sqrt(size(MI1_asynadapt,1)) nanstd(MI2_asynadapt-MI2_noadapt_asyn)./sqrt(size(MI2_asynadapt,1))],'-o');
% hold on
% errorbar([0.25 0.5], [nanmean(MI1_contadapt-MI1_noadapt_cont) nanmean(MI2_contadapt-MI2_noadapt_cont)], [nanstd(MI1_contadapt-MI1_noadapt_cont)./sqrt(size(MI1_contadapt,1)) nanstd(MI2_contadapt-MI2_noadapt_cont)./sqrt(size(MI2_contadapt,1))],'-o');
% legend('','Asynchronous','Contingent');
% xlim([0 1]);
% ylim([-1 1])
% xlabel('Preferred stim contrast')
% ylabel('Modulation Index (Post-Pre)')
% suptitle('Only cells responsive at max preferred contrast')
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ModulationSummary_RespCells_PrefTestOnly.pdf']),'-dpdf','-bestfit');


% noadapt_allcell_resp = zeros(nCells,nMask,nTest);
% asynadapt_allcell_resp = zeros(nCells,nMask,nTest);
% contadapt_allcell_resp = zeros(nCells,nMask,nTest);
% for im = 1:nMask
%     for it = 1:nTest
%         noadapt_allcell_resp(:,im,it) = nanmean(noadapt_resp_cell{im,it},2);
%         asynadapt_allcell_resp(:,im,it) = nanmean(asynadapt_resp_cell{im,it},2);
%         contadapt_allcell_resp(:,im,it) = nanmean(contadapt_resp_cell{im,it},2);
%     end
% end
            
% for iCell = 1:10
%     figure;
%     subplot(1,3,1)
%     plot(repmat(testCons', [1 3]), squeeze(noadapt_allcell_resp(iCell,:,:))');
%     ylim([-0.1 0.4])
%     subplot(1,3,2)
%     plot(repmat(testCons', [1 3]), squeeze(asynadapt_allcell_resp(iCell,:,:))');
%     title('Asynchronous')
%     ylim([-0.1 0.4])
%     subplot(1,3,3)
%     plot(repmat(testCons', [1 3]), squeeze(contadapt_allcell_resp(iCell,:,:))');
%     ylim([-0.1 0.4])
%     title('Contingent')
%     legend(num2str(maskCons'))
% end
