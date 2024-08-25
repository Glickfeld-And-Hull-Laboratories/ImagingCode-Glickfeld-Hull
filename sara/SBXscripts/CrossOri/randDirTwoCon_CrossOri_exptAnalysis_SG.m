
clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirTwoCon_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
iexp = 3;
frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';

fprintf([mouse ' ' date '\n'])

%% Pref direction analysis
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%%
if doRedChannel == 0
    red_cells = [];
end

if ~exist('sz','var')
    nCells = size(data_tc,2);
    nFrames = size(data_tc,1);
else
end


SF_all = celleqel2mat_padded(input.tStimOneGratingSpatialFreqCPD);
SFs = unique(SF_all);
nSF = length(SFs);
if nSF>1
    doSF = 1;
    nF = nSF;
else 
    doSF = 0;
    nF = 1;
end
TF_all = celleqel2mat_padded(input.tMaskOneGratingTemporalFreqCPS);
TFs = unique(TF_all);
nTF = length(TFs);
if nTF>1
    doTF = 1;
    nF = nTF;
else 
    doTF = 0;
    nF = 1;
end

prewin_frames = frame_rate;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

if exist('npSub_tc','var') == 0
    npSub_tc = data_tc;
    clear data_tc
end

if exist('sz', 'var')
    for itrial = 1:nTrials
        if cStimOn(itrial) + postwin_frames < sz(3)  %sz(3) should be = to nFrames
            data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        end
    end
else
    for itrial = 1:nTrials
        if cStimOn(itrial) + postwin_frames < nFrames
            data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
        end
    end
end


data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);

nMaskCon = length(maskCons);
nStimCon = length(stimCons);

ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
ind_maskAlone = intersect(find(stimCon_all==0),find(maskCon_all));
ind_plaid = intersect(find(stimCon_all>0),find(maskCon_all>0));
ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));

base_cell = cell(nStimDir,nMaskCon,nStimCon);
resp_cell = cell(nStimDir,nMaskCon,nStimCon);
trialInd = cell(nStimDir,nMaskCon,nStimCon);
trialsperstim = zeros(nStimDir,nMaskCon,nStimCon);
h_resp =zeros(nCells,nStimDir,nMaskCon,nStimCon);
p_resp =zeros(nCells,nStimDir,nMaskCon,nStimCon);
base_win = prewin_frames-ceil(frame_rate/2):prewin_frames;
resp_win = prewin_frames+5:prewin_frames+nFramesOn;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
data_dfof_con_ph_tc_avg = nan(prewin_frames+postwin_frames, nCells, nStimDir, nMaskCon, nStimCon);

avg_resp_dir = zeros(nCells,nStimDir,nMaskCon,nStimCon);
all_resp_dir = [];
all_resp_plaid = cell(1,nMaskCon,nStimCon);
all_dir = [];
all_plaid = cell(1,nMaskCon,nStimCon);
nStim = nStimDir;


for iDir = 1:nStimDir
    ind_stimdir = find(stimDir_all == stimDirs(iDir));
    ind_maskdir = find(maskDir_all == stimDirs(iDir));
    ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
    ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
    all_resp_dir = [all_resp_dir squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone),1))];
    all_dir = [all_dir iDir.*ones(size(ind_diralone))];
    for im = 1:nMaskCon
        ind_mask = find(maskCon_all == maskCons(im));
        for it = 1:nStimCon
            if maskCons(im) == stimCons(it)
                ind_stim = find(stimCon_all == stimCons(it));
                ind = intersect(ind_dirplaid, intersect(ind_stim,ind_mask));
                trialsperstim(iDir,im,it) = length(ind);
                resp_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                base_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
                [h_resp(:,iDir,im,it), p_resp(:,iDir,im,it)] = ttest2(resp_cell{iDir,im,it},base_cell{iDir,im,it},'dim',2,'tail','right','alpha', 0.05./nStim);
                avg_resp_dir(:,iDir,im,it,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind),1),3));
                avg_resp_dir(:,iDir,im,it,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind),1),[],3)./sqrt(length(ind)));
                trialInd{iDir,im,it} = ind;
            elseif maskCons(im)==0;
                ind_stim = find(stimCon_all == stimCons(it));
                ind = intersect(ind_diralone, intersect(ind_stim,ind_mask));
                trialsperstim(iDir,im,it) = length(ind);
                resp_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
                base_cell{iDir,im,it} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
                data_dfof_con_ph_tc_avg(:,:,iDir,im,it,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
                
                avg_resp_dir(:,iDir,im,it,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind),1),3));
                avg_resp_dir(:,iDir,im,it,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind),1),[],3)./sqrt(length(ind)));
                trialInd{iDir,im,it} = ind;
            end
        end
    end
end
 

respgratlow_ind = find(sum(h_resp(:,:,1,1),2));
respgrathigh_ind = find(sum(h_resp(:,:,1,2),2));
respplaidlow_ind = find(sum(h_resp(:,:,2,1),2));
respplaidhigh_ind = find(sum(h_resp(:,:,3,2),2));

resp_ind = unique([respgratlow_ind; respgrathigh_ind; respplaidlow_ind; respplaidhigh_ind]);

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof_con_ph_tc_avg', 'resp_cell','base_cell', 'data_dfof_tc', 'resp_ind', 'respgratlow_ind', 'respgrathigh_ind', 'respplaidlow_ind', 'respplaidhigh_ind', 'frame_rate');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'base_win', 'stimCon_all', 'maskCon_all', 'maskPhas_all', 'stimCons', 'maskCons', 'maskPhas', 'nStimCon', 'nMaskCon', 'nMaskPhas','TF_all','nTF','TFs');

