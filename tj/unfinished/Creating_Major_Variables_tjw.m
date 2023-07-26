clear all
clear global
%% 
mouse = 'i1803';
ref_date = '210802';
day2 = '210804';
day3 = '210806';
ImgFolder = strvcat('003');
ImgFolder2 = strvcat('003');
ImgFolder3 = strvcat('003');
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str3 = catRunName(ImgFolder3, nrun);
run_str2 = catRunName(ImgFolder2, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
% code works differently if there is a day3. Newer datasets have a day 3
% but some older ones don't. Set day3x = 1 if there is a day3
day3x = 1;

%% load data
oriTuning_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_oriTuningAndFits.mat']));
oriTuning_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_oriTuningAndFits.mat']));
oriTuning_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_oriTuningAndFits.mat']));

TCs_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
TCs_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_TCs.mat']));
TCs_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_TCs.mat']));

trialData_D1 = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_trialData.mat']));
trialData_D2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_trialData.mat']));
trialData_D3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_trialData.mat']));

%% loading stuff 

% day 1
load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']));
tGratingDir1 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs1 = unique(tGratingDir1);
nDir1 = length(dirs1);
nOn1 = input.nScansOn;
nOff1 = input.nScansOff;
mask = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
mask_cell = mask.mask_cell;
data_dfof1 = trialData_D1.data_dfof;
nTrials = size(data_dfof1,3);
npSub_tc1 = TCs_D1.npSub_tc;
nFrames1 = size(npSub_tc1,1);
nCells1 = size(npSub_tc1,2);
dfof = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']));
dfof_max = dfof.data_dfof_max;
pixel = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_pixel.mat']));
pix = pixel.pix;

% day 2
load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_input.mat']));
tGratingDir2 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs2 = unique(tGratingDir1);
nDir2 = length(dirs2);
nOn2 = input.nScansOn;
nOff2 = input.nScansOff;
mask2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_mask_cell.mat']));
mask_cell2 = mask2.mask_cell;
data_dfof2 = trialData_D2.data_dfof;
nTrials2 = size(data_dfof2,3);
np2 = TCs_D2.cellTCs_match;
npSub_tc2 = np2{2};
nFrames2 = size(npSub_tc2,1);
nCells2 = size(npSub_tc2,2);
cell_list2 = TCs_D2.match_ind;

dir_mat = tGratingDir1;
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
nOri = length(oris);
dirs = unique(dir_mat);
nDir = length(dirs);

dir_mat2 = tGratingDir2;
ori_mat2 = dir_mat2;
ori_mat2(find(dir_mat2>=180)) = ori_mat2(dir_mat2>=180)-180;
oris2 = unique(ori_mat2);
nOri2 = length(oris2);
dirs2 = unique(dir_mat2);
nDir2 = length(dirs2);

if day3x
% other stuff from day1 and day2
align = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_multiday_alignment.mat']));
redChImg = align.redChImg;
align2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_multiday_alignment.mat']));
redChImg2 = align2.redChImg;
trans2 = align2.fitGeoTAf;
dfof_max2 = align2.dfmax;
dfof_max2 = dfof_max2{3};
pix2 = align2.corrmap;
pix2 = pix2{3};

% day 3
load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_input.mat']));
tGratingDir3 = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs3 = unique(tGratingDir1);
nDir3 = length(dirs3);
nOn3 = input.nScansOn;
nOff3 = input.nScansOff;
mask3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_mask_cell.mat']));
mask_cell3 = mask3.mask_cell;
data_dfof3 = trialData_D3.data_dfof;
nTrials3 = size(data_dfof3,3);
np3 = TCs_D3.cellTCs_match;
npSub_tc3 = np3{2};
nFrames3 = size(npSub_tc3,1);
nCells3 = size(npSub_tc3,2);
cell_list3 = TCs_D3.match_ind;
align3 = load(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_multiday_alignment.mat']));
trans3 = align3.fitGeoTAf;
dfof_max3 = align3.dfmax;
dfof_max3 = dfof_max3{3};
pix3 = align3.corrmap;
pix3 = pix3{3};
redChImg3 = align3.redChImg;

dir_mat3 = tGratingDir3;
ori_mat3 = dir_mat3;
ori_mat3(find(dir_mat3>=180)) = ori_mat3(dir_mat3>=180)-180;
oris3 = unique(ori_mat3);
nOri3 = length(oris3);
dirs3 = unique(dir_mat3);
nDir3 = length(dirs3);
end
%% define well fit cells
if day3x
cells = intersect(cell_list2,cell_list3);
goodfit_D1 = find(oriTuning_D1.fitReliability<22.5);
goodfit_D2 = find(oriTuning_D2.fitReliability<22.5);
goodfit_D3 = find(oriTuning_D3.fitReliability<22.5);
goodfit1 = intersect(goodfit_D1,goodfit_D2);
goodfit = intersect(goodfit1,goodfit_D3);
cell_list = intersect(cells,goodfit);
nCells = length(cell_list);
else
goodfit_D1 = find(oriTuning_D1.fitReliability<22.5);
goodfit_D2 = find(oriTuning_D2.fitReliability<22.5);
goodfit = intersect(goodfit_D1,goodfit_D2);
cell_list = intersect(cell_list2,goodfit);
nCells = length(cell_list);
end
%% define variables
avgResp_D1 = oriTuning_D1.avgResponseEaOri(cell_list,:);
avgResp_D2 = oriTuning_D2.avgResponseEaOri(cell_list,:);

semResp_D1 = oriTuning_D1.semResponseEaOri(cell_list,:);
semResp_D2 = oriTuning_D2.semResponseEaOri(cell_list,:);

vonMisesFit_D1 = oriTuning_D1.vonMisesFitAllCellsAllBoots(:,:,cell_list);
vonMisesFit_D2 = oriTuning_D2.vonMisesFitAllCellsAllBoots(:,:,cell_list);

[maxResp_D1 prefOri_D1] = max(squeeze(oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);
[maxResp_D2 prefOri_D2] = max(squeeze(oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);

dfof_dir1 = trialData_D1.dir_resp(cell_list,:);
dfof_dir2 = trialData_D2.dir_resp(cell_list,:);

if day3x
avgResp_D3 = oriTuning_D3.avgResponseEaOri(cell_list,:);
semResp_D3 = oriTuning_D3.semResponseEaOri(cell_list,:);
vonMisesFit_D3 = oriTuning_D3.vonMisesFitAllCellsAllBoots(:,:,cell_list);
[maxResp_D3 prefOri_D3] = max(squeeze(oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);
dfof_dir3 = trialData_D3.dir_resp(cell_list,:);
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_prefOri.mat']),'prefOri_D1','prefOri_D2','prefOri_D3','avgResp_D1','avgResp_D2','avgResp_D3','vonMisesFit_D1','vonMisesFit_D2','vonMisesFit_D3');
else
prefOri_D3 = [];
avgResp_D3 = [];
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_prefOri.mat']),'prefOri_D1','prefOri_D2','prefOri_D3','avgResp_D1','avgResp_D2','avgResp_D3','vonMisesFit_D1','vonMisesFit_D2','vonMisesFit_D3');
end

%% creating variables for average tuning curve - what is this??
[maxResp prefOri_ind] = max(avgResp_D1,[],2);
newAvgD1 = zeros(nCells,nOri);
for iCell = 1:nCells
    if prefOri_ind(iCell)<4
        index = 4 - prefOri_ind(iCell);
        newAvgD1(iCell,:) = avgResp_D1(iCell,[8-index+1:8 1:8-index]);
    elseif prefOri_ind(iCell)>4
        index = prefOri_ind(iCell) - 4;
        newAvgD1(iCell,:) = avgResp_D1(iCell,[index+1:8 1:index]);
    elseif prefOri_ind(iCell)==4
        newAvgD1(iCell,:) = avgResp_D1(iCell,:);
    end
end
[maxResp prefOri_ind] = max(avgResp_D2,[],2);
newAvgD2 = zeros(nCells,nOri);
for iCell = 1:nCells
    if prefOri_ind(iCell)<4
        index = 4 - prefOri_ind(iCell);
        newAvgD2(iCell,:) = avgResp_D2(iCell,[8-index+1:8 1:8-index]);
    elseif prefOri_ind(iCell)>4
        index = prefOri_ind(iCell) - 4;
        newAvgD2(iCell,:) = avgResp_D2(iCell,[index+1:8 1:index]);
    elseif prefOri_ind(iCell)==4
        newAvgD2(iCell,:) = avgResp_D2(iCell,:);
    end
end
if day3x
[maxResp prefOri_ind] = max(avgResp_D3,[],2);
newAvgD3 = zeros(nCells,nOri);
for iCell = 1:nCells
    if prefOri_ind(iCell)<4
        index = 4 - prefOri_ind(iCell);
        newAvgD3(iCell,:) = avgResp_D3(iCell,[8-index+1:8 1:8-index]);
    elseif prefOri_ind(iCell)>4
        index = prefOri_ind(iCell) - 4;
        newAvgD3(iCell,:) = avgResp_D3(iCell,[index+1:8 1:index]);
    elseif prefOri_ind(iCell)==4
        newAvgD3(iCell,:) = avgResp_D3(iCell,:);
    end
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_newAvg.mat']),'newAvgD1','newAvgD2','newAvgD3');
else
newAvgD3 = [];
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_newAvg.mat']),'newAvgD1','newAvgD2','newAvgD3');
end

%% DSI/OSI
orthoOri_D1 = zeros(nCells,1);
orthoResp_D1 = zeros(nCells,1);
orthoOri_D2 = zeros(nCells,1);
orthoResp_D2 = zeros(nCells,1);
orthoOri_D3 = zeros(nCells,1);
orthoResp_D3 = zeros(nCells,1);
[maxResp prefOri_ind] = max(avgResp_D1,[],2);
for iCell = 1:nCells
if prefOri_ind(iCell) <= 4
  orthoOri_D1(iCell) = prefOri_ind(iCell) + 4;
elseif prefOri_ind(iCell) > 4
  orthoOri_D1(iCell) = prefOri_ind(iCell) - 4;
end
orthoResp_D1(iCell) = avgResp_D1(iCell,orthoOri_D1(iCell));
end
orthoResp_D1(find(orthoResp_D1<0)) = 0;
OSI1 = (maxResp_D1'-orthoResp_D1)./(maxResp_D1'+orthoResp_D1);

[maxResp prefOri_ind] = max(avgResp_D2,[],2);
for iCell = 1:nCells
if prefOri_ind(iCell) <= 4
  orthoOri_D2(iCell) = prefOri_ind(iCell) + 4;
elseif prefOri_ind(iCell) > 4
  orthoOri_D2(iCell) = prefOri_ind(iCell) - 4;
end
orthoResp_D2(iCell) = avgResp_D2(iCell,orthoOri_D2(iCell));
end
orthoResp_D2(find(orthoResp_D2<0)) = 0;
OSI2 = (maxResp_D2'-orthoResp_D2)./(maxResp_D2'+orthoResp_D2);

if day3x
[maxResp prefOri_ind] = max(avgResp_D3,[],2);
for iCell = 1:nCells
if prefOri_ind(iCell) <= 4
  orthoOri_D3(iCell) = prefOri_ind(iCell) + 4;
elseif prefOri_ind(iCell) > 4
  orthoOri_D3(iCell) = prefOri_ind(iCell) - 4;
end
orthoResp_D3(iCell) = avgResp_D3(iCell,orthoOri_D3(iCell));
end
orthoResp_D3(find(orthoResp_D3<0)) = 0;
OSI3 = (maxResp_D3'-orthoResp_D3)./(maxResp_D3'+orthoResp_D3);
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_OSI.mat']),'OSI1','OSI2','OSI3');
else
OSI3 = [];
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_OSI.mat']),'OSI1','OSI2','OSI3');
end

oppDir_D1 = zeros(nCells,1);
oppResp_D1 = zeros(nCells,1);
oppDir_D2 = zeros(nCells,1);
oppResp_D2 = zeros(nCells,1);
oppDir_D3 = zeros(nCells,1);
oppResp_D3 = zeros(nCells,1);
[maxResp1 prefDir_ind] = max(dfof_dir1,[],2);
for iCell = 1:nCells
if prefDir_ind(iCell) <= 8
  oppDir_D1(iCell) = prefDir_ind(iCell) + 8;
elseif prefDir_ind(iCell) > 8
  oppDir_D1(iCell) = prefDir_ind(iCell) - 8;
end
oppResp_D1(iCell) = dfof_dir1(iCell,oppDir_D1(iCell));
end
oppResp_D1(find(oppResp_D1<0)) = 0;
DSI1 = (maxResp1-oppResp_D1)./(maxResp1+oppResp_D1);

[maxResp2 prefDir_ind] = max(dfof_dir2,[],2);
for iCell = 1:nCells
if prefDir_ind(iCell) <= 8
  oppDir_D2(iCell) = prefDir_ind(iCell) + 8;
elseif prefDir_ind(iCell) > 8
  oppDir_D2(iCell) = prefDir_ind(iCell) - 8;
end
oppResp_D2(iCell) = dfof_dir2(iCell,oppDir_D2(iCell));
end
oppResp_D2(find(oppResp_D2<0)) = 0;
DSI2 = (maxResp2-oppResp_D2)./(maxResp2+oppResp_D2);

if day3x
[maxResp3 prefDir_ind] = max(dfof_dir3,[],2);
for iCell = 1:nCells
if prefDir_ind(iCell) <= 8
  oppDir_D3(iCell) = prefDir_ind(iCell) + 8;
elseif prefDir_ind(iCell) > 8
  oppDir_D3(iCell) = prefDir_ind(iCell) - 8;
en
oppResp_D3(iCell) = dfof_dir3(iCell,oppDir_D3(iCell));
end
oppResp_D3(find(oppResp_D3<0)) = 0;
DSI3 = (maxResp3-oppResp_D3)./(maxResp3+oppResp_D3);
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_DSI.mat']),'DSI1','DSI2','DSI3');
else
DSI3 = [];
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_DSI.mat']),'DSI1','DSI2','DSI3');
end

%% Getting response average for each direction for three groups of trials
resp_wind = nOff1+10:nOff1+nOn1;
resp = squeeze(mean(data_dfof1(resp_wind,:,:),1));
resp2 = squeeze(mean(data_dfof2(resp_wind,:,:),1));
resp3 = squeeze(mean(data_dfof3(resp_wind,:,:),1));
ind{1} = randperm(nTrials,nTrials/2);
ind{2} = setdiff(1:nTrials, ind{1});
ind{3} = 1:nTrials;
dir_resp_avg = zeros(nCells,nDir,3);
dir_resp_avg2 = zeros(nCells,nDir,3);
dir_resp_avg3 = zeros(nCells,nDir,3);
dir_resp_max = zeros(nCells,nDir,3);
dir_resp_max2 = zeros(nCells,nDir,3);
dir_resp_max3 = zeros(nCells,nDir,3);
for i = 1:3
  for iDir = 1:nDir
    ind_dir = intersect(ind{i},find(dir_mat == dirs(iDir)));
    dir_resp_avg(:,iDir,i) = squeeze(mean(resp(cell_list,ind_dir),2));
    dir_resp_max(:,iDir,i) = squeeze(max(resp(cell_list,ind_dir),[],2));
  end
end
for i = 1:3
  for iDir = 1:nDir
    ind_dir = intersect(ind{i},find(dir_mat2 == dirs(iDir)));
    dir_resp_avg2(:,iDir,i) = squeeze(mean(resp2(cell_list,ind_dir),2)); 
    dir_resp_max2(:,iDir,i) = squeeze(max(resp2(cell_list,ind_dir),[],2));
  end
end
if day3x
    for i = 1:3
      for iDir = 1:nDir
        ind_dir = intersect(ind{i},find(dir_mat3 == dirs(iDir)));
        dir_resp_avg3(:,iDir,i) = squeeze(mean(resp3(cell_list,ind_dir),2)); 
        dir_resp_max3(:,iDir,i) = squeeze(max(resp3(cell_list,ind_dir),[],2));
      end
    end
else
dir_resp_avg3 = [];
dir_resp_max3 = [];
end

%% selectivity
% transform all of the negative responses to zero
dir_resp_avg_rect = dir_resp_avg(:,:,:);
dir_resp_avg_rect(find(dir_resp_avg(:,:,:)<0))=0;
dir_resp_avg2_rect = dir_resp_avg2(:,:,:);
dir_resp_avg2_rect(find(dir_resp_avg2(:,:,:)<0))=0;
dir_resp_avg3_rect = dir_resp_avg3(:,:,:);
dir_resp_avg3_rect(find(dir_resp_avg3(:,:,:)<0))=0;
f_max = squeeze(max(dir_resp_avg_rect,[],2));
f_max2 = squeeze(max(dir_resp_avg2_rect,[],2));
f_max3 = squeeze(max(dir_resp_avg3_rect,[],2));
c_pass = NaN(nCells,100);
c_pass2 = NaN(nCells,100);
c_pass3 = NaN(nCells,100);
S = NaN(nCells,1);
S2 = NaN(nCells,1);
S3 = NaN(nCells,1);
t = NaN(nCells,100);
t2 = NaN(nCells,100);
t3 = NaN(nCells,100);
[n n2] = subplotn(nCells);
% create 100 threshold values and quantify the responses that passes each
% thresh
figure;
for iCell = 1:nCells
    t(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t(iCell,it);
        c_pass(iCell,it) = sum(dir_resp_avg_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S(iCell) = 1 - (2*(sum(c_pass(iCell,:),2)/100));
    subplot(n,n2,iCell)
    bar(t(iCell,:),c_pass(iCell,:),'r')
    title(['cell ' num2str(iCell) ' S ' num2str(S(iCell))])
    xlabel(['Threshold'])
    ylabel(['Responses'])
    print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['threshold.pdf']),'-dpdf', '-bestfit')
end
figure;
for iCell = 1:nCells
    t2(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t2(iCell,it);
    c_pass2(iCell,it) = sum(dir_resp_avg2_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S2(iCell) = 1 - (2*(sum(c_pass2(iCell,:),2)/100));
    subplot(n,n2,iCell)
    bar(t2(iCell,:),c_pass2(iCell,:),'r')
    title(['cell ' num2str(iCell) ' S ' num2str(S2(iCell))])
    xlabel(['Threshold'])
    ylabel(['Responses'])
    print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], ['threshold2.pdf']),'-dpdf', '-bestfit')
end
if day3x
figure;
for iCell = 1:nCells
    t3(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t3(iCell,it);
    c_pass3(iCell,it) = sum(dir_resp_avg3_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S3(iCell) = 1 - (2*(sum(c_pass3(iCell,:),2)/100));
    subplot(n,n2,iCell)
    bar(t3(iCell,:),c_pass3(iCell,:),'r')
    title(['cell ' num2str(iCell) ' S ' num2str(S3(iCell))])
    xlabel(['Threshold'])
    ylabel(['Responses'])
    print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['threshold3.pdf']),'-dpdf', '-bestfit')
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_threshold.mat']),'t','c_pass','t2','c_pass2','t3','c_pass3','S','S2','S3');
else
    t3 = [];
    c_pass3 = [];
    S3 = [];
    save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_threshold.mat']),'t','c_pass','t2','c_pass2','t3','c_pass3','S','S2','S3');
end

%% ttest
base_wind = 1+nOff1-nOn1:nOff1;
dfof_base = squeeze(mean(data_dfof1(base_wind,cell_list,:),1))';
resp_wind = nOff1+10:nOff1+nOn1;
dfof_resp = squeeze(mean(data_dfof1(resp_wind,cell_list,:),1))';
h1 = NaN(nDir,nCells);
p1 = NaN(nDir,nCells);
for iDir = 1:nDir
    ind = find(tGratingDir1==dirs(iDir));
    x1 = dfof_base(ind,:);
    y1 = dfof_resp(ind,:);
    [h1(iDir,:),p1(iDir,:)] = ttest(x1,y1,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
end
dfof_base = squeeze(mean(data_dfof2(base_wind,cell_list,:),1))';
dfof_resp = squeeze(mean(data_dfof2(resp_wind,cell_list,:),1))';
h2 = NaN(nDir,nCells);
p2 = NaN(nDir,nCells);
for iDir = 1:nDir
    ind = find(tGratingDir2==dirs(iDir));
    x2 = dfof_base(ind,:);
    y2 = dfof_resp(ind,:);
    [h2(iDir,:),p2(iDir,:)] = ttest(x2,y2,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
end
if day3x
dfof_base = squeeze(mean(data_dfof3(base_wind,cell_list,:),1))';
dfof_resp = squeeze(mean(data_dfof3(resp_wind,cell_list,:),1))';
h3 = NaN(nDir,nCells);
p3 = NaN(nDir,nCells);
for iDir = 1:nDir
    ind = find(tGratingDir3==dirs(iDir));
    x3 = dfof_base(ind,:);
    y3 = dfof_resp(ind,:);
    [h3(iDir,:),p3(iDir,:)] = ttest(x3,y3,'dim',1,'Alpha',0.05./(nDir-1),'tail','left');
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_ttest.mat']),'h1','h2','h3','p1','p2','p3');
else
    h3 = [];
    p3 = []
    save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_ttest.mat']),'h1','h2','h3','p1','p2','p3');
end
%% signal correlation
signal_corr_half1 = NaN(nCells,1);
signal_corr_half2 = NaN(nCells,1);
signal_corr_half3 = NaN(nCells,1);
signal_corr_day2 = NaN(nCells,1);
signal_corr_day3 = NaN(nCells,1);
signal_corr_day_half2 = NaN(nCells,1);
signal_corr_day_half3 = NaN(nCells,1);
for iCell = 1:nCells
    for i = 1:3
        signal_corr_half1(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1),dir_resp_avg(iCell,:,2)));
        signal_corr_half2(iCell,1) = triu2vec(corrcoef(dir_resp_avg2(iCell,:,1),dir_resp_avg2(iCell,:,2)));
        signal_corr_day2(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,3),dir_resp_avg2(iCell,:,3)));
        signal_corr_day_half2(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1),dir_resp_avg2(iCell,:,1)));
        if day3x
        signal_corr_half3(iCell,1) = triu2vec(corrcoef(dir_resp_avg3(iCell,:,1),dir_resp_avg3(iCell,:,2)));
        signal_corr_day3(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,3),dir_resp_avg3(iCell,:,3)));
        signal_corr_day_half3(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1),dir_resp_avg3(iCell,:,1)));
        else
            signal_corr_half3 = [];
            signal_corr_day3 = [];
            signal_corr_day_half3 = []
        end
    end
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_signalCorr.mat']),'signal_corr_half1','signal_corr_half2','signal_corr_half3','signal_corr_day2','signal_corr_day3','signal_corr_day_half2','signal_corr_day_half3');

%% vonmises fit
theta = deg2rad(dirs); 
thetafine = deg2rad(1:360); 
for iCell = 1:nCells
    try
    [b_hat(iCell),k1_hat(iCell),R1_hat(iCell),R2_hat(iCell),u1_hat(iCell),u2_hat(iCell),sse(iCell),R_square(iCell)] = miaovonmisesfit_dir(theta,dir_resp_avg(iCell,:,3));
    y_fit(iCell,:) = b_hat(iCell)+R1_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u1_hat(iCell))-1))+R2_hat(iCell).*exp(k1_hat(iCell).*(cos(thetafine-u2_hat(iCell))-1));
    catch
        b_hat(iCell) = NaN;
        k1_hat(iCell) = NaN;
        R1_hat(iCell) = NaN;
        R2_hat(iCell) = NaN;
        u1_hat(iCell) = NaN;
        u2_hat(iCell) = NaN;
        sse_hat(iCell) = NaN;
        R_square_hat(iCell)= NaN; 
        y_fit(iCell,:) = NaN(1,length(thetafine));
    end
end

for iCell = 1:nCells
    try
    [b_hat2(iCell),k1_hat2(iCell),R1_hat2(iCell),R2_hat2(iCell),u1_hat2(iCell),u2_hat2(iCell),sse2(iCell),R_square2(iCell)] = miaovonmisesfit_dir(theta,dir_resp_avg2(iCell,:,3));
    y_fit2(iCell,:) = b_hat2(iCell)+R1_hat2(iCell).*exp(k1_hat2(iCell).*(cos(thetafine-u1_hat2(iCell))-1))+R2_hat2(iCell).*exp(k1_hat2(iCell).*(cos(thetafine-u2_hat2(iCell))-1));
    catch
        b_hat2(iCell) = NaN;
        k1_hat2(iCell) = NaN;
        R1_hat2(iCell) = NaN;
        R2_hat2(iCell) = NaN;
        u1_hat2(iCell) = NaN;
        u2_hat2(iCell) = NaN;
        sse_hat2(iCell) = NaN;
        R_square_hat2(iCell)= NaN; 
        y_fit2(iCell,:) = NaN(1,length(thetafine));
    end
end
if day3x
for iCell = 1:nCells
    try
    [b_hat3(iCell),k1_hat3(iCell),R1_hat3(iCell),R2_hat3(iCell),u1_hat3(iCell),u2_hat3(iCell),sse3(iCell),R_square3(iCell)] = miaovonmisesfit_dir(theta,dir_resp_avg3(iCell,:,3));
    y_fit3(iCell,:) = b_hat3(iCell)+R1_hat3(iCell).*exp(k1_hat3(iCell).*(cos(thetafine-u1_hat3(iCell))-1))+R2_hat3(iCell).*exp(k1_hat3(iCell).*(cos(thetafine-u2_hat3(iCell))-1));
    catch
        b_hat3(iCell) = NaN;
        k1_hat3(iCell) = NaN;
        R1_hat3(iCell) = NaN;
        R2_hat3(iCell) = NaN;
        u1_hat3(iCell) = NaN;
        u2_hat3(iCell) = NaN;
        sse_hat3(iCell) = NaN;
        R_square_hat3(iCell)= NaN; 
        y_fit3(iCell,:) = NaN(1,length(thetafine));
    end
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_vonmises.mat']),'k1_hat','k1_hat2','k1_hat3');
else
    k1_hat3 = [];
    save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_vonmises.mat']),'k1_hat','k1_hat2','k1_hat3');
end

%% tuning curves for 3 days overlapped
figure;
[n n2] = subplotn(nCells);
start = 1;
x = 0;
for iCell = 1:nCells
            subplot(n,n2,iCell-(x.*49))
            errorbar(oris,avgResp_D1(iCell,:), semResp_D1(iCell,:),'-o','color',[1 0 0])
            hold on
            plot(0:180,vonMisesFit_D1(:,1,iCell));
            hold on
            errorbar(oris,avgResp_D2(iCell,:), semResp_D2(iCell,:),'-o','color',[1 .54 0])
            hold on
            plot(0:180,vonMisesFit_D2(:,1,iCell));
            hold on
            errorbar(oris,avgResp_D3(iCell,:), semResp_D3(iCell,:),'-o','color',[1 .83 0])
            hold on
            plot(0:180,vonMisesFit_D3(:,1,iCell));
            hold on
            title(num2str(iCell))
            start = start+1;
            axis square
end
legend({'day1','day2','day3'})
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_tuning_curves.pdf']),'-dpdf', '-bestfit')

%% create a grid of cell squares showing pixel and dfof images
cell_stats = regionprops(mask_cell);
cell_stats = cell_stats(cell_list);
cell_stats2 = regionprops(mask_cell2);
cell_stats2 = cell_stats2(cell_list);
cell_stats3 = regionprops(mask_cell3);
cell_stats3 = cell_stats3(cell_list);
figure;
start = 1;
for iCell = 1:nCells
    width = 24; height = 24;
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xLeft(iCell) = (xCenter - width/2);
    yBottom(iCell) = (yCenter - height/2);
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft2(iCell) = (xCenter2 - width/2);
    yBottom2(iCell) = (yCenter2 - height/2);
    xCenter3 = round(cell_stats3(iCell).Centroid(2));
    yCenter3 = round(cell_stats3(iCell).Centroid(1));
    xLeft3(iCell) = (xCenter3 - width/2);
    yBottom3(iCell) = (yCenter3 - height/2);
    cell_mask = mask_cell(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
    cell_mask2 = mask_cell2(xLeft2(iCell):(xLeft2(iCell)+width),yBottom2(iCell):(height+yBottom2(iCell)));
    cell_mask3 = mask_cell3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
    if xLeft(iCell) > 12 && xLeft(iCell) < 488 && yBottom(iCell) > 12 && yBottom(iCell) < 772
    subplot(nCells,6,start);
        x = pix(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
        imagesc(x)
        pos = get(gca, 'Position');
        pos(1) = 0.025;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('pix')
    subplot(nCells,6,start+1);
        y = dfof_max(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
        imagesc(y)
        pos = get(gca, 'Position');
        pos(1) = 0.075;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('dfof')
    subplot(nCells,6,start+2);
        z = pix2(xLeft2(iCell):(xLeft2(iCell)+width),yBottom2(iCell):(height+yBottom2(iCell)),:);
        imagesc(z)
        pos = get(gca, 'Position');
        pos(1) = 0.125;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask2(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('pix')
    subplot(nCells,6,start+3);
        g = dfof_max2(xLeft2(iCell):(xLeft2(iCell)+width),yBottom2(iCell):(height+yBottom2(iCell)));
        imagesc(g)
        pos = get(gca, 'Position');
        pos(1) = 0.175;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask2(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('dfof')
    subplot(nCells,6,start+4);
        h = pix3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
        imagesc(h)
        pos = get(gca, 'Position');
        pos(1) = 0.225;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask3(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('red')
    subplot(nCells,6,start+5);
        v = dfof_max3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
        imagesc(v)
        pos = get(gca, 'Position');
        pos(1) = 0.275;
        pos(3) = 0.05;
        set(gca, 'Position', pos)
        axis off
        axis square
        hold on
        bound = cell2mat(bwboundaries(cell_mask3(:,:,1)));
        plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
        title('dfof')
   
    start = start+6;
    end
end
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_3day_cell_squares.pdf']),'-dpdf', '-bestfit')

%% dfof image with cells boxed and numbered
figure; imagesc(dfof_max)
hold on
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    text(cell_stats(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
    'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end
title('day1')
figure; imagesc(dfof_max2)
hold on
bound = cell2mat(bwboundaries(mask_cell2(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    text(cell_stats2(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
    'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
end
title('day2')
figure; imagesc(dfof_max3)
hold on
bound = cell2mat(bwboundaries(mask_cell3(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    text(cell_stats3(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
    'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
title('day3')
end


%% 4 main scatter plots
sz = 100;
figure;scatter(signal_corr_half1,signal_corr_day_half2,sz,'filled')
hold on
scatter(signal_corr_half1,signal_corr_day_half3,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('signal correlation within D1')
ylabel('signal correlation D1 vs DX')
refline(1,0)
legend({'within D1 session' 'between sessions'})
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['signal_corr_scat1.pdf']),'-dpdf', '-bestfit')

sz = 100;
figure;scatter(signal_corr_day2,signal_corr_day3,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('signal correlation D1 vs D2')
ylabel('signal correlation D1 vs D3')
refline(1,0)
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['signal_corr_scat2.pdf']),'-dpdf', '-bestfit')

figure;scatter(maxResp_D1,maxResp_D2,sz,'filled')
hold on;scatter(maxResp_D1,maxResp_D3,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('max resp D1')
ylabel('max resp DX')
refline(1,0)
legend({'day 2' 'day 3'})
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['maxresp_scat.pdf']),'-dpdf', '-bestfit')

figure;scatter(prefOri_D1,prefOri_D2,sz,'filled')
hold on;
scatter(prefOri_D1,prefOri_D3,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('pref ori D1')
ylabel('pref ori DX')
legend({'Day 2' 'Day 3'})
refline(1,0)
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['prefOri_scat.pdf']),'-dpdf', '-bestfit')

%% noise correlation
% day1
dirs = unique(dir_mat);
ind = cell(1,3);
ind{1} = randperm(nTrials,nTrials/2);
ind{2} = setdiff(1:nTrials, ind{1});
ind{3} = 1:nTrials;
half1_noise = zeros(nCells1,nOri,2,3);
for i = 1:3
    for idir = 1:nDir
    ind_dir = intersect(ind{i}, find(dir_mat == dirs(idir)));
    for iI = 1:length(ind_dir)
    iind = ind_dir(iI);
    half1_noise(:,iind,1,i) = resp(:,iind) - mean(resp(:,ind_dir),2);
    half1_noise(:,iind,2,i) = squeeze(std(resp(:,ind_dir),[],2))./sqrt(length(ind_dir));
    end
    end
end
% noise correlations between pairs of cells on the same day
half_noise_corr1 = NaN(nCells1,nCells1);
half_noise_corr2 = NaN(nCells1,nCells1);
full_noise_corr1 = NaN(nCells1,nCells1);
for iCell = 1:nCells1
    for iC = 1:nCells1
    half_noise_corr1(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,1),half1_noise(iC,:,1,1)));
    half_noise_corr2(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,2),half1_noise(iC,:,1,2)));
    full_noise_corr1(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,3),half1_noise(iC,:,1,3)));
    end
end
% day2
dirs = unique(dir_mat2);
ind{1} = randperm(nTrials2,nTrials2/2);
ind{2} = setdiff(1:nTrials2, ind{1});
ind{3} = 1:nTrials2;
half2_noise = zeros(nCells1,nOri,2,3);
for i = 1:3
    for idir = 1:nOri
    ind_dir = intersect(ind{i}, find(dir_mat2 == dirs(idir)));
    for iI = 1:length(ind_dir)
    iind = ind_dir(iI);
    half2_noise(:,iind,1,i) = resp2(:,iind) - mean(resp2(:,ind_dir),2);
    half1_noise(:,iind,2,i) = squeeze(std(resp(:,ind_dir),[],2))./sqrt(length(ind_ori));
    end
    end
end
half_noise_corr3 = NaN(nCells1,nCells1);
half_noise_corr4 = NaN(nCells1,nCells1);
full_noise_corr2 = NaN(nCells1,nCells1);
for iCell = 1:nCells1
    for iC = 1:nCells1
    half_noise_corr3(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,1),half2_noise(iC,:,1,1)));
    half_noise_corr4(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,2),half2_noise(iC,:,1,2)));
    full_noise_corr2(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,3),half2_noise(iC,:,1,3)));
    end
end
% day3
dirs = unique(dir_mat3);
ind{1} = randperm(nTrials2,nTrials2/2);
ind{2} = setdiff(1:nTrials2, ind{1});
ind{3} = 1:nTrials2;
half2_noise = zeros(nCells1,nOri,2,3);
for i = 1:3
    for idir = 1:nOri
    ind_dir = intersect(ind{i}, find(dir_mat3 == dirs(idir)));
    for iI = 1:length(ind_dir)
    iind = ind_dir(iI);
    half2_noise(:,iind,1,i) = resp2(:,iind) - mean(resp2(:,ind_dir),2);
    half1_noise(:,iind,2,i) = squeeze(std(resp(:,ind_dir),[],2))./sqrt(length(ind_ori));
    end
    end
end
half_noise_corr5 = NaN(nCells1,nCells1);
half_noise_corr6 = NaN(nCells1,nCells1);
full_noise_corr3 = NaN(nCells1,nCells1);
for iCell = 1:nCells1
    for iC = 1:nCells1
    half_noise_corr5(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,1),half2_noise(iC,:,1,1)));
    half_noise_corr6(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,2),half2_noise(iC,:,1,2)));
    full_noise_corr3(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,3),half2_noise(iC,:,1,3)));
    end
end

%% SNR and Skew
trial_tc1 = reshape(npSub_tc1,[nFrames1 nTrials1 nCells1]);
trial_tc2 = reshape(npSub_tc2,[nFrames2 nTrials2 nCells2]);
trial_tc3 = reshape(npSub_tc3,[nFrames3 nTrials3 nCells3]);

off_only1 = [ ];
for iTrial = 1:nTrials1
    off_only1 = [off_only1; squeeze(trial_tc1(nOff1/2:nOff1,iTrial,:))];
end
meanf1 = mean(off_only1,1);
stdf1 = std(off_only1,[],1);
SNR1 = meanf1./stdf1;
save(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], [ref_date '_' mouse '_' run_str '_off_only.mat']), 'off_only1')
skew1 = skewness(off_only1);

off_only2 = [ ];
for iTrial = 1:nTrials2
    off_only2 = [off_only2; squeeze(trial_tc2(nOff2/2:nOff2,iTrial,:))];
end
meanf2 = mean(off_only2,1);
stdf2 = std(off_only2,[],1);
SNR2 = meanf2./stdf2;
SNR2(isnan(SNR2)) = 0;
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_off_only.mat']), 'off_only2')
skew2 = skewness(off_only2);

off_only3 = [ ];
for iTrial = 1:nTrials3
    off_only3 = [off_only3; squeeze(trial_tc3(nOff3/2:nOff3,iTrial,:))];
end
meanf3 = mean(off_only3,1);
stdf3 = std(off_only3,[],1);
SNR3 = meanf3./stdf3;
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_off_only.mat']), 'off_only3')
skew3 = skewness(off_only3);


