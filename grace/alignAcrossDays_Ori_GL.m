clear all
clear global
%% 
mouse = 'i1800';
ref_date = '210505';
day2 = '210507';
day3 = '210512';
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
align = load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_multiday_alignment.mat']));
redChImg = align.redChImg;
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
align2 = load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_multiday_alignment.mat']));
trans2 = align2.fitGeoTAf;
dfof_max2 = align2.dfmax;
dfof_max2 = dfof_max2{3};
pix2 = align2.corrmap;
pix2 = pix2{3};
redChImg2 = align2.redChImg;

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

dir_mat3 = tGratingDir3;
ori_mat3 = dir_mat3;
ori_mat3(find(dir_mat3>=180)) = ori_mat3(dir_mat3>=180)-180;
oris3 = unique(ori_mat3);
nOri3 = length(oris3);
dirs3 = unique(dir_mat3);
nDir3 = length(dirs3);

cells = intersect(cell_list2,cell_list3);
% [out,day2_cells] = intersect(cell_list2,day1_cells);
% [out,day3_cells] = intersect(cell_list3,day1_cells);
goodfit_D1 = find(oriTuning_D1.fitReliability<22.5);
goodfit_D2 = find(oriTuning_D2.fitReliability<22.5);
goodfit_D3 = find(oriTuning_D3.fitReliability<22.5);
goodfit1 = intersect(goodfit_D1,goodfit_D2);
goodfit = intersect(goodfit1,goodfit_D3);
cell_list = intersect(cells,goodfit);
nCells = length(cell_list);

%% define variables
avgResp_D1 = oriTuning_D1.avgResponseEaOri(cell_list,:);
avgResp_D2 = oriTuning_D2.avgResponseEaOri(cell_list,:);
avgResp_D3 = oriTuning_D3.avgResponseEaOri(cell_list,:);

semResp_D1 = oriTuning_D1.semResponseEaOri(cell_list,:);
semResp_D2 = oriTuning_D2.semResponseEaOri(cell_list,:);
semResp_D3 = oriTuning_D3.semResponseEaOri(cell_list,:);

vonMisesFit_D1 = oriTuning_D1.vonMisesFitAllCellsAllBoots(:,:,cell_list);
vonMisesFit_D2 = oriTuning_D2.vonMisesFitAllCellsAllBoots(:,:,cell_list);
vonMisesFit_D3 = oriTuning_D3.vonMisesFitAllCellsAllBoots(:,:,cell_list);

[maxResp_D1 prefOri_D1] = max(squeeze(oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);
[maxResp_D2 prefOri_D2] = max(squeeze(oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);
[maxResp_D3 prefOri_D3] = max(squeeze(oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,cell_list)),[],1);

% save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_prefOri.mat']),'prefOri_D1','prefOri_D2','prefOri_D3','avgResp_D1','avgResp_D2','avgResp_D3');

%% START HERE withuot frequency 
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
h1 = NaN(nDir,nCells);
p1 = NaN(nDir,nCells);
h2 = NaN(nDir,nCells);
p2 = NaN(nDir,nCells);
h3 = NaN(nDir,nCells);
p3 = NaN(nDir,nCells);
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
for i = 1:3
  for iDir = 1:nDir
    ind_dir = intersect(ind{i},find(dir_mat3 == dirs(iDir)));
    dir_resp_avg3(:,iDir,i) = squeeze(mean(resp3(cell_list,ind_dir),2)); 
    dir_resp_max3(:,iDir,i) = squeeze(max(resp3(cell_list,ind_dir),[],2));
  end
end

%% selectivity
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
% figure;
for iCell = 1:nCells
    t(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t(iCell,it);
    c_pass(iCell,it) = sum(dir_resp_avg_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S(iCell) = 1 - (2*(sum(c_pass(iCell,:),2)/100));
%     subplot(n,n2,iCell)
%     bar(t(iCell,:),c_pass(iCell,:),'r')
%     title(['cell ' num2str(iCell) ' S ' num2str(S(iCell))])
%     xlabel(['Threshold'])
%     ylabel(['Responses'])
%     print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['threshold.pdf']),'-dpdf', '-bestfit')
end
% figure;
for iCell = 1:nCells
    t2(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t2(iCell,it);
    c_pass2(iCell,it) = sum(dir_resp_avg2_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S2(iCell) = 1 - (2*(sum(c_pass2(iCell,:),2)/100));
%     subplot(n,n2,iCell)
%     bar(t2(iCell,:),c_pass2(iCell,:),'r')
%     title(['cell ' num2str(iCell) ' S ' num2str(S2(iCell))])
%     xlabel(['Threshold'])
%     ylabel(['Responses'])
%     print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], ['threshold2.pdf']),'-dpdf', '-bestfit')
end
% figure;
for iCell = 1:nCells
    t3(iCell,:) = linspace(0.0,f_max(iCell,3),100);
    for it = 1:100
        iT = t3(iCell,it);
    c_pass3(iCell,it) = sum(dir_resp_avg3_rect(iCell,:,3)>=iT)/size(dir_resp_avg_rect,2);
    end
    S3(iCell) = 1 - (2*(sum(c_pass3(iCell,:),2)/100));
%     subplot(n,n2,iCell)
%     bar(t3(iCell,:),c_pass3(iCell,:),'r')
%     title(['cell ' num2str(iCell) ' S ' num2str(S3(iCell))])
%     xlabel(['Threshold'])
%     ylabel(['Responses'])
%     print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['threshold3.pdf']),'-dpdf', '-bestfit')
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_threshold.mat']),'t','c_pass','t2','c_pass2','t3','c_pass3','S','S2','S3');

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
% figure;
% cdfplot(sum(h1,1));hold on;cdfplot(sum(h2,1));hold on;cdfplot(sum(h3,1))
% legend({'day1' 'day2' 'day3'})
%% 
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
        signal_corr_half3(iCell,1) = triu2vec(corrcoef(dir_resp_avg3(iCell,:,1),dir_resp_avg3(iCell,:,2)));
        signal_corr_day2(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,3),dir_resp_avg2(iCell,:,3)));
        signal_corr_day3(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,3),dir_resp_avg3(iCell,:,3)));
        signal_corr_day_half2(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1),dir_resp_avg2(iCell,:,1)));
        signal_corr_day_half3(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1),dir_resp_avg3(iCell,:,1)));
    end
end
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_signalCorr.mat']),'signal_corr_half1','signal_corr_half2','signal_corr_half3','signal_corr_day2','signal_corr_day3','signal_corr_day_half2','signal_corr_day_half3');

%% vonmises without frequency 
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

%% tuning curves
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
%% 
cell_stats = regionprops(mask_cell);
cell_stats = cell_stats(cell_list);
cell_stats2 = regionprops(mask_cell2);
cell_stats2 = cell_stats2(cell_list);
% cell_stats3 = regionprops(mask_cell3);
% cell_stats3 = cell_stats3(cell_list);
figure;
start = 1;
for iCell = 1:10
    width = 24; height = 24;
    xCenter = round(cell_stats(iCell).Centroid(2));
    yCenter = round(cell_stats(iCell).Centroid(1));
    xLeft(iCell) = (xCenter - width/2);
    yBottom(iCell) = (yCenter - height/2);
    xCenter2 = round(cell_stats2(iCell).Centroid(2));
    yCenter2 = round(cell_stats2(iCell).Centroid(1));
    xLeft2(iCell) = (xCenter2 - width/2);
    yBottom2(iCell) = (yCenter2 - height/2);
%     xCenter3 = round(cell_stats3(iCell).Centroid(2));
%     yCenter3 = round(cell_stats3(iCell).Centroid(1));
%     xLeft3(iCell) = (xCenter3 - width/2);
%     yBottom3(iCell) = (yCenter3 - height/2);
    cell_mask = mask_cell(xLeft(iCell):(xLeft(iCell)+width),yBottom(iCell):(height+yBottom(iCell)));
    cell_mask2 = mask_cell2(xLeft2(iCell):(xLeft2(iCell)+width),yBottom2(iCell):(height+yBottom2(iCell)));
%     cell_mask3 = mask_cell3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
%     if xLeft(iCell) > 12 && xLeft(iCell) < 488 && yBottom(iCell) > 12 && yBottom(iCell) < 772
    subplot(10,4,start);
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
    subplot(10,4,start+1);
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
    subplot(10,4,start+2);
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
    subplot(10,4,start+3);
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
%     subplot(10,6,start+4);
%     h = pix3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
%     imagesc(h)
%     pos = get(gca, 'Position');
%     pos(1) = 0.225;
%     pos(3) = 0.05;
%     set(gca, 'Position', pos)
%     axis off
%     axis square
%     hold on
%     bound = cell2mat(bwboundaries(cell_mask3(:,:,1)));
%     plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
%     title('red')
%     subplot(10,6,start+5);
%     v = dfof_max3(xLeft3(iCell):(xLeft3(iCell)+width),yBottom3(iCell):(height+yBottom3(iCell)));
%     imagesc(v)
%     pos = get(gca, 'Position');
%     pos(1) = 0.275;
%     pos(3) = 0.05;
%     set(gca, 'Position', pos)
%     axis off
%     axis square
%     hold on
%     bound = cell2mat(bwboundaries(cell_mask3(:,:,1)));
%     plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
%     title('dfof')
   
    start = start+4;
    end
% end
% print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_3day_cell_squares.pdf']),'-dpdf', '-bestfit')

%% numbered maps 
figure; imagesc(dfof_max)
hold on
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell == iCell))
        text(cell_stats(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
            'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
    else
        day1_cells(iC) = NaN;
    end
end
title('day1')
figure; imagesc(dfof_max2)
hold on
bound = cell2mat(bwboundaries(mask_cell2(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell2 == iCell))
        text(cell_stats2(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
            'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
    else
        cell_list3(iC) = NaN;
    end
end
title('day2')
figure; imagesc(dfof_max3)
hold on
bound = cell2mat(bwboundaries(mask_cell3(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if length(find(mask_cell3 == iCell))
        text(cell_stats3(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
            'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
    else
        cell_list3(iC) = NaN;
    end
title('day3')
end
%% figures
figure;
filler = zeros(size(mask_cell));
imshow(cat(3,mask_cell3,mask_cell,filler))

figure;bar([max(maxResp_D1) mean(maxResp_D1); max(maxResp_D2) mean(maxResp_D2); max(maxResp_D3) mean(maxResp_D3)])
legend({'max resp for all cells' 'mean max resp for all cells'})
xlabel('Day','FontSize',20);ylabel('dfof','FontSize',20)
title([mouse])
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], [day3 '_' mouse '_' run_str3 '_max_dfof_bar.pdf']),'-dpdf', '-bestfit')

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
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str3], ['signal_corr_scat.pdf']),'-dpdf', '-bestfit')

figure;scatter(maxResp_D1,maxResp_D2,sz,'filled')
hold on;scatter(maxResp_D1,maxResp_D3,sz,'filled')
ax = gca;
ax.FontSize = 20;
% xlim([0 1])
% ylim([0 1])
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


% without SF
indc{1} = 'D1 half1';
indc{2} = 'D1 half2';
indc{3} = 'day1';
indc{4} = 'day2';
in{1} = signal_corr_half1;
in{2} = signal_corr_day_half;
iCell = 2;
figure;
suptitle([mouse ' Cell ' num2str(iCell)])
start = 1;
for iCell = 2
    for iTrial = 1
    subplot(2,2,start)
    imagesc(dir_resp_avg(iCell,:,iTrial))
    title([indc{iTrial} ' ' round(num2str(in{iTrial}(iCell)))])
    subplot(2,2,start+1)
    imagesc(dir_resp_avg(iCell,:,iTrial+1))
    title([indc{iTrial+1}])
    subplot(2,2,start+2)
    imagesc(dir_resp_avg(iCell,:,iTrial))
    title([indc{iTrial+2} ' ' round(num2str(in{iTrial+1}(iCell)))])
    subplot(2,2,start+3)
    imagesc(dir_resp_avg2(iCell,:,iTrial))
    title([indc{iTrial+3}])
    end
end

%% histograms
figure;histogram(corrcoef(avgResp_D1(:,:)'),20)

same_oriD1 = zeros(nCells1,nCells1);
for icell = 1:length(goodfit)
    iCell = goodfit(icell);
     for ic = 1:length(goodfit)
        iC = goodfit(ic);
        if iCell~= iC
        if abs(prefOri_D1(iCell)-prefOri_D1(iC))<11
            same_oriD1(iCell,iC) = triu2vec(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)));
%         if prefOri_D1(iCell)>prefOri_D1(iC)+11||prefOri_D1(iCell)<prefOri_D1(iC)-11
%             hold on
%             histogram(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)))
        else
            same_oriD1(iCell,iC) = NaN;
        end
        end
     end
end

diff_oriD1 = zeros(nCells1,nCells1);
for icell = 1:length(goodfit)
    iCell = goodfit(icell);
     for ic = 1:length(goodfit)
        iC = goodfit(ic);
        if iCell~= iC
        if abs(prefOri_D1(iCell)-prefOri_D1(iC))>11
            diff_oriD1(iCell,iC) = triu2vec(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)));
%         if prefOri_D1(iCell)>prefOri_D1(iC)+11||prefOri_D1(iCell)<prefOri_D1(iC)-11
%             hold on
%             histogram(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)))
        else
            diff_oriD1(iCell,iC) = NaN;
        end
        end
     end
end

ortho_oriD1 = zeros(nCells1,nCells1);
for icell = 1:length(goodfit)
    iCell = goodfit(icell);
     for ic = 1:length(goodfit)
        iC = goodfit(ic);
        if iCell~= iC
        if abs(prefOri_D1(iCell)-prefOri_D1(iC))>79&&abs(prefOri_D1(iCell)-prefOri_D1(iC))<101
            ortho_oriD1(iCell,iC) = triu2vec(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)));
%         if prefOri_D1(iCell)>prefOri_D1(iC)+11||prefOri_D1(iCell)<prefOri_D1(iC)-11
%             hold on
%             histogram(corrcoef(avgResp_D1(iCell,:),avgResp_D1(iC,:)))
        else
            ortho_oriD1(iCell,iC) = NaN;
        end
        end
     end
end

same_oriD1 = same_oriD1(~isnan(same_oriD1));
diff_oriD1 = diff_oriD1(~isnan(diff_oriD1));
ortho_oriD1 = ortho_oriD1(~isnan(ortho_oriD1));
same_oriD1(same_oriD1==0)=[];
diff_oriD1(diff_oriD1==0)=[];
ortho_oriD1(ortho_oriD1==0)=[];

figure;
histogram(same_oriD1,'BinWidth',0.1,'FaceColor',[0 0.3906 0])
hold on
histogram(diff_oriD1,'BinWidth',0.1,'FaceColor',[0.5430 0 0.5430])
hold on
histogram(ortho_oriD1,'BinWidth',0.1,'FaceColor',[1 0 0])
xlabel('signal correlation')
ylabel('cell pairs')
legend({'same pref ori','diff pref ori','ortho pref ori'},'northeast')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['pref_ori_hist.pdf']),'-dpdf', '-bestfit')


% comparing max dfof between days
ori_resp_max = zeros(nCells2,1);
ori_resp_max2 = zeros(nCells2,1);
for iCell = 1:nCells2
    ori_resp_max(iCell,1) = max(dir_resp_avg(iCell,:,3));
    ori_resp_max2(iCell,1) = max(dir_resp_avg2(iCell,:,3));
end
figure;scatter(ori_resp_max,ori_resp_max2)
% xlim([0 3]);ylim([0 3])
refline(1,0)
xlabel('dfof of pref ori D1')
ylabel('dfof of pref ori D2')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['max_dfof_scatter.pdf']),'-dpdf', '-bestfit')
figure;histogram(ori_resp_max./ori_resp_max2,'BinWidth',.5,'BinLimits',[0 20])
hold on;histogram(ori_resp_max-ori_resp_max2,'BinWidth',.5,'BinLimits',[0 20])
legend({'max dfof D1/D2','max dfof D1 - D2'},'northeast')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str2], ['dfof_hist.pdf']),'-dpdf', '-bestfit')

% diff in signal corr between trial halves
[B,I] = sort(goodfit_PO1);
D1half1 = corrcoef(dir_resp_avg(I,:,1,1)');
D1half2 = corrcoef(dir_resp_avg(I,:,1,2)');
half_diff = abs(D1half1-D1half2);
color_axis_limit = [0 2];
figure;
title('abs difference in signal cor between randomized trial halves')
colormap(brewermap([],'*Blues'))
imagesc(half_diff)
colorbar
clim(color_axis_limit)
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['same_day_colormap.pdf']),'-dpdf', '-bestfit')

figure;
histogram(D1half1,'BinWidth',0.1,'FaceColor',[0 0.543 0.543])
hold on
histogram(D1half2,'BinWidth',0.1,'FaceColor',[1 0.8398 0])
legend({'half 1','half 2'},'northeast')
xlabel('signal correlation')
ylabel('cell pairs')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['half_hist.pdf']),'-dpdf', '-bestfit')

pair_list = cell(length(goodfit),1);
pair_list2 = cell(length(goodfit),1);
for iCell = goodfit
    valid_cells = [];
    valid_cells2 = [];
    for iC = goodfit
        half1 = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,1)',dir_resp_avg(iC,:,1,1)'));
        half2 = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,2)',dir_resp_avg(iC,:,1,2)'));
        if abs(half1-half2)>0.7
%              if isempty(find(pair_list==iCell)) && isempty(find(pair_list==iC))
            valid_cells = [valid_cells,iC];
        end
        if abs(half1-half2)<0.2
            valid_cells2 = [valid_cells2,iC];
        end
    end
    pair_list{iCell} = valid_cells;
    pair_list2{iCell} = valid_cells2;
end
  
% are corr coefs between halves higher than between cells
nCells = nCells1;
half_corr = NaN(nCells,1);
cell_corr1 = NaN(nCells,nCells);
cell_corr2 = NaN(nCells,nCells);
avg_corr = NaN(nCells,nCells);
avg_corr2 = NaN(nCells,nCells);
two_day_corr = NaN(nCells,nCells);
for iCell = 1:nCells
        half_corr(iCell,1) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,1)',dir_resp_avg(iCell,:,1,2)'));
        for iC = 1:nCells
        cell_corr1(iCell,iC) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,1)',dir_resp_avg(iC,:,1,1)'));
        cell_corr2(iCell,iC) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,2)',dir_resp_avg(iC,:,1,2)'));
        avg_corr(iCell,iC) = triu2vec(corrcoef(dir_resp_avg(iCell,:,1,3),dir_resp_avg(iC,:,1,3)));
        avg_corr2(iCell,iC) = triu2vec(corrcoef(dir_resp_avg2(iCell,:,1,3)',dir_resp_avg2(iC,:,1,3)'));
        end
end
half_dif = abs(cell_corr1(cell_list,cell_list)-cell_corr2(cell_list,cell_list));
half_dif = half_dif(~isnan(half_dif));
day_dif = abs(avg_corr(goodfit,goodfit)-avg_corr2(goodfit,goodfit));
day_dif = day_dif(~isnan(day_dif));
figure;histogram(half_dif,'BinWidth',0.1,'FaceColor',[0 0 .5430]);hold on
histogram(day_dif,'BinWidth',0.1,'FaceColor',[0.2500 0.8750 0.8125])
title('differences in signal correlation for reliably fit cell pairs')
legend({'between halves','between days'},'Location','northeast')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['halves_days_hist.pdf']),'-dpdf', '-bestfit')

% pref ori vs signal correlation
diff_ori = zeros(length(cell_list),length(cell_list));
diff_ori2 = zeros(length(cell_list),length(cell_list));
fit = zeros(nCells1,nCells1);
for icell = 1:length(cell_list)
    iCell = cell_list(icell);
    for ic = 1:length(cell_list)
       iC = cell_list(ic); 
        diff_ori(iCell,iC) = abs(prefOri_D1(iC) - prefOri_D1(iCell));
        diff_ori2(iCell,iC) = abs(prefOri_D2(iC)- prefOri_D2(iCell));
%         fit(iCell,iC) = abs(fitR1(iC) - fitR1(iCell));
    end
end
diff_ori1 = reshape(diff_ori(goodfit_D1,goodfit_D1),[],1);
diff_ori2 = reshape(diff_ori2(goodfit_D2,goodfit_D2),[],1);
corr = reshape(avg_corr(goodfit_D1,goodfit_D1),[],1);
fit = reshape(fit,[],1);
cell_corr_diff = reshape(abs(cell_corr1-cell_corr2),[],1);
figure;scatter(diff_ori1,corr)
xlabel('delta pref ori')
ylabel('signal corr')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['ori_vs_corr.pdf']),'-dpdf', '-bestfit')
figure;scatter(fit,corr)
ylabel('delta reliability of tuning')
xlabel('signal corr')
yticks([0:1:20])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['fit_vs_corr.pdf']),'-dpdf', '-bestfit')

suptitle(['Cell 123, Cell N, Half1, Half2, all trials'])
for iCell = 5
    figure;
    [n n2] = subplotn(size(pair_list{iCell},2));
    for e = 1:8
        sample = pair_list{iCell};
        iC = sample(e);
%         x = length(e);
        subplot(n,n2,e) 
        y = [half_corr(iCell,:),half_corr(iC,:),cell_corr1(iCell,iC),cell_corr2(iCell,iC),avg_corr(iCell,iC),avg_corr2(iCell,iC)];
        bar(y)
    end
end

suptitle(['cell' num2str(iCell) ' pairs with delta sig corr > 0.7 btwn halves'])
start = 1;
for iCell = 5
    figure;
    [n n2] = subplotn(size(pair_list{iCell},2));
    for e = 1:8
        sample = pair_list{iCell};
        iC = sample(e);
for i = 1:2
    subplot(n,n2,e)
    title(num2str(iC))
    errorbar(oris, dir_resp_avg(iCell,:,1,i), dir_resp_avg(iCell,:,2,i), '-o')
    hold on
    errorbar(oris, dir_resp_avg(iC,:,1,i), dir_resp_avg(iC,:,2,i), '-o')
    hold on
end
    start = start+1;
    end
%     legend({'C1H1','C2H1','C1H2','C2H2','overlap'},'Location','northeast')
end

% tuning curves
figure;
start = 1;
for iC = 1:length(goodfit)
    iCell = goodfit(iC);
  for i = 1:2
    if start>25
      figure;
      start = 1;
    end
    subplot(5,5,start)
    errorbar(oris, dir_resp_avg(iCell,:,1,i), dir_resp_avg(iCell,:,2,i), '-o')
    hold on
  end
  start = start+1;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['half_tuning.pdf']),'-dpdf', '-bestfit')


% avg dfof response for cell pairs in each trial half
figure;
suptitle('Cell1 Half1, Cell1 Half2, Cell2 Half1, Cell2 Half2, similar ori')
start = 1;
for icell = 1:length(goodfit)
    iCell = rand(goodfit(icell));
     for ic = 1:length(goodfit)
        iC = rand(goodfit(ic));
        if iCell~= iC
            if prefOri_D1(iCell)<prefOri_D1(iC)+11&&prefOri_D1(iCell)>prefOri_D1(iC)-11
            subplot(4,4,start)
            plot(1:8,dir_resp_avg(iCell,:,1,1),'color',[1 0 0])
            hold on
            plot(1:8,dir_resp_avg(iCell,:,1,2),'color',[1 0.5469 0])
            hold on
            plot(1:8,dir_resp_avg(iC,:,1,1),'color',[1 0.8398 0])
            hold on
            plot(1:8,dir_resp_avg(iC,:,1,2),'color',[0.1328 0.5430 0.1328])
            title([num2str(iCell) 'and' num2str(iC)])
            start = start+1;
            axis square
            end
        end
     end
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['same_ori_AR.pdf']),'-dpdf', '-bestfit')

figure;
suptitle('Cell1 Half1, Cell1 Half2, Cell2 Half1, Cell2 Half2, diff ori')
start = 1;
for icell = 1:length(goodfit)
    iCell = goodfit(icell);
     for ic = 1:length(goodfit)
        iC = goodfit(ic);
        if iCell~=iC
        if prefOri_D1(iCell)>prefOri_D1(iC)+11||prefOri_D1(iCell)<prefOri_D1(iC)-11
            subplot(4,4,start)
            plot(1:8,dir_resp_avg(iCell,:,1,1),'color',[0.6016 0.8008 0.1953])
            hold on
            plot(1:8,dir_resp_avg(iCell,:,1,2),'color',[0.2500 0.8750 0.8125])
            hold on
            plot(1:8,dir_resp_avg(iC,:,1,1),'color',[0.9297 0.5078 0.9297])
            hold on
            plot(1:8,dir_resp_avg(iC,:,1,2),'color',[0.4648 0.5313 0.5977])
            title([num2str(iCell) 'and' num2str(iC)])
            start = start+1;
            axis square
        end
        end
     end
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['diff_ori_AR.pdf']),'-dpdf', '-bestfit')

% scatter
nCells = length(goodfit);
diffOri1 = zeros(nCells,nCells);
diffOri2 = zeros(nCells,nCells);
diffOri3 = zeros(nCells,nCells);
ddeltaOri = zeros(nCells,nCells);
diffSigCor = zeros(nCells,nCells);
for icell = 1:length(goodfit)
    iCell = (goodfit(icell));
    for ic = 1:length(goodfit)
        iC = goodfit(ic);
    diffOri1(iC,iCell) = prefOri_D1(iC) - prefOri_D1(iCell);
    diffOri2(iC,iCell) = prefOri_D2(iC) - prefOri_D2(iCell);
    diffOri3(iC,iCell) = prefOri_D3(iC) - prefOri_D3(iCell);
    ddeltaOri(iC,iCell) = diffOri1(iC,iCell)-diffOri2(iC,iCell);
    D1sigcor = triu2vec(corrcoef(avgResp_D1(iC,:),avgResp_D1(iCell,:)));
    D2sigcor = triu2vec(corrcoef(avgResp_D2(iC,:),avgResp_D2(iCell,:)));
    diffSigCor(iC,iCell) = D1sigcor - D2sigcor;
    end
end
% diffOri1 = tril(diffOri1);
% diffOri1(diffOri1==0)=[];
% diffOri2 = tril(diffOri2);
% diffOri2(diffOri2==0)=[];
% diffOri3 = tril(diffOri3);
% diffOri3(diffOri3==0)=[];
% PO1 = find(diffOri1<22.5&diffOri1>-22.5);
% PO2 = find(diffOri2<22.5&diffOri2>-22.5);
% PO3 = find(diffOri3<22.5&diffOri3>-22.5);
% prefOri2 = intersect(PO1,PO2);
% prefOri3 = intersect(PO1,PO3);
deltPrefOri = tril(ddeltaOri);
deltPrefOri(deltPrefOri==0)=[];
deltSigCor = tril(diffSigCor);
deltSigCor(deltSigCor==0)=[];
figure;scatter(deltPrefOri,deltSigCor(1:344))
xlabel('delta pref ori D1 - delta pref ori D2','FontSize',20)
ylabel('difference in D1 vs D2 signal corr','FontSize',20)
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['signal_corr_scatter.pdf']),'-dpdf', '-bestfit')


% all cell pairs
D1 = tril(triu2vec(corrcoef(avgResp_D1,avgResp_D1)));
D2 = tril(triu2vec(corrcoef(avgResp_D2,avgResp_D2)));
D3 = tril(triu2vec(corrcoef(avgResp_D3,avgResp_D3)));
D1(D1==0)=[];
D2(D2==0)=[];
D3(D3==0)=[];
b = [D1(1:100)' D2(1:100)' D3(1:100)'];
x = [1 2 3];
figure;
plot(x,b)

D12same = find(D2<D1+0.2 & D2>D1-0.2);
D12large = find(D2>D1+0.2);
D12small = find(D2<D1-0.2);
D13same = find(D3<D1+0.2 & D3>D1-0.2);
D13large = find(D3>D1+0.2);
D13small = find(D3<D1-0.2);
overlapSame = intersect(D12same,D13same);
overlapLarge = intersect(D12large,D13large);
overlapSmall = intersect(D12small,D13small);
b = [length(D12same) length(D13same) length(overlapSame);length(D12large) length(D13large) length(overlapLarge);length(D12small) length(D13small) length(overlapSmall)];
figure;
c = bar(b);
cellnames = {'within 0.2','greater than 0.2','less than 0.2'};
set(gca,'xticklabel',cellnames,'FontSize',20)
ylabel('n cell pairs')
legend({'Days 1 & 2', 'Days 1 & 3','overlap'},'Location','northeast')
c(1).FaceColor = [0 0.3906 0];
c(2).FaceColor = [0.1797 0.5430 0.3398];
c(3).FaceColor = [0.5625 0.9297 0.5625];
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['signal_corr_flux.pdf']),'-dpdf', '-bestfit')


% noise correlation
% ind is the trials where the orientations are the same
nTrials1 = 160;
dir_mat = tGratingDir1;
dir_mat = dir_mat(1:160);
ori_mat = dir_mat(1:160);
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
nOri = length(oris);
ind = cell(1,3);
ind{1} = randperm(nTrials1,nTrials1/2);
ind{2} = setdiff(1:nTrials1, ind{1});
ind{3} = 1:nTrials1;
half1_noise = zeros(nCells1,nOri,2,3);
for i = 1:3
    for idir = 1:nOri
    ind_dir = intersect(ind{i}, find(ori_mat == oris(idir)));
    for iI = 1:length(ind_dir)
    iind = ind_dir(iI);
    half1_noise(:,iind,1,i) = resp(:,iind) - mean(resp(:,ind_dir),2);
%     half1_noise(:,iind,2,i) = squeeze(std(resp(:,ind_dir),[],2))./sqrt(length(ind_ori));
    end
    end
end
half_noise_corr1 = NaN(nCells1,nCells1);
half_noise_corr2 = NaN(nCells1,nCells1);
full_noise_corr1 = NaN(nCells1,nCells1);
dfof_pair_diff = NaN(nCells1,nCells1);
for iCell = 1:nCells1
    for iC = 1:nCells1
    half_noise_corr1(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,1),half1_noise(iC,:,1,1)));
    half_noise_corr2(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,2),half1_noise(iC,:,1,2)));
    full_noise_corr1(iCell,iC) = triu2vec(corrcoef(half1_noise(iCell,:,1,3),half1_noise(iC,:,1,3)));
    dfof_pair_diff(iCell,iC) = ori_resp_max(iCell) - ori_resp_max(iC);
    end
end
half_diff_noise = abs(half_noise_corr2 - half_noise_corr1);
half_diff_noise1 = reshape(half_diff_noise(goodfit_D1,goodfit_D1),[],1);
half_diff_signal1 = half_dif;
figure;histogram(half_diff_signal1,'BinWidth',0.1);hold on;histogram(half_diff_noise,'BinWidth',0.1)
legend({'delta signal corr','delta noise corr'},'Location','northeast')
title('difference in coorelations between halves')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['halfdif_sig_vs_noise.pdf']),'-dpdf', '-bestfit')
noise_corr1 = reshape(full_noise_corr1(goodfit,goodfit),[],1);
figure;scatter(diff_ori1,noise_corr1)
ylabel('noise correlation')
xlabel('delta pref ori')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['ori_vs_noise.pdf']),'-dpdf', '-bestfit')
figure;scatter(corr,noise_corr1)
xlabel('signal correlation')
ylabel('noise correlation')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['signal_vs_noise.pdf']),'-dpdf', '-bestfit')

dir_mat = tGratingDir2;
dir_mat = dir_mat(1:160);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
nOri = length(oris);
ind{1} = randperm(nTrials2,nTrials2/2);
ind{2} = setdiff(1:nTrials2, ind{1});
ind{3} = 1:nTrials2;
half2_noise = zeros(nCells1,nOri,2,3);
for i = 1:3
    for idir = 1:nOri
    ind_dir = intersect(ind{i}, find(ori_mat == oris(idir)));
    for iI = 1:length(ind_dir)
    iind = ind_dir(iI);
    half2_noise(:,iind,1,i) = resp2(:,iind) - mean(resp2(:,ind_dir),2);
%     half1_noise(:,iind,2,i) = squeeze(std(resp(:,ind_dir),[],2))./sqrt(length(ind_ori));
    end
    end
end
half_noise_corr3 = NaN(nCells1,nCells1);
half_noise_corr4 = NaN(nCells1,nCells1);
full_noise_corr2 = NaN(nCells1,nCells1);
dfof_pair_diff2 = NaN(nCells1,nCells1);
for iCell = 1:nCells1
    for iC = 1:nCells1
    half_noise_corr3(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,1),half2_noise(iC,:,1,1)));
    half_noise_corr4(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,2),half2_noise(iC,:,1,2)));
    full_noise_corr2(iCell,iC) = triu2vec(corrcoef(half2_noise(iCell,:,1,3),half2_noise(iC,:,1,3)));
    dfof_pair_diff2(iCell,iC) = ori_resp_max2(iCell) - ori_resp_max2(iC);
    end
end
noise_corr2 = reshape(full_noise_corr2(goodfit_D2,goodfit_D2),[],1);
figure;scatter(diff_ori2,noise_corr2)
ylabel('noise correlation')
xlabel('delta pref ori')
    
half_noise_diff = abs(half_noise_corr1(goodfit,goodfit)-half_noise_corr2(goodfit,goodfit));
day_noise_diff = tril(full_noise_corr1(goodfit,goodfit)-full_noise_corr2(goodfit,goodfit));
day_noise_diff(day_noise_diff==0)=[];
day_dif_noise = reshape(day_noise_diff,[],1);
figure;histogram(half_noise_diff,'BinWidth',0.1,'BinLimits',[0 2],'FaceColor',[0 0.3 0.6])
hold on;histogram(day_noise_diff,'BinWidth',0.1,'BinLimits',[0 2],'FaceColor',[0.6016 0.8008 0.1953])
title('differences in noise correlation for all cells')
legend({'between halves','between days'},'Location','northeast')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['noise_halves_days_hist.pdf']),'-dpdf', '-bestfit')

day_dif_signal = reshape(day_dif,[],1);
figure;histogram(day_dif_signal,'BinWidth',0.1);hold on;histogram(day_dif_noise,'BinWidth',0.1)
legend({'delta signal corr','delta noise corr'},'Location','northeast')
title('difference in correlations between days')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['daydif_sig_vs_noise.pdf']),'-dpdf', '-bestfit')

dfof_day_diff = tril(dfof_pair_diff(cell_list,cell_list) - dfof_pair_diff2(cell_list,cell_list));
dfof_day_diff(dfof_day_diff==0)=[];
figure;histogram(day_noise_diff,'BinWidth',0.1)
hold on;histogram(dfof_day_diff,'BinWidth',0.1)
legend({'delta noise correlation','delta dfof between pairs'},'Location','northeast')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['noise_vs_dfof_days.pdf']),'-dpdf', '-bestfit')

%% original plots 
figure; 
subplot(1,2,1)
scatter(prefOri_D1,prefOri_D2,'filled','MarkerFaceColor',[0 0.5430 0.5430])
hold on
scatter(prefOri_D1(cells_all2),prefOri_D2(cells_all2),'filled','MarkerFaceColor',[1 0.8398 0])
% scatter(maxResp_D1(goodfit),maxResp_D2(goodfit))
axis square
% xlim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
% ylim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
refline(1,0)
set(gca,'FontSize',14)
xlabel('D1 Preferred Ori (deg)')
ylabel('D2 Preferred Ori (deg)')
subplot(1,2,2)
scatter(prefOri_D1,prefOri_D3,'filled','MarkerFaceColor',[0 0.5430 0.5430])
hold on
scatter(prefOri_D1(cells_all3),prefOri_D3(cells_all3),'filled','MarkerFaceColor',[1 0.8398 0])
% scatter(maxResp_D1(goodfit),maxResp_D2(goodfit))
axis square
% xlim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
% ylim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
refline(1,0)
set(gca,'FontSize',14)
xlabel('D1 Preferred Ori (deg)')
ylabel('D3 Preferred Ori (deg)')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' run_str], ['prefOriScat.pdf']),'-dpdf', '-bestfit')


prefOri_diff2 = abs(prefOri_D1-prefOri_D2);
prefOri_diff2(find(prefOri_diff2>90)) = 180-prefOri_diff2(find(prefOri_diff2>90));
% subplot(2,2,2)
hist(prefOri_diff2)
xlabel('D1 vs D2 Pref ori')
ylabel('Number of cells')

subplot(2,2,3)
scatter(maxResp_D1,prefOri_diff2)
xlabel('Day 1 peak dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

maxResp_diff = abs((maxResp_D1-maxResp_D2)./(maxResp_D1+maxResp_D2));
subplot(2,2,4)
scatter(maxResp_diff,prefOri_diff2)
xlabel('D1 vs D2 peak dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

suptitle([mouse ' ' ref_date ' vs ' day2 '- n = ' num2str(nCells1) ' cells'])
print(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_compAcrossDaysTuning.pdf']),'-dpdf', '-bestfit')
save(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str2], [day2 '_' mouse '_' run_str2 '_compAcrossDaysTuning.mat']), 'goodfit', 'prefOri_D1', 'maxResp_D1', 'prefOri_D2', 'maxResp_D2', 'prefOri_diff2', 'maxResp_diff');

figure; 
% subplot(2,2,1)
scatter(prefOri_D1,prefOri_D3,'filled','MarkerFaceColor',[0 0.5430 0.5430])
hold on
scatter(prefOri_D1(cells_all3),prefOri_D3(cells_all3),'filled','MarkerFaceColor',[1 0.8398 0])
% scatter(maxResp_D1(goodfit),maxResp_D2(goodfit))
axis square
% xlim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
% ylim([0 max(max(maxResp_D1(goodfit),maxResp_D2(goodfit)))+0.05])
refline(1,0)
set(gca,'FontSize',16)
xlabel('D1 Preferred Ori (deg)')
ylabel('D3 Preferred Ori (deg)')

prefOri_diff3 = abs(prefOri_D1-prefOri_D3);
prefOri_diff3(find(prefOri_diff3>90)) = 180-prefOri_diff3(find(prefOri_diff3>90));
subplot(2,2,2)
hist(prefOri_diff3)
xlabel('D1 vs D3 Pref ori')
ylabel('Number of cells')

subplot(2,2,3)
scatter(maxResp_D1,prefOri_diff3)
xlabel('Day 1 peak dF/F')
ylabel('D1 vs D3 Pref ori')
axis square

maxResp_diff = abs((maxResp_D1-maxResp_D3)./(maxResp_D1+maxResp_D3));
subplot(2,2,4)
scatter(maxResp_diff,prefOri_diff3)
xlabel('D1 vs D3 peak dF/F')
ylabel('D1 vs D3 Pref ori')
axis square

suptitle([mouse ' ' ref_date ' vs ' day3 '- n = ' num2str(nCells1) ' cells'])
print(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_compAcrossDaysTuning.pdf']),'-dpdf', '-bestfit')
save(fullfile(fnout, [day3 '_' mouse], [day3 '_' mouse '_' run_str], [day3 '_' mouse '_' run_str '_compAcrossDaysTuning.mat']), 'goodfit', 'prefOri_D1', 'maxResp_D1', 'prefOri_D2', 'maxResp_D2', 'prefOri_diff3', 'maxResp_diff');

%% Identify cells that are significantly responsive to at least 1 direction

% Define analysis windows
resp_wind1 = nOff1+1:nOff1+nOn1;
base_wind1 = 1+nOff1-nOn1:nOff1;

resp_wind2 = nOff2+1:nOff2+nOn2;
base_wind2 = 1+nOff2-nOn2:nOff2;

resp_wind3 = nOff3+1:nOff3+nOn3;
base_wind3 = 1+nOff3-nOn3:nOff3;
% 
% resp_wind4 = nOff4+1:nOff4+nOn4;
% base_wind4 = 1+nOff4-nOn4:nOff4;

% Create two matrices
dfof_resp1 = squeeze(mean(trial_dfof1(resp_wind1,:,:),1));
dfof_base1 = squeeze(mean(trial_dfof1(base_wind1,:,:),1));
dfof_subtract1 = dfof_resp1 - dfof_base1;

dfof_resp2 = squeeze(mean(trial_dfof2(resp_wind2,:,:),1));
dfof_base2 = squeeze(mean(trial_dfof2(base_wind2,:,:),1));
dfof_subtract2 = dfof_resp2 - dfof_base2;

dfof_resp3 = squeeze(mean(trial_dfof3(resp_wind3,:,:),1));
dfof_base3 = squeeze(mean(trial_dfof3(base_wind3,:,:),1));
dfof_subtract3 = dfof_resp3 - dfof_base3;

% dfof_resp4 = squeeze(mean(trial_dfof4(resp_wind4,:,:),1));
% dfof_base4 = squeeze(mean(trial_dfof4(base_wind4,:,:),1));
% dfof_subtract4 = dfof_resp4 - dfof_base4;

% ttest for response significance
% nCells = intersect(cells_all2, cells_all3);
% nCells1 = length(nCells);
% nCells2 = nCells1;
% nCells3 = nCells1;
% day 1
dfof_dir1 = zeros(nDir1, nCells1, 2);
h1 = zeros(nDir1, nCells2);
p1 = zeros(nDir1, nCells2);
for idir = 1:nDir1
    ind = find(tGratingDir1==dirs(idir));
    x1 = dfof_base(ind,:);
    y1 = dfof_resp(ind,:);
    [h1(idir,:),p1(idir,:)] = ttest(x1,y1,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
%     dfof_dir1(idir,:,1) = mean(y1-x1,1);
%     dfof_dir1(idir,:,2) = std(y1-x1,[],1)./sqrt(length(ind));
end
h1sum = zeros(1, nCells1);
h1over1 = zeros(1, nCells1);
h1sum = sum(h1(:,:));
h1dirsum = sum(h1);
for iCell = 1:nCells1
    h1over1(:,iCell) = h1sum(:,iCell) >= 1;
end
sig_resp_cells1 = sum(h1over1(:,:));

h_1 = h1';
h1sum_dirs = zeros(1,nDir1);
h1sum_dirs = sum(h_1(:,:));

% day 2
dfof_dir1 = zeros(nDir2, nCells2, 2);
h2 = zeros(nDir2, nCells2);
p2 = zeros(nDir2, nCells2);
for idir = 1:nDir2
    ind = find(tGratingDir2==dirs2(idir));
    x2 = dfof_base2(ind,:);
    y2 = dfof_resp2(ind,:);
    [h2(idir,:),p2(idir,:)] = ttest(x2,y2,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir2(idir,:,1) = mean(y2-x2,1);
    dfof_dir2(idir,:,2) = std(y2-x2,[],1)./sqrt(length(ind));
end
h2sum = zeros(1, nCells2);
h2over1 = zeros(1, nCells2);
h2sum = sum(h2(:,:));
h2dirsum = sum(h2);
for iCell = 1:nCells2
    h2over1(:,iCell) = h2sum(:,iCell) >= 1;
end
sig_resp_cells2 = sum(h2over1(:,:));

h_2 = h2';
h2sum_dirs = zeros(1,nDir2);
h2sum_dirs = sum(h_2(:,:));

% day 3
dfof_dir3 = zeros(nDir3, nCells3, 2);
h3 = zeros(nDir3, nCells3);
p3 = zeros(nDir3, nCells3);
for idir = 1:nDir3
    ind = find(tGratingDir3==dirs3(idir));
    x3 = dfof_base3(ind,:);
    y3 = dfof_resp3(ind,:);
    [h3(idir,:),p3(idir,:)] = ttest(x3,y3,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
    dfof_dir3(idir,:,1) = mean(y3-x3,1);
    dfof_dir3(idir,:,2) = std(y3-x3,[],1)./sqrt(length(ind));
end
h3sum = zeros(1, nCells3);
h3over1 = zeros(1, nCells3);
h3sum = sum(h3(:,:));
h3dirsum = sum(h3);
for iCell = 1:nCells3
    h3over1(:,iCell) = h3sum(:,iCell) >= 1;
end
sig_resp_cells3 = sum(h3over1(:,:));

h_3 = h3';
h3sum_dirs = zeros(1,nDir3);
h3sum_dirs = sum(h_3(:,:));

% day 4
% dfof_dir4 = zeros(nDir4, nCells4, 2);
% h4 = zeros(nDir4, nCells4);
% p4 = zeros(nDir4, nCells4);
% for idir = 1:nDir4
%     ind = find(tGratingDir4==dirs4(idir));
%     x4 = dfof_base4(ind,:);
%     y4 = dfof_resp4(ind,:);
%     [h4(idir,:),p4(idir,:)] = ttest(x4,y4,'dim',1,'Alpha',0.05./(nDir1-1),'tail','left');
%     dfof_dir4(idir,:,1) = mean(y4-x4,1);
%     dfof_dir4(idir,:,2) = std(y4-x4,[],1)./sqrt(length(ind));
% end
% h4sum = zeros(1, nCells4);
% h4over1 = zeros(1, nCells4);
% h4sum = sum(h4(:,:));
% h4dirsum = sum(h4);
% for iCell = 1:nCells4
%     h4over1(:,iCell) = h4sum(:,iCell) >= 1;
% end
% sig_resp_cells4 = sum(h4over1(:,:));
% 
% h_4 = h4';
% h4sum_dirs = zeros(1,nDir4);
% h4sum_dirs = sum(h_4(:,:));

intersect1 = intersect(sig_resp_cells1,sig_resp_cells2);
intersect2 = intersect(intersect1, sig_resp_cells3);
% all_sig_cells = intersect(intersect2, sig_resp_cells4);

%% cell classifying figs 
figure; 
x = [1 2 3];
y = [sig_resp_cells1 sig_resp_cells2 sig_resp_cells3];
cellnames = {'Day 1'; 'Day 2'; 'Day 3'};
bar(x,y);
set(gca,'xticklabel',cellnames)
ylim([0 nCells1])
title('Significantly Responsive Cells')
ylabel('Number of cells')
% imagesc(h');
% xlabel('Direction')
% ylabel('Cells')
% figAxForm
if exist(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_cells'],'-dpdf') 

else mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays'])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_cells'],'-dpdf') 
end

figure;
x = [1 2 3];
vals = [sig_resp_cells1 gdft_D1; sig_resp_cells2 gdft_D2; sig_resp_cells3 gdft_D3]; 
cellnames = {'Day 1'; 'Day 2'; 'Day 3'};
b11 = bar(x,vals);
b11(1).FaceColor = [0.6953 0.1328 0.1328];
b11(2).FaceColor = [0.9102 0.5859 0.4766];
set(gca,'xticklabel',cellnames,'FontSize',16)
% title('Significantly Responsive Cells vs Goodfit Cells')
ylabel('Number of cells')
legend({'Sig Cells', 'Reliably Fit Cells'},'Location','northeast')
ylim([0 55])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_vs_goodfit'],'-dpdf') 

figure;
b = [h1dirsum; h2dirsum; h3dirsum]';
b1 = bar(b,'stacked');
b1(1).FaceColor = [0.0977 0.0977 0.5];
b1(2).FaceColor = [0.2334 0.4678 0.700];
b1(3).FaceColor = [0.6558 0.8238 0.8984];
b1(4).FaceColor = [0.9375 0.9708 1.0000];
legend({'Day 1', 'Day 2', 'Day 3'},'Location','northeast')
xlabel('Cell #')
ylabel('nDirs')
ylim([0 max(b)+2])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\nDirs_eachCell'],'-dpdf') 

figure;
b1 = bar(h1dirsum);
b1(1).FaceColor = [0.2334 0.4678 0.700];
xlabel('Cell #')
ylabel('nDirs')
legend({'Day 1', 'Day 2', 'Day 3'},'Location','northeast')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\nDirs_eachCellD1'],'-dpdf') 

figure;
y = [h1sum_dirs; h2sum_dirs; h3sum_dirs]';
y1 = bar(y,'stacked');
y1(1).FaceColor = [0.0977 0.0977 0.5];
y1(2).FaceColor = [0.2334 0.4678 0.700];
y1(3).FaceColor = [0.6558 0.8238 0.8984];
% y1(4).FaceColor = [0.9375 0.9708 1.0000];
% cellnames = {'0';'22.5';'45';'67.5';'90';'112.5';'135';'157.5';'180';'202.5';'225';'247.5';'270';'292.5';'315';'337.5'};
cellnames = {'0';'45';'90';'135';'180';'225';'270';'315'};
xticks([1:2:16])
set(gca,'xticklabel',cellnames)
xlabel('Direction')
ylabel('Significantly Responsive Cells')
legend({'Day 1', 'Day 2', 'Day 3'},'Location','northeast')
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\sig_resp_16dir'],'-dpdf') 

good2th2 = length(intersect(goodfit_D2,cells_all2'));
good2th3 = length(intersect(goodfit_D3,cells_all3'));
figure;
x = [1 2];
vals = [length(cells_all2) gdft_D2 good2th2; length(cells_all3) gdft_D3 good2th3]; 
cellnames = {'Day 2'; 'Day 3'};
b1 = bar(x,vals);
b1(1).FaceColor = [0.2334 0.4678 0.700];
b1(2).FaceColor = [0.6558 0.8238 0.8984];
b1(3).FaceColor = [0.9375 0.9708 1.0000];
set(gca,'xticklabel',cellnames,'FontSize',16)
% title('Significantly Responsive Cells vs Goodfit Cells')
ylabel('Number of cells')
legend({'Thresh Cells', 'Reliably Fit Cells', 'Overlap'},'Location','northwest')
ylim([0 50])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\goodfit_thresh_overlap'],'-dpdf') 

% scatter of pref ori diff
figure;
subplot(1,2,1)
scatter(prefOri_diff2,fitR2,'filled','MarkerFaceColor',[0 1 0])
hold on
scatter(prefOri_diff2(cells_all2),fitR2(cells_all2),'filled','MarkerFaceColor',[1 0 1])
refline(1,0)
axis square
xlabel('D1 vs D2 Pref Ori')
ylabel('D2 90th Per Pref Ori Fit')
set(gca,'FontSize',10)
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\prefOri_dif_D2scatter'],'-dpdf') 
subplot(1,2,2)
scatter(prefOri_diff3,fitR3,'filled','MarkerFaceColor',[0 1 0])
hold onsc
scatter(prefOri_diff3(cells_all3),fitR3(cells_all3),'filled','MarkerFaceColor',[1 0 1])
refline(1,0)
axis square
xlabel('D1 vs D3 Pref Ori')
ylabel('D3 90th Per Pref Ori Fit')
set(gca,'FontSize',10)
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' ref_date '_' mouse '\' ref_date '_' mouse '_' ref_str '\AcrossAllDays\prefOri_dif_D3scatter'],'-dpdf') 

%% dfof traces
% trace max dfof of each cell over 4 days
figure;
x = [1 2 3];
xticks([1 2 3]);
maxResp_all = [maxResp_D1' maxResp_D2' maxResp_D3'];
plot(x,maxResp_all)
cellnames = {'Day1'; 'Day2'; 'Day3'};
set(gca,'xticklabel',cellnames,'FontSize',16)
ylabel('max dF/F')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_maxDfofAcrossDays.pdf']),'-dpdf', '-bestfit')

% trace avg max dfof over 4 days
avgdfof1 = mean(maxResp_D1);
avgdfof2 = mean(maxResp_D2);
avgdfof3 = mean(maxResp_D3);
figure;
x = [1 2 3];
avgdfof = [avgdfof1 avgdfof2 avgdfof3];
cellnames = {'Day1'; 'Day2'; 'Day3'};
bar(x,avgdfof)
set(gca,'xticklabel',cellnames,'Fontsize',16)
ylabel('mean max dF/F')
title(mouse)
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_AvgMaxDfofAcrossDays.pdf']),'-dpdf', '-bestfit')

%% SNR
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

off_only4 = [ ];
for iTrial = 1:nTrials4
    off_only4 = [off_only4; squeeze(trial_tc4(nOff4/2:nOff4,iTrial,:))];
end
meanf4 = mean(off_only4,1);
stdf4 = std(off_only4,[],1);
SNR4 = meanf4./stdf4;
save(fullfile(fnout, [day4 '_' mouse], [day4 '_' mouse '_' run_str], [day4 '_' mouse '_' run_str '_off_only.mat']), 'off_only4')
skew4 = skewness(off_only4);

figure;
subplot(1,3,1);
scatter(SNR1,SNR2);
axis square
xlim([0 max(SNR1)+1])
ylim([0 max(SNR2)+1])
refline(1,0)
xlabel('SNR Day 1')
ylabel('SNR Day 2')
R = corrcoef(SNR1,SNR2);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_SNR1-2.pdf']),'-dpdf', '-bestfit')

avgSNR1 = mean(SNR1);
avgSNR2 = mean(SNR2);
avgSNR3 = mean(SNR3);
avgSNR4 = mean(SNR4);
stdSNR1 = std(SNR1);
stdSNR2 = std(SNR2);
stdSNR3 = std(SNR3);
stdSNR4 = std(SNR4);
figure;
x = [1 2 3 4];
y = [avgSNR1 avgSNR2 avgSNR3 avgSNR4];
err = [stdSNR1/sqrt(nCells1) stdSNR2/sqrt(nCells1) stdSNR3/sqrt(nCells1) stdSNR4/sqrt(nCells1)];
bar(x,y);
hold on
er = errorbar(x,y,err);
er.Color = [0 0 0];
er.LineStyle = 'none';
cellnames = {'Day1'; 'Day2'; 'Day3'; 'Day4'};
set(gca, 'xticklabel', cellnames)
ylabel('average SNR')
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_avgSNR.pdf']),'-dpdf', '-bestfit')

days = {ref_date; day2; day3; day4};
sz = size(off_only1);
% off_only = nan(sz(1),sz(2),length(days));
% off_only_mean = nan(sz(1),sz(2),length(days));
off_only_df = nan(sz(1),sz(2),length(days));
run_strs = strvcat(run_str, run_str2, run_str3, run_str);
for iday = 1:length(days)
    iDay = cell2mat(days(iday));
    run_str = run_strs(iday,:);
    off_only = load(fullfile(fnout, [iDay '_' mouse], [iDay '_' mouse '_' run_str], [iDay '_' mouse '_' run_str '_off_only.mat']));
    off_only = struct2array(off_only);
    off_only_mean = mean(off_only,1);
    off_only_prctile = prctile(off_only,10,1);
    sz = size(off_only,1);
    off_only_df(1:sz,:,iday) = (off_only - off_only_prctile)./off_only_mean;
end

[n, n2] = subplotn(length(goodfit3));
start = 1;
figure;
for iC = 1:length(goodfit3)
    iCell = goodfit3(iC);
    subplot(n,n2,start)
    tcOffsetPlot(squeeze(off_only_df(:,iCell,:)));
    xlim([2500 3500])
    til_str = ([num2str(chop(skew1(iCell),2)) ', ' num2str(chop(skew2(iCell),2)) ', ' num2str(chop(skew3(iCell),2)) ', ' num2str(chop(skew4(iCell),2))]);
    title(til_str)
    start = start + 1;
end
% [hleg, hobj, hout, mout] = legend({'Day 1', 'Day 2', 'Day 3', 'Day 4'},'Position',[0.5 0.9 0.2 0.1],'Orientation','horizontal');
% set(hobj,'linewidth',1.5);
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_TCoffset-mean.pdf']),'-dpdf', '-bestfit')

off_only_df1 = off_only1 - prctile(off_only1,10,1);
snr1 = mean(off_only_df1,1)./std(off_only_df1,[],1);
off_only_df2 = off_only2 - prctile(off_only2,10,1);
snr2 = mean(off_only_df2,1)./std(off_only_df2,[],1);
off_only_df3 = off_only3 - prctile(off_only3,10,1);
snr3 = mean(off_only_df3,1)./std(off_only_df3,[],1);
off_only_df4 = off_only4 - prctile(off_only4,10,1);
snr4 = mean(off_only_df4,1)./std(off_only_df4,[],1);

[n, n2] = subplotn(length(goodfit2));
start = 1;
figure;
for iC = 1:length(goodfit2)
    iCell = goodfit2(iC);
    subplot(n,n2,start)
    x = [1 2 3];
    y = [snr1(:,iCell) snr2(:,iCell) snr3(:,iCell)];
%     err = [stdSNR1/sqrt(nCells1) stdSNR2/sqrt(nCells1) stdSNR3/sqrt(nCells1) stdSNR4/sqrt(nCells1)];
    bar(x,y);
    cellnames = {'Day1'; 'Day2'; 'Day3'};
    set(gca, 'xticklabel', cellnames)
    ylabel('SNR')
    ylim([0 max(y)+.1])
    title(iCell)
    start = start + 1;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_goodfitSNR-allDays.pdf']),'-dpdf', '-bestfit')

%% Skewness
figure;
subplot(1,2,1);
scatter(skew1,skew2);
axis square
xlim([0 max(skew1)+.25])
ylim([0 max(skew2)+.25])
refline(1,0)
xlabel('Skew Day 1')
ylabel('Skew Day 2')
R = corrcoef(skew1,skew2);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

subplot(1,2,2)
scatter(skew1,skew3);
axis square
xlim([0 max(skew1)+.25])
ylim([0 max(skew3)+.25])
refline(1,0)
xlabel('Skew Day 1')
ylabel('Skew Day 3')
R = corrcoef(skew1,skew3);
disp(R(1,2));
str = ['    r = ',num2str(R(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% subplot(1,3,3)
% scatter(skew1,skew4);
% axis square
% xlim([0 max(skew1)+.25])
% ylim([0 max(skew4)+.25])
% refline(1,0)
% xlabel('Skew Day 1')
% ylabel('Skew Day 4')
% R = corrcoef(skew1,skew4);
% disp(R(1,2));
% str = ['    r = ',num2str(R(1,2))];
% T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
% set(T, 'fontsize', 10, 'verticalalignment', 'top', 'horizontalalignment', 'left');
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_skewScatter.pdf']),'-dpdf', '-bestfit')

[n, n2] = subplotn(length(goodfit2));
start = 1;
figure;
for iC = 1:length(goodfit2)
    iCell = goodfit2(iC);
    subplot(n,n2,start)
    x = [1 2 3];
    y = [skew1(:,iCell) skew2(:,iCell) skew3(:,iCell)];
%     err = [stdSNR1/sqrt(nCells1) stdSNR2/sqrt(nCells1) stdSNR3/sqrt(nCells1) stdSNR4/sqrt(nCells1)];
    bar(x,y);
    cellnames = {'Day1'; 'Day2'; 'Day3'};
    set(gca, 'xticklabel', cellnames)
    ylabel('Skew')
    ylim([0 max(y)+.1])
    title(iCell)
    start = start + 1;
end
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], ['SNR'], [ref_date '_' mouse '_' ref_str '_goodfitSkew.pdf']),'-dpdf', '-bestfit')

nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D2.avgResponseEaOri(iCell,:), oriTuning_D2.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
%         subplot(5, 4, start+2)
%         plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
%         hold on
%         errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
%         title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D3(iCell)) ' deg'])
%         subplot(5, 4, start+3)
%         plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
%         hold on
%         errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
%         title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D4(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 2 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-2.pdf']),'-dpdf', '-bestfit')

cell_list = intersect(goodfit3to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 3 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-3.pdf']),'-dpdf', '-bestfit')

cell_list = intersect(goodfit4to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(6, 2, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(6, 2, start+1)
        plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        start = start+2;
    end
end
suptitle([mouse ' day 1 to day 4 tuning'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTunin1-4.pdf']),'-dpdf', '-bestfit')

%% how each day individually aligns to day 1
cell_list = intersect(goodfit_D1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,1)
imagesc(data_dfof_max_ref);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 1')

cell_list = intersect(goodfit, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,2)
imagesc(reg2ref_dfof2);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 2')

cell_list = intersect(goodfit3to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,3)
imagesc(reg2ref_dfof3);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 3')

cell_list = intersect(goodfit4to1, unique(mask_cell2(ROI(1,:),ROI(2,:))));
subplot(1,4,4)
imagesc(reg2ref_dfof4);
for iC = 1:length(cell_list)
    iCell = cell_list(iC);
    if ~isnan(iCell)
        text(cell_stats(iCell).Centroid(1), cell_stats(iCell).Centroid(2), num2str(iCell), 'Color', 'white',...
            'Fontsize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle') 
        hold on
    end
end
xlim([ROI(2,1) ROI(2,end)])
ylim([ROI(1,1) ROI(1,end)])
axis square
axis off
title('Day 4')

suptitle([mouse ' days 2-4 goodfit with day 1'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_goodfitToDay1.pdf']),'-dpdf', '-bestfit')


%% tuning across all days
cell_list = intersect(goodfit3, unique(mask_cell2(ROI(1,:),ROI(2,:))));
nC = sum(~isnan(cell_list));
figure; 
start = 1;
for iC = 1:length(cell_list)
    if ~isnan(cell_list(iC))
        iCell = cell_list(iC);
        subplot(5, 4, start)
        plot(0:180, oriTuning_D1.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D1.avgResponseEaOri(iCell,:), oriTuning_D1.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D1(iCell)) ' deg'])
        subplot(5, 4, start+1)
        plot(0:180, oriTuning_D2.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D2.avgResponseEaOri(iCell,:), oriTuning_D2.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D2(iCell)) ' deg'])
        subplot(5, 4, start+2)
        plot(0:180, oriTuning_D3.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D3.avgResponseEaOri(iCell,:), oriTuning_D3.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D3(iCell)) ' deg'])
        subplot(5, 4, start+3)
        plot(0:180, oriTuning_D4.vonMisesFitAllCellsAllBoots(:,1,iCell))
        hold on
        errorbar(0:22.5:180-22.5, oriTuning_D4.avgResponseEaOri(iCell,:), oriTuning_D4.semResponseEaOri(iCell,:))
        title(['Cell ' num2str(iCell) '- ' num2str(prefOri_D4(iCell)) ' deg'])
        start = start+4;
    end
end
suptitle([mouse ' tuning across all days'])
print(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], ['AcrossAllDays'], [ref_date '_' mouse '_' ref_str '_alignTuning.pdf']),'-dpdf', '-bestfit')


