%% get path names
close all;clear all;clc;

ds = 'CrossOriRandPhase_15Hz_ExptList';
eval(ds)
nexp = length(expt);
for iexp = 6
rc = behavConstsAV;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).dirFolder;
coFolder = expt(iexp).coFolder;
time = expt(iexp).dirTime;
nrun = length(ImgFolder);
frameRateHz = params.frameRate;

run_str = catRunName(cell2mat(ImgFolder), nrun);
co_run_str = catRunName(cell2mat(coFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging Dir analysis \nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = ntrials*(nOn+nOff);
    
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,nframes);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    fprintf('Complete')
end
input = concatenateDataBlocks(temp);
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp
clear temp

toc

% register to cross-ori experiment

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_reg_shifts.mat']))
[out, data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data
% test stability
figure; 
subplot(2,2,1);
imagesc(data_reg_avg);
title('Direction run avg')
subplot(2,2,2);
imagesc(data_avg)
title('Cross-ori run avg')
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
rgb(:,:,2) = data_avg./max(data_avg(:));
subplot(2,2,3);
image(rgb)
title('Overlay')

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

% use cross-ori mask to get TCs

fprintf(['Loading masks from cross-ori runs: ' cell2mat(coFolder) '\n'])

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

%neuropil subtraction
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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

fprintf('\nNeuropil subtraction complete\n')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%% Direction analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(input.tGratingDirectionDeg);
nCells = size(npSub_tc,2);

data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff./2:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);

base_win = [nOff-5:nOff];
resp_win = [nOff+4:nOff+nOn];
tt_dir = (1-nOff:nOn).*(1000/frameRateHz);

data_dfof_dir = zeros(nOn+nOff, nCells, nDir);
h_dir = zeros(nCells,nDir);
p_dir = zeros(nCells,nDir);
dir_resp_avg = zeros(nCells,nDir,2);
dir_resp_mat = [];
dir_list = [];
for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
    data_dfof_dir(:,:,iDir) = mean(data_dfof(:,:,ind),3);
    dir_resp_avg(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
    dir_resp_avg(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,:,ind),1),[],3))./sqrt(length(ind));
    dir_resp_mat = [dir_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
    dir_list = [dir_list; iDir.*ones(length(ind),1)];
    [h_dir(:,iDir), p_dir(:,iDir)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
end
h_dir_all = zeros(1,nCells);
p_diranova = zeros(1,nCells); 
b_dir = nan(1,nCells);
k1_dir = nan(1,nCells);
R1_dir = nan(1,nCells);
R2_dir = nan(1,nCells);
u1_dir = nan(1,nCells);
sse_dir = nan(1,nCells);
R_square_dir = nan(1,nCells);
range = 1:1:360;
y_dir_fit = nan(nCells,length(range));
figure;
start = 1;
n = 0;
for iCell = 1:nCells
    if find(p_dir(iCell,:)<0.05/16)
        h_dir_all(iCell) = 1;
    elseif length(find(p_dir(iCell,:)<0.05/8))>=2
        h_dir_all(iCell) = 1;
    elseif length(find(p_dir(iCell,:)<0.05/4))>=4
        h_dir_all(iCell) = 1;
    end
    p_diranova(iCell) = anova1(dir_resp_mat(:,iCell),dir_list,'off');
    if start>25
        suptitle([mouse ' ' date ' Dir'])
     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuningAndFitsCells' num2str(n) '.pdf']),'-dpdf', '-bestfit')
        figure;
        n = n+1;
        start = 1;
    end
    subplot(5,5,start)
    errorbar(dirs, dir_resp_avg(iCell,:,1), dir_resp_avg(iCell,:,2), 'ok')
    hold on
    data = [dir_resp_avg(iCell,:,1) dir_resp_avg(iCell,1,1)];
    theta = [deg2rad(dirs) 2.*pi];
    [b_dir(:,iCell),k1_dir(:,iCell),R1_dir(:,iCell),R2_dir(:,iCell),u1_dir(:,iCell),u2_dir(:,iCell),sse_dir(:,iCell),R_square_dir(:,iCell)] ...
        = miaovonmisesfit_dir(theta,data);
    y_dir_fit(iCell,:) = b_dir(:,iCell)+R1_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell))-1))+R2_dir(:,iCell).*exp(k1_dir(:,iCell).*(cos(deg2rad(range)-u1_dir(:,iCell)-pi)-1));
    plot(range,y_dir_fit(iCell,:))
    title(num2str(chop(p_diranova(iCell),2)))
    start = start+1;
end
suptitle([mouse ' ' date ' Dir'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuningAndFitsCells' num2str(n) '.pdf']),'-dpdf', '-bestfit')
dirresp_ind = find(h_dir_all);
dirtuned_ind = find(p_diranova<0.05);
[max_val, max_ind_dir] = max(dir_resp_avg(:,:,1),[],2);
null_ind_dir = max_ind_dir-(nDir/2);
null_ind_dir(find(null_ind_dir<1)) = null_ind_dir(find(null_ind_dir<1))+nDir;
null_val = zeros(nCells,1);
for iC = 1:nCells
    null_val(iC,:) = dir_resp_avg(iC,null_ind_dir(iC),1);
end
null_val(find(null_val<0)) = 0;
DSI_resp = (max_val-null_val)./(max_val+null_val);
[max_val, max_ind_dir] = max(y_dir_fit,[],2);
null_ind_ori = max_ind_dir-(length(range)./2);
null_ind_ori(find(null_ind_ori<1)) = null_ind_ori(find(null_ind_ori<1))+length(range);
for iC = 1:nCells
    null_val(iC,:) = y_dir_fit(iC,null_ind_ori(iC));
end
null_val(find(null_val<0)) = 0;
DSI_fit = (max_val-null_val)./(max_val+null_val);

%ori tuning
ori_mat = dir_mat;
ori_mat(find(dir_mat>179)) = ori_mat(find(dir_mat>179))-180;
oris = unique(ori_mat);
nOri = length(oris);

data_dfof_ori = zeros(nOn+nOff, nCells, nOri);
h_ori = zeros(nCells,nOri);
p_ori = zeros(nCells,nOri);
ori_resp_avg = zeros(nCells,nOri,2);
ori_resp_mat = [];
ori_list = [];
for iOri = 1:nOri
    ind = find(ori_mat == dirs(iOri));
    data_dfof_ori(:,:,iOri) = mean(data_dfof(:,:,ind),3);
    ori_resp_avg(:,iOri,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
    ori_resp_avg(:,iOri,2) = squeeze(std(mean(data_dfof(resp_win,:,ind),1),[],3))./sqrt(length(ind));
    ori_resp_mat = [ori_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
    ori_list = [ori_list; iOri.*ones(length(ind),1)];
    [h_ori(:,iOri), p_ori(:,iOri)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
end
h_ori_all = zeros(1,nCells);
p_orianova = zeros(1,nCells);
b_ori = nan(1,nCells);
k1_ori = nan(1,nCells);
R1_ori = nan(1,nCells);
u1_ori = nan(1,nCells);
sse_ori = nan(1,nCells);
R_square_ori = nan(1,nCells);
range = 1:1:180;
y_ori_fit = nan(nCells,length(range));
figure;
start = 1;
n = 0;
for iCell = 1:nCells
    if find(p_ori(iCell,:)<0.05/8)
        h_ori_all(iCell) = 1;
    elseif length(find(p_ori(iCell,:)<0.05/4))>=2
        h_ori_all(iCell) = 1;
    elseif length(find(p_ori(iCell,:)<0.05/2))>=4
        h_ori_all(iCell) = 1;
    end
    p_orianova(iCell) = anova1(ori_resp_mat(:,iCell),ori_list,'off');
    if start>25
        suptitle([mouse ' ' date ' Ori'])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFitsCells' num2str(start+((n-1).*25)) '-' num2str(25+((n-1).*25)) '.pdf']))
        figure;
        n = n+1;
        start = 1;
    end
    subplot(5,5,start)
    errorbar(oris, ori_resp_avg(iCell,:,1), ori_resp_avg(iCell,:,2), 'ok')
    hold on
    data = [ori_resp_avg(iCell,:,1) ori_resp_avg(iCell,1,1)];
    theta = [deg2rad(oris) pi];
    [b_ori(:,iCell),k1_ori(:,iCell),R1_ori(:,iCell),u1_ori(:,iCell),sse_ori(:,iCell),R_square_ori(:,iCell)] ...
        = miaovonmisesfit_ori(theta,data);
    y_ori_fit(iCell,:) = b_ori(:,iCell)+R1_ori(:,iCell).*exp(k1_ori(:,iCell).*(cos(2.*(deg2rad(range)-u1_ori(:,iCell)))-1));
    plot(range,y_ori_fit(iCell,:))
    title(num2str(chop(p_orianova(iCell),2)))
    start = start+1;
end
suptitle([mouse ' ' date ' Ori'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFitsCells' num2str(start+((n-1).*25)) '-' num2str(25+((n-1).*25)) '.pdf']))
oriresp_ind = find(h_ori_all);
orituned_ind = find(p_orianova<0.05);
[max_val, max_ind_ori] = max(ori_resp_avg(:,:,1),[],2);
opp_ind_ori = max_ind_ori-(nOri/2);
opp_ind_ori(find(opp_ind_ori<1)) = opp_ind_ori(find(opp_ind_ori<1))+nOri;
opp_val = zeros(nCells,1);
for iC = 1:nCells
    opp_val(iC,:) = ori_resp_avg(iC,opp_ind_ori(iC),1)';
end
opp_val(find(opp_val<0)) = 0;
OSI_resp = (max_val-opp_val)./(max_val+opp_val);
[max_val, max_ind_ori] = max(y_ori_fit,[],2);
opp_ind_ori = max_ind_ori-(length(range)./2);
opp_ind_ori(find(opp_ind_ori<1)) = opp_ind_ori(find(opp_ind_ori<1))+length(range);
for iC = 1:nCells
    opp_val(iC,:) = y_ori_fit(iC,opp_ind_ori(iC))';
end
opp_val(find(opp_val<0)) = 0;
OSI_fit = (max_val-opp_val)./(max_val+opp_val);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'max_ind_ori', 'max_ind_dir', 'data_dfof', 'resp_win', 'base_win', 'tt_dir', 'data_dfof_ori', 'data_dfof_dir', 'h_dir', 'dirresp_ind', 'dirtuned_ind', 'h_ori', 'oriresp_ind', 'orituned_ind', 'b_dir', 'k1_dir', 'R1_dir', 'R2_dir', 'u1_dir', 'R_square_dir','sse_dir', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori','sse_ori', 'R_square_ori','DSI_resp', 'DSI_fit', 'OSI_resp', 'OSI_fit')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'dir_mat', 'dirs', 'nDir', 'ori_mat', 'oris', 'nOri', 'nOn', 'nOff','frameRateHz')

end
