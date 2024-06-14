%% get path names
close all;clear all;clc;
ds = 'DART_behavior_ExptList';
experimentFolder = 'SST_behavior';
iexp = 3; 
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

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [database mouse '\' date '\' ImgFolder{irun}];
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

data_avg = mean(data(:,:,20001:20500),3);
%% Register data
[database mouse '\' date '\' ImgFolder{irun}]
if exist(fullfile(base, mouse,date, string(ImgFolder),'\regOuts&Img.mat'))
    load(fullfile(base, mouse,date, string(ImgFolder),'regOuts&Img.mat'))
    save(fullfile(base, mouse,date, string(ImgFolder),'input.mat'), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, mouse,date, string(ImgFolder)))
    save(fullfile(base, mouse,date, string(ImgFolder),'regOuts&Img.mat'), 'out', 'data_avg')
    save(fullfile(base, mouse,date, string(ImgFolder),'input.mat'), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*step):500+((i-1)*step)),3)); title([num2str(1+((i-1)*step)) '-' num2str(500+((i-1)*step))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

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
        if cStart(itrial)+19 < sz(3)
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
b1 = celleqel2mat_padded(input.tBlock2TrialNumber);
targSpeed = celleqel2mat_padded(input.tDotSpeedDPS);
spds = cell(1,2);
spds{1} = unique(targSpeed(find(b1)));
spds{2} = unique(targSpeed(find(b1==0)));
nSpd = length(spds{2});
baseCon = celleqel2mat_padded(input.tBaseDotContrast);
cons = unique(baseCon);
nCon = length(cons);
data_dfof_spd = zeros(sz(1),sz(2),nSpd,nCon);
[n n2] = subplotn(nSpd);
for iCon = 1:nCon
    figure;
    ind_con = find(baseCon == cons(iCon));
    for ispd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(ispd)));
        data_dfof_spd(:,:,ispd,iCon) = nanmean(data_targ_dfof(:,:,ind_spd),3);
        subplot(n,n2,ispd)
        imagesc(data_dfof_spd(:,:,ispd,iCon))
        title(spds{iCon}(ispd))
    end
    suptitle([mouse ' ' date '- Base con = ' num2str(cons(iCon))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVbySpd_Con' num2str(cons(iCon)) '.pdf']),'-dpdf','-bestfit')
end

ind_con = find(baseCon == 1);
data_dfof_base = nanmean(data_base_dfof(:,:,ind_con),3);

data_dfof = cat(3, reshape(data_dfof_spd,[sz(1), sz(2), nCon*nSpd]), data_dfof_base);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)
data_dfof = cat(3, data_dfof, data_dfof_max);

down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
pixelCorr = getPixelCorrelationImage(data_reg_down);
figure; 
imagesc(pixelCorr)
data_dfof = cat(3, data_dfof, pixelCorr);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellSelect.mat']), 'data_dfof')
%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
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

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

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
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% target analysis
data_targ = nan(120,nCells,nTrials);
data_base = nan(120,nCells,nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        if cStart(itrial)+99 <= size(npSub_tc,1)
            data_base(:,:,itrial) = npSub_tc(cStart(itrial)-20:cStart(itrial)+99,:);
        end
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+99 <= size(npSub_tc,1)
            data_targ(:,:,itrial) = npSub_tc(cTarget(itrial)-20:cTarget(itrial)+99,:);
        end
    end
end
data_f_targ = nanmean(data_targ(1:20,:,:),1);
data_dfof_targ = bsxfun(@rdivide,bsxfun(@minus,data_targ,data_f_targ),data_f_targ);
data_f_base = nanmean(data_base(1:20,:,:),1);
data_dfof_base = bsxfun(@rdivide,bsxfun(@minus,data_base,data_f_base),data_f_base);

frameRateHz = double(input.frameRateHz);

base_win =19:21;
resp_win =25:29; 

figure;
if nCells<25
    ii = nCells;
else
    ii = 25;
end
for i = 1:ii
    subplot(5,5,i)
    plot(squeeze(nanmean(mean(data_dfof_targ(20:50,i,:),2),3)))
    vline(base_win-19)
    vline(resp_win-19)
end

tt = [-20:99].*(1000/frameRateHz);
figure;
subplot(2,1,1)
plot(tt,squeeze(nanmean(mean(data_dfof_base,2),3)));
vline((base_win-20).*(1000/frameRateHz))
vline((resp_win-20).*(1000/frameRateHz))
title('Baseline')
subplot(2,1,2)
plot(tt,squeeze(nanmean(mean(data_dfof_targ,2),3)));
vline((base_win-20).*(1000/frameRateHz))
vline((resp_win-20).*(1000/frameRateHz))
title('Target')


%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof_base', 'data_dfof_targ')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'baseCon', 'cons', 'nCon', 'targSpeed','spds','nSpd', 'nCells', 'frameRateHz', 'nTrials', 'cTarget', 'cStart', 'base_win','resp_win')

%% Speed analysis
late_win = resp_win+(2.*frameRateHz);
full_win = resp_win(1): resp_win(1)+(2.*frameRateHz);
figure;
spd_mat = zeros(nSpd,nCon,nCells,2);
spd_mat_late = zeros(nSpd,nCon,nCells,2);
h_early = zeros(nCon,nSpd,nCells);
p_early = zeros(nCon,nSpd,nCells);
h_late = zeros(nCon,nSpd,nCells);
p_late = zeros(nCon,nSpd,nCells);
h_supp = zeros(nCon,nSpd,nCells);
p_supp = zeros(nCon,nSpd,nCells);
data_dfof_spd = zeros(size(data_dfof_targ,1), nCells, nCon,nSpd);
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    Ax(iCon+(iCon-1)) = subplot(2,2,iCon+(iCon-1));
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,:,ind_spd),3),2))
        data_dfof_spd(:,:,iCon,iSpd) = nanmean(data_dfof_targ(:,:,ind_spd),3);
        hold on
        for iCell = 1:nCells
            [h_early(iCon,iSpd,iCell), p_early(iCon,iSpd,iCell)] = ttest(mean(permute(data_dfof_targ(resp_win,iCell,ind_spd),[1 3 2]),1),mean(permute(data_dfof_targ(base_win,iCell,ind_spd),[1 3 2]),1),'tail','right','alpha',0.05./(nSpd*nCon));
            [h_late(iCon,iSpd,iCell), p_late(iCon,iSpd,iCell)] = ttest(mean(permute(data_dfof_targ(late_win,iCell,ind_spd),[1 3 2]),1),mean(permute(data_dfof_targ(base_win,iCell,ind_spd),[1 3 2]),1),'tail','right','alpha',0.05./(nSpd*nCon));
            [h_supp(iCon,iSpd,iCell), p_supp(iCon,iSpd,iCell)] = ttest(mean(data_dfof_targ(late_win,iCell,ind_spd),1),mean(data_dfof_targ(base_win,iCell,ind_spd),1),'tail','left','alpha',0.05./(nSpd*nCon));

            spd_mat(iSpd,iCon,iCell,1) = mean(nanmean(data_dfof_targ(resp_win,iCell,ind_spd),3),1);
            spd_mat(iSpd,iCon,iCell,2) = nanstd(nanmean(data_dfof_targ(resp_win,iCell,ind_spd),1),[],3)./sqrt(length(ind_spd));
            spd_mat_late(iSpd,iCon,iCell,1) = mean(nanmean(data_dfof_targ(late_win,iCell,ind_spd),3),1);
            spd_mat_late(iSpd,iCon,iCell,2) = nanstd(nanmean(data_dfof_targ(late_win,iCell,ind_spd),1),[],3)./sqrt(length(ind_spd));
        end
    end
    y(iCon+(iCon-1),:) = Ax(iCon+(iCon-1)).YLim;
    title(['BaseCon = ' num2str(cons(iCon))]) 
    ylabel('dF/F')
    xlabel('Time from target (ms)')
    if iCon == 1
        legend(num2str(chop(spds{iCon},2)'))
    end
    Ax(iCon+(iCon-1)+1) = subplot(2,2,iCon+(iCon-1)+1);
    errorbar(spds{iCon}, mean(spd_mat(:,iCon,:,1),3), std(spd_mat(:,iCon,:,1),[],3)./sqrt(nCells))
    hold on
    errorbar(spds{iCon}, mean(spd_mat_late(:,iCon,:,1),3), std(spd_mat_late(:,iCon,:,1),[],3)./sqrt(nCells))
    y(iCon+(iCon-1)+1,:) = Ax(iCon+(iCon-1)+1).YLim;
    ylabel('dF/F')
    xlabel('Speed (DPS)')
    if iCon == 1
        legend({'early', 'late'})
    end
    title(['BaseCon = ' num2str(cons(iCon))]) 
end
min_y = min(y(:,1),[],1);
max_y = max(y(:,2),[],1);
set(Ax(1:4),'YLim',[min_y max_y])
suptitle([date ' ' mouse '- Avg speed tuning- n = ' num2str(nCells)])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgSpeedResp.pdf']),'-dpdf', '-bestfit')

good_ind = unique([find(sum(sum(h_early,1),2)); find(sum(sum(h_late,1),2))]);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellResp.mat']),'data_dfof_spd', 'good_ind', 'tt')

for iCell = 1:length(good_ind)
    figure;
    iC = good_ind(iCell);
    y = zeros(4,2);
    clear Ax
    for iCon = 1:nCon
        ind_con = find(baseCon == cons(iCon));
        Ax(iCon+(iCon-1)) = subplot(2,2,iCon+(iCon-1));
        for iSpd = 1:nSpd
            ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(iSpd)));
            plot(tt, nanmean(data_dfof_targ(:,iC,ind_spd),3))
            hold on
        end
        y(iCon+(iCon-1),:) = Ax(iCon+(iCon-1)).YLim;
        title(['BaseCon = ' num2str(cons(iCon))]) 
        ylabel('dF/F')
        xlabel('Time from target (ms)')
        if iCon == 1
            legend(num2str(chop(spds{iCon},2)'))
        end
        Ax(iCon+(iCon-1)+1) = subplot(2,2,iCon+(iCon-1)+1);
        errorbar(spds{iCon}, spd_mat(:,iCon,iC,1), spd_mat(:,iCon,iC,2))
        hold on
        errorbar(spds{iCon}, spd_mat_late(:,iCon,iC,1), spd_mat_late(:,iCon,iC,2))
        y(iCon+(iCon-1)+1,:) = Ax(iCon+(iCon-1)+1).YLim;
        ylabel('dF/F')
        xlabel('Speed (DPS)')
        if iCon == 1
            legend({'early', 'late'})
        end
        title(['BaseCon = ' num2str(cons(iCon))]) 
    end
    min_y = min(y(:,1),[],1);
    max_y = max(y(:,2),[],1);
    set(Ax(1:4),'YLim',[min_y max_y])
    suptitle([date ' ' mouse '- Cell #' num2str(iC)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgSpeedResp_Cell#' num2str(iC) '.pdf']),'-dpdf', '-bestfit')
end

supp_ind = find(h_supp(end,end,:));
high_resp_ind = find(h_late(end,end,:));
low_resp_ind = find(h_late(end,1,:));

figure;
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,low_resp_ind,ind_spd),3),2))
        hold on
    end
    title('Resp to low speed increment')
end
title(['n = ' num2str(size(low_resp_ind,1))])
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon+2)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,high_resp_ind,ind_spd),3),2))
        hold on
    end
    title('Resp to high speed increment')
end
title(['n = ' num2str(size(high_resp_ind,1))])
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon+4)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds{iCon}(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,supp_ind,ind_spd),3),2))
        hold on
    end
    title('Supp to high speed increment')
end
title(['n = ' num2str(size(supp_ind,1))])

%% distribution of speed prefereances
% load((fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat'])))
% load((fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat'])))
% load((fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat'])))
% load((fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellResp.mat'])))
close all
stimOn_fr = 23;
resp_win_sec = [200 700];
resp_win_fr = stimOn_fr + ceil(resp_win_sec./(1000/frameRateHz));
base_win_fr = resp_win_fr-resp_win_fr(1)+1;
h_ttest = zeros(nSpd+1,nCon,nCells);

anova_spd_test = zeros(1,nCells);
for iCell = 1:nCells
    anova_spd_mat = [];
    group_mat = [];
    for icon = 1:nCon
        ind_con = find(baseCon == cons(icon));
        for iSpd = 1:nSpd
            ind_spd = intersect(ind_con, find(targSpeed==spds{icon}(iSpd)));
            base_temp = squeeze(mean(data_dfof_targ(base_win_fr(1):base_win_fr(2),iCell,ind_spd),1));
            resp_temp = squeeze(mean(data_dfof_targ(resp_win_fr(1):resp_win_fr(2),iCell,ind_spd),1));
            [h_ttest(iSpd,icon,iCell)] = ttest(base_temp,resp_temp,'tail','both','alpha',0.05./(nSpd));
            if icon == 2
                anova_spd_mat = [anova_spd_mat; squeeze(mean(data_dfof_targ(resp_win_fr(1):resp_win_fr(2),iCell,ind_spd),1))];
                group_mat = [group_mat; iSpd.*ones(length(ind_spd),1)];
            end
        end
        base_temp = squeeze(mean(data_dfof_targ(base_win_fr(1):base_win_fr(2),iCell,:),1));
        resp_temp = squeeze(mean(data_dfof_targ(resp_win_fr(1):resp_win_fr(2),iCell,:),1));
        [h_ttest(nSpd+1,icon,iCell)] = ttest(base_temp,resp_temp,'tail','both','alpha',0.05./(nSpd));
        if icon == 2
            [anova_spd_test(1,iCell)] = anova1(anova_spd_mat,group_mat,'off');
        end
    end
end
anova_ind = find(anova_spd_test<0.05);
data_dfof_spd_inc_resp = squeeze(mean(data_dfof_spd(resp_win_fr(1):resp_win_fr(2),:,2,:),1));
[max_val max_ind] = max(data_dfof_spd_inc_resp(anova_ind,:),[],2);
spds{1} = [double(input.block2BaseDotSpeedDPS) spds{1}];
spds{2} = [double(input.baseDotSpeedDPS) spds{2}];
spds_log = cell(1,2);
for i=1:2
    spds_log{i} = log2(double(spds{i}));
end
max_ind(find(max_val<0)) = 0;

figure; 
subplot(1,2,1)
plot(tt,squeeze(mean(data_dfof_spd(:,anova_ind,2,:),2)))
ylabel('dF/F')
xlabel('Time (ms)')
legend(num2str(spds{2}(2:end)'))
title([num2str(length(anova_ind)) ' tuned cells'])
temp = ylim;
subplot(1,2,2)
plot(tt,squeeze(mean(data_dfof_spd(:,setdiff(good_ind,anova_ind),2,:),2)))
ylabel('dF/F')
xlabel('Time (ms)')
ylim(temp)
title([num2str(length(setdiff(good_ind,anova_ind))) ' responsive untuned cells'])
suptitle([date ' ' mouse])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_tunedVuntuned.pdf']),'-dpdf', '-bestfit')

figure; 
[n, edges, bin] = histcounts(max_ind);
[n n2] = subplotn(length(spds{2})+1);
subplot(n,n2,1)
hist(spds_log{2}(max_ind+1))
ylabel('Number of cells')
xlabel('log2(Speed)')


for i = 1:length(spds{2})
    subplot(n,n2,1+i)
    ind = find(max_ind==i-1);
    plot(tt,squeeze(mean(data_dfof_spd(:,anova_ind(ind),2,:),2)))
    ylabel('dF/F')
    xlabel('Time (ms)')
    title([num2str(spds{2}(i)) ' deg pref cells- n = ' num2str(length(ind))])
    ylim([-.15 .25])
end

suptitle([date ' ' mouse '- ' num2str(length(anova_ind)) ' tuned cells'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_tunedCellSummary.pdf']),'-dpdf', '-bestfit')

%%
spd_step = [spds{2}-input.baseDotSpeedDPS];
spd_resp_mat = nan(nCells,nCon,nSpd);
for icon = 1:nCon
    tr_resp_mat = cell(1,nSpd+1);
    ind_con = find(baseCon == cons(icon));
    targetTr = [];
    resp = [];
    for iSpd = 2:nSpd+1
        ind_spd = intersect(ind_con, find(targSpeed==spds{icon}(iSpd)));
        if length(ind_spd)>0
            resp = [resp squeeze(mean(data_dfof_targ(resp_win_fr(1):resp_win_fr(2),:,ind_spd),1))];
            targetTr = [targetTr spd_step(iSpd).*ones(1, length(ind_spd))];
            resp = [resp squeeze(mean(data_dfof_targ(base_win_fr(1):base_win_fr(2),:,ind_spd),1))];
            targetTr = [targetTr zeros(1, length(ind_spd))];
            tr_resp_mat{iSpd} = squeeze(mean(data_dfof_targ(resp_win_fr(1):resp_win_fr(2),:,ind_spd),1));
            tr_resp_mat{1} = [tr_resp_mat{1} squeeze(mean(data_dfof_targ(base_win_fr(1):base_win_fr(2),:,ind_spd),1))];
            spd_resp_mat(:,icon,iSpd-1) = nanmean(tr_resp_mat{iSpd},2);
        end
    end
end
targetTrInd = zeros(size(targetTr));
targetTrInd(find(targetTr)) = 1;

resp_ind = find(nansum(nansum(h_ttest,1),2));
inc_resp_ind = find(nansum(h_ttest(:,2,:),1));

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trResp.mat']),'spd_resp_mat', 'tr_resp_mat', 'spds', 'anova_ind', 'resp_ind', 'inc_resp_ind', 'resp_win_fr','base_win_fr', 'resp', 'targetTr')

%% single cell decode
trInd_train_all = 1:length(targetTr);
pctCorrectTarget_ho_singlecell = zeros(1,nCells);
yhat_singlecell = zeros(length(targetTr),nCells);
cellW = zeros(1,nCells);
for i = 1:nCells
    [pctCorrectTarget_ho_singlecell(i), yhat_singlecell(:,i), targetGLM(i)] = getPctCorr_hoData_train(resp(i,:)',targetTrInd',trInd_train_all,trInd_train_all,0.5);
    cellW(1,i) = targetGLM(i).beta(2);
    fprintf([num2str(i) '/' num2str(nCells) '\n'])
end        

data_dfof_spd_inc_resp = squeeze(mean(data_dfof_spd(resp_win_fr(1):resp_win_fr(2),:,2,:),1));
data_dfof_spd_inc_base = squeeze(nanmean(mean(data_dfof_spd(base_win_fr(1):base_win_fr(2),:,2,:),1),4));
[max_val max_ind] = max(data_dfof_spd_inc_resp,[],2);
[min_val min_ind] = min(data_dfof_spd_inc_resp,[],2);
val_diff = (max_val-min_val)./max_val;

data_dfof_spd_resp = squeeze(mean(data_dfof_spd(resp_win_fr(1):resp_win_fr(2),:,1,:),1));
[max_val_nobase max_ind_nobase] = max(data_dfof_spd_resp,[],2);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_singleCellDecode.mat']),'pctCorrectTarget_ho_singlecell', 'max_ind', 'max_val', 'max_ind_nobase', 'max_val_nobase', 'val_diff', 'yhat_singlecell', 'spd_step', 'cellW')

%% cell group decode


resp_high_ind = find(max(data_dfof_spd_inc_resp,[],2)>0.05);
n = [4 8 16 32];
figure;
p = [];
ci =[];
for i = 1:length(n)
    n_tot = length(resp_ind);
    for ir = 1:ceil(n_tot./n(i))
        rand_ind = randperm(n_tot,n(i));
        use_ind = resp_ind(rand_ind);
        [pctCorrectTarget_ho, yhat_temp, targetGLM] = getPctCorr_hoData_train(resp(use_ind,:)',targetTrInd',trInd_train_all,trInd_train_all,0.5);
        
        for iS = 1:nSpd+1
            ind = find(targetTr == spd_step(iS));
            [p(iS,ir,i) ci(iS,:,ir,i)] = binofit(sum(yhat_temp(ind)>0.5),length(ind));
        end
    end
    %subplot(3,3,i)
    %errorbar(spd_step, p, p'-ci(:,1), ci(:,2)-p')
    p(find(p==0))=NaN;
    errorbar(spd_step, nanmean(p(:,:,i),2), nanstd(p(:,:,i),[],2))
    hold on
    title(num2str(n(i)))
end

figure;
[n n2] = subplotn(length(use_ind));
for iC = 1:length(use_ind)
    iCell = use_ind(iC);
    subplot(n,n2,iC)
    for i = 1:nSpd+1
        errorbar(i,nanmean(tr_resp_mat{i}(iCell,:),2), nanstd(tr_resp_mat{i}(iCell,:),[],2),'k')
        hold on
    end
    title(num2str(targetGLM.beta(iC+1)))
    %ylim([-.1 0.5])
    xlim([0 6])
    hline(0)
end