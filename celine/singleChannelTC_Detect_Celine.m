%% get path names
close all;clear all;clc;
ds = 'DART_behavior_ExptList';
iexp = 12; 
eval(ds);

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).behavior_runs;
time = expt(iexp).behavior_time;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);
data_folder = expt(iexp).data_loc;
tquit = expt(iexp).trial_quit;
redrun = expt(iexp).redChannelRun;

if strcmp(data_folder,'ACh')
    data_base = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\',data_folder,'Data','2P_data');
    out_base = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\',data_folder,'Analysis','2p_analysis');
else
    data_base = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\',data_folder,'Data','2P_images');
    out_base = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\',data_folder,'Analysis','2P');
end

fprintf(['2P imaging behavior analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
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
    CD = fullfile(data_base, mouse,date,ImgFolder{irun});
    %CD = [out_base '\Data\2P_images\' [date '_' mouse] '\' ImgFolder{irun}];
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

data_avg = mean(data(:,:,35001:35500),3);
%% Register data

if exist(fullfile(out_base, mouse, date, run_str))
    load(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(out_base, mouse, date, run_str))
    save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*step):500+((i-1)*step)),3)); title([num2str(1+((i-1)*step)) '-' num2str(500+((i-1)*step))]); end
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(out_base, mouse, date, run_str,[date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

CD = fullfile(data_base, mouse,date,redrun);
cd(CD);
imgMatFile = [redrun '_000_000.mat'];
load(imgMatFile);
nframes = info.config.frames;
data_rg = sbxread(imgMatFile(1,1:11),0,nframes);
data_green = squeeze(data_rg(1,:,:,:));
[reg_out data_green_reg] = stackRegister(data_green,data_avg);
figure; subplot(2,1,1); imagesc(mean(data_green_reg,3)); title('Green'); axis off
data_red = squeeze(data_rg(2,:,:,:));
[reg_out data_red_reg] = stackRegister_MA(data_red,[],[],reg_out);
red_avg = mean(data_red_reg,3);
subplot(2,1,2); imagesc(mean(data_red_reg,3)); title('Red'); axis off
print(fullfile(out_base, mouse, date, run_str,[date '_' mouse '_' run_str '_red&green_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

cTarget = celleqel2mat_padded(input.cTargetOn);
cStart = celleqel2mat_padded(input.cLeverDown);
nTrials = length(cTarget);
sz = size(data_reg);
data_f_targ = zeros(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
%data_targ_tc = nan(50,nTrials);
for itrial = 1:nTrials    
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+19 < sz(3)
            data_f_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)-20:cTarget(itrial)-1),3);
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+5:cTarget(itrial)+8),3);
            %data_targ_tc(:,itrial) = squeeze(mean(mean(data_reg(:,:,cTarget(itrial)-20:cTarget(itrial)+29),1),2));
        end
    end
end
data_targ_dfof = (data_targ-data_f_targ)./data_f_targ;
b1 = celleqel2mat_padded(input.tBlock2TrialNumber);
SIx = strcmp(input.trialOutcomeCell,'success');
MIx = strcmp(input.trialOutcomeCell,'ignore');
tContrast = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tContrast);
nCon = length(cons);
data_dfof_con = zeros(sz(1),sz(2),nCon);
[n n2] = subplotn(nCon);
figure;
for iCon = 1:nCon
    ind_con = find(tContrast == cons(iCon));
    data_dfof_con(:,:,iCon) = nanmean(data_targ_dfof(:,:,ind_con),3);
    subplot(n,n2,iCon)
    imagesc(data_dfof_con(:,:,iCon))
    title(cons(iCon))
end
suptitle([mouse ' ' date])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_FOVbySpd_Con' num2str(cons(iCon)) '.pdf']),'-dpdf','-bestfit')

data_dfof = data_dfof_con;
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

save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_cellSelect.mat']), 'data_dfof')
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
shadeimg = imShade(red_avg,mask_cell);
figure; imagesc(shadeimg)
print(fullfile(out_base, mouse, date, run_str,[date '_' mouse '_' run_str '_maskRedCells.pdf']),'-dpdf', '-bestfit')
%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(out_base,mouse, date, run_str, [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np', 'red_avg')

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

save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% target analysis
data_targ = nan(120,nCells,nTrials);
for itrial = 1:nTrials
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+99 <= size(npSub_tc,1)
            data_targ(:,:,itrial) = npSub_tc(cTarget(itrial)-20:cTarget(itrial)+99,:);
        end
    end
end
data_f_targ = nanmean(data_targ(1:20,:,:),1);
data_dfof_targ = bsxfun(@rdivide,bsxfun(@minus,data_targ,data_f_targ),data_f_targ);

frameRateHz = double(input.frameRateHz);

base_win =19:21;
resp_win =24:26; 

tt = [-20:99].*(1000/frameRateHz);
figure;
plot(tt,squeeze(nanmean(mean(data_dfof_targ,2),3)));
vline((base_win-20).*(1000/frameRateHz))
vline((resp_win-20).*(1000/frameRateHz))
title('Target')
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_avgTC.pdf']),'-dpdf','-bestfit')

%%

save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof_targ')
save(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_stimData.mat']),'b1','tContrast','cons', 'nCon', 'nCells', 'frameRateHz', 'nTrials', 'cTarget', 'cStart', 'base_win','resp_win')

%% Contrast analysis
nblock = length(unique(b1));
con_resp_mat = zeros(nCon,nCells,2);
con_block_resp_mat = zeros(nCon,nCells,nblock,2);
data_dfof_con = zeros(size(data_dfof_targ,1), nCells, nCon);
data_dfof_con_block = zeros(size(data_dfof_targ,1), nCells, nCon, nblock);
h = zeros(nCon,nCells);
p = zeros(nCon,nCells);
trialn = zeros(nCon, nblock);

for iCon = 1:nCon
    ind_con = find(tContrast == cons(iCon));
    data_dfof_con(:,:,iCon) = nanmean(data_dfof_targ(:,:,ind_con),3);
    con_resp_mat(iCon,:,1) = mean(nanmean(data_dfof_targ(resp_win,:,ind_con)-data_dfof_targ(base_win,:,ind_con),3),1);
    con_resp_mat(iCon,:,2) = nanstd(nanmean(data_dfof_targ(resp_win,:,ind_con)-data_dfof_targ(base_win,:,ind_con),1),[],3)./sqrt(length(ind_con));
    for iCell = 1:nCells
        [h(iCon,iCell), p(iCon,iCell)] = ttest(mean(permute(data_dfof_targ(resp_win,iCell,ind_con),[1 3 2]),1),mean(permute(data_dfof_targ(base_win,iCell,ind_con),[1 3 2]),1),'tail','right','alpha',0.05./(nCon));
    end
    for ib = 1:nblock
        ind_b = intersect(find(MIx+SIx),intersect(1:tquit,intersect(ind_con,find(b1==ib-1))));
        data_dfof_con_block(:,:,iCon,ib) = nanmean(data_dfof_targ(:,:,ind_b),3);
        con_block_resp_mat(iCon,:,ib,1) = mean(nanmean(data_dfof_targ(resp_win,:,ind_b)-data_dfof_targ(base_win,:,ind_b),3),1);
        con_block_resp_mat(iCon,:,ib,2) = nanstd(nanmean(data_dfof_targ(resp_win,:,ind_b)-data_dfof_targ(base_win,:,ind_b),1),[],3)./sqrt(length(ind_b));
        trialn(iCon,ib) = length(ind_b);
    end
end

good_ind = find(h(end,:));
figure;
[n n2] = subplotn(nCon);
for iCon = 1:nCon
    subplot(n,n2,iCon)
    plot(tt, mean(data_dfof_con(:,good_ind,iCon),2))
    ylabel('dF/F')
    xlabel('Time from target (ms)')
    title(num2str(cons(iCon)))
    ylim([-0.02 .1])
end
suptitle([date ' ' mouse '- Avg Contrast Response- All trials- n = ' num2str(length(good_ind))])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_contrastRespTC_allTrials.pdf']),'-dpdf', '-bestfit')


figure;
[n n2] = subplotn(nCon);
for iCon = 1:nCon
    for ib = 1:nblock
        subplot(n,n2,iCon)
        plot(tt, mean(data_dfof_con_block(:,good_ind,iCon,ib),2))
        hold on
        ylabel('dF/F')
        xlabel('Time from target (ms)')
        title(num2str(cons(iCon)))
        ylim([-0.02 .1])
    end
end
suptitle([date ' ' mouse '- Avg Contrast Response- By block- n = ' num2str(length(good_ind))])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_contrastRespTC_byBlock.pdf']),'-dpdf', '-bestfit')

figure;
[n n2] = subplotn(length(good_ind));
for i = 1:length(good_ind)
    subplot(n,n2,i)
    plot(tt, data_dfof_con(:,good_ind(i),iCon))
    title(num2str(good_ind(i)))
end
suptitle([date ' ' mouse '- Avg Response- high con- n = ' num2str(length(good_ind))])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_highConRespTC_allCells.pdf']),'-dpdf', '-bestfit')



figure;
[n n2] = subplotn(length(good_ind));
for i = 1:length(good_ind)
    subplot(n,n2,i)
    for ib = 1:nblock
        errorbar(cons, con_block_resp_mat(:,good_ind(i),ib,1), con_block_resp_mat(:,good_ind(i),ib,2))
        hold on
    end
    xlabel('Contrast')
    title(num2str(good_ind(i)))
end
suptitle([date ' ' mouse '- Avg Contrast Response- By block- n = ' num2str(length(good_ind))])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_contrastResp_byBlock.pdf']),'-dpdf', '-bestfit')

figure; 
for ib = 1:nblock
    errorbar(cons, mean(con_block_resp_mat(:,good_ind,ib,1),2), std(con_block_resp_mat(:,good_ind,ib,1),[],2)./sqrt(length(good_ind)))
    hold on
end
xlabel('Contrast')
title([date ' ' mouse '- Avg Contrast Response- By block- n = ' num2str(length(good_ind))])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_AvgContrastResp_byBlock.pdf']),'-dpdf', '-bestfit')


save(fullfile(out_base, mouse, date, run_str,  [date '_' mouse '_' run_str '_cellResp.mat']),'trialn', 'data_dfof_con', 'data_dfof_con_block', 'con_block_resp_mat', 'con_resp_mat', 'good_ind', 'tt')

good_ind_mask = zeros(size(mask_cell));
for i = 1:length(good_ind)
    good_ind_mask(find(mask_cell==good_ind(i))) = 1;
end
shadeimg = imShade(red_avg,good_ind_mask);
figure; imagesc(shadeimg)
%% Behavior
hr = zeros(nCon,nblock);
ci95 = zeros(2,nCon,nblock);
for iCon = 1:nCon
    ind_con = find(tContrast==cons(iCon));
    for ib = 1:nblock
        ind_b = find(b1==ib-1);
        ind_correct = intersect(1:tquit,intersect(find(SIx),intersect(ind_con,ind_b)));
        ind_miss = intersect(1:tquit,intersect(find(MIx),intersect(ind_con,ind_b)));
        ntrials = length([ind_correct ind_miss]);
        [hr(iCon,ib) ci(:,iCon,ib)] = binofit(length(ind_correct),ntrials);
    end
end
ci_flip = permute(ci,[2 3 1]);
figure; 
for ib = 1:nblock
    errorbar(cons, hr(:,ib),hr(:,ib)-ci_flip(:,ib,1),ci_flip(:,ib,2)-hr(:,ib))
    hold on
end
xlabel('Contrast')
ylabel('Hit rate')

title([date ' ' mouse '- hit rate by block'])
print(fullfile(out_base, mouse, date, run_str, [date '_' mouse '_' run_str '_HitRateByBlock.pdf']),'-dpdf', '-bestfit')

