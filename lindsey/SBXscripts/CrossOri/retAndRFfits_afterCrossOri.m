%% retOnly script by kevin 
% modified from singleChannelTCScript.m
% then modified from newScript.m

% reads in Tuning Curves from single channel imaging data

%% get path names
clear all; close all; clc;

ds = 'TwoStimTogether_ExptList';
iexp = 14; 
rc = behavConstsAV;
eval(ds)

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = cell2mat(expt(iexp).retFolder);
time = cell2mat(expt(iexp).retTime);
nrun = size(ImgFolder,1);
frame_rate = params.frameRate;
        
run_str = catRunName(ImgFolder, nrun);
datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];

fprintf(['2P imaging retinotopy analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
data_fn = fullfile(lg_fn, 'Data\2P_images');
mworks_fn = fullfile(fn_base, 'Behavior\Data');
fnout = fullfile(lg_fn, 'Analysis\2P');

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    % load 2p imaging data
    CD = fullfile(data_fn, mouse, date, ImgFolder(irun,:));
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    %imgMatFile = ['001_000_' ImgFolder(irun,:) '.mat']; % for bad file name
    load(imgMatFile);
    
    % load behavior/experimental data
    fName = fullfile(mworks_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);
    
    % read in frames with sbxread
    nframes = info.config.frames;
    fprintf(['\nReading run ' num2str(irun) ' - ' num2str(nframes) ' frames \n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    %data_temp = sbxread(['001_000_' ImgFolder(irun,:)],0,nframes); % for bad file name
    toc
    
    fprintf('Reshaping data...\n')
    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Register data
ImgFolder = cell2mat(expt(iexp).coFolder);
nrun = size(ImgFolder,1);

run_str = catRunName(ImgFolder, nrun);
datemouserefrun = [date '_' mouse '_' run_str];

load(fullfile(fnout, datemouse, datemouserefrun, [datemouserefrun '_reg_shifts.mat']))

% register
fprintf('stackRegister_MA, using target from previous registration\n')
[outs, data_reg]=stackRegister(data,data_avg);
fprintf('Registration complete, now saving...\n')
mkdir(fullfile(fnout, datemouse, datemouserun))
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_input.mat']), 'input')
clear data % depending on memory

%% test stability
% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
regIntv = 5000;
nep = floor(size(data_reg,3)./regIntv);
[n n2] = subplotn(nep);
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

% figure 3 shows imagesq (scaled?) of the registered data 1-10000, then prints to FOV_avg.mat
fprintf('Examine FOV\n')
figure(3);
if nframes>10000
    imagesq(mean(data_reg(:,:,1:10000),3));
else
    imagesq(mean(data_reg,3));
end
truesize;
% print to file
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_FOV_avg.pdf']),'-dpdf')

%% Use cross-ori mask to get TCs

fprintf(['Loading masks from cross-ori runs: ' cell2mat(expt(iexp).coFolder) '\n'])

load(fullfile(fnout, [datemouse], [datemouserefrun], [datemouserefrun '_mask_cell.mat']))
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

save(fullfile(fnout, datemouse, [datemouserun '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

fprintf('\nNeuropil subtraction complete\n')

clear data_tc data_tc_down np_tc np_tc_down

%% find activated cells
% calculate dF/F
fprintf('\nBegin image analysis...\n')

% max by trial type
if isfield(input, 'nScansOn')
    % nScansOn -> passive vis stim ret
    fprintf('nScansOn method - get dF/F\n')
    
    % load defining variables
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
   
    
    % store vectors for azimuth and elevation
    % what is celleqel2mat_padded?
    Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
    El = celleqel2mat_padded(input.tGratingElevationDeg);
    Azs = unique(Az);
    Els = unique(El);
    if min(Els,[],2)<0
        Els = fliplr(Els);
    end
    
    nStim = length(Azs)*length(Els);
    fprintf([num2str(nStim) ' unique Az+El stimuli\n'])
    
    
    Stims = [];
    %zeros(nStim, 2);
    
    for iEl = 1:length(Els)
        % runs through all unique Els to find indices of all same El
        indE = find(El == Els(iEl));
        for iAz = 1:length(Azs)
            % runs through all unique Azs to find indices of all same Az
            % then choose common El and Az indices
            % average over selected indices for all pixels and frames
            indA = find(Az == Azs(iAz));
            ind = intersect(indE,indA);
            
            % stores combination to Stims and iterates start
            Stims = [Stims; Els(iEl) Azs(iAz)];
        end
    end
end

%% calculate tuning mat and plot
%get dF/F
nCells = size(npSub_tc,2);
tc_mat = zeros(nOn+nOff, nCells, ntrials);
fulltc_mat = zeros(nOn+nOff, ntrials); % 1 for overall image, 2 for cells average, 3 for neuropil average
% reshape into trials, could do this in one line with reshape?
for itrial = 1:ntrials
    tc_mat(:,:,itrial) = npSub_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)),:);
    fulltc_mat(:,itrial) = mean(npSub_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)),:),2);
end
tc_f = mean(tc_mat(nOff/2:nOff,:,:),1);
fulltc_f = mean(fulltc_mat(nOff/2:nOff,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
fulltc_dfof = (fulltc_mat - fulltc_f) ./ fulltc_f;
clear tc_mat tc_f fulltc_mat fulltc_f

dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);
fprintf('\nPlotting timecourses and measuring stimOn response\n')
tuning_mat_all = zeros(nStim, 2, nCells);
tuning_mat_dir = zeros(nStim, 2,nDir, nCells);
fulltuning_mat = zeros(nStim, 2, 4);
Ind_struct = [];
tt= (1-nOff:nOn)*(1000./frame_rate);
%figure;
for iCell = 1:nCells
    for iStim = 1:nStim
        indA = find(Az == Stims(iStim,2));
        indE = find(El == Stims(iStim,1));
        ind = intersect(indE,indA);
        tuning_mat_all(iStim,1,iCell) = mean(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),1),3);
        tuning_mat_all(iStim,2,iCell) = std(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
        for idir = 1:nDir
            ind_dir = intersect(ind,find(dir_mat==dirs(idir)));
            tuning_mat_dir(iStim,1,idir,iCell) = mean(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind_dir),1),3);
            tuning_mat_dir(iStim,2,idir,iCell) = std(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind_dir),1),[],3)./sqrt(length(ind_dir));
            Ind_struct(iStim).dir_trials{idir} = ind_dir;
        end
    end
%     ylim([-0.05 0.25])
%     vline(nOff)
%     start = start + 1;
end


fprintf('Plotting tuning maps\n')
figure;
start = 1;
f = 1;
if nCells<36
    [n n2] = subplotn(nCells);
else
    [n n2] = subplotn(36);
end
for iCell = 1:nCells
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    ret_mat = reshape(tuning_mat_all(:,1,iCell), [length(Azs) length(Els)]);
    ret_mat = ret_mat';
    imagesc(ret_mat)
    colormap gray
    %clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
    %clim([0 chop(max(tuning_mat(:,1,iCell),[],1),2)])
    title(num2str(chop(max(tuning_mat_all(:,1,iCell),[],1),2)))
    start = start +1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile(fnout, datemouse, datemouserun, [datemouserun '_Tuning.mat']), 'tc_dfof', 'tuning_mat_all','tuning_mat_dir', 'Stims', 'Ind_struct')

% plot tc and ret_mat for full and cell avg
figure;
ret_mat = reshape(mean(tuning_mat_all(:,1,:),3), [length(Azs) length(Els)]);
ret_mat = ret_mat';
imagesc(ret_mat)
colormap gray
title(['All Cell Average dF/F tuning map, max:' num2str(chop(max(fulltuning_mat(:,1,2),[],1),2))])
%% retinotopy for these cells

% ADD IN: normalize each cell by its own max

% examine heat map of cell responses over all trials
% this allows examination of long-term effects over trials (ex: eye goop)
% responses are average of all stim on time (nOff+1:nOn+nOff) for a trial
fprintf('Examine cell responses across all trials\n')
trialRespMap = squeeze(mean(tc_dfof(nOff+1:nOn+nOff,:,:),1));
maxR = max(trialRespMap(:));
trialRespMap = [trialRespMap;El/max(El)*maxR;Az/max(Az)*maxR]; %add rows of normalized El and Az
figure;
imagesc(trialRespMap);
title('Avg dF/F resp by cell,trial')
xlabel('Trial')
ylabel('Cell')
colorbar

fprintf('Re-order by stimulus\n')
trialRespMapSort = sortrows(trialRespMap',[nCells+1 nCells+2])';
figure;
imagesc(trialRespMapSort);
title('Avg dF/F resp by cell,stimulus')
xlabel('Stim (El,Az)')
ylabel('Cell')
colorbar

%% fit retinotopy data

close all

fprintf('\nBegin fitting retinotopy data...\n')

fprintf('Plot tc_dfof for all stims of cell 10\n')
figure;
for iCell = 10
    for iCond = 1:nStim
        subplot(7,7,iCond)
        ind_all = Ind_struct(iCond).all_trials;
        plot(squeeze(tc_dfof(:,iCell,ind_all)))
        ylim([-0.1 0.4])
    end
end

Fit_struct = [];
[AzAz, ElEl] = meshgrid(Azs,Els);
grid2.AzAz = AzAz;
grid2.ElEl = ElEl;

dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;
Nshuf = 500;
fprintf(['Nshuf = ' num2str(Nshuf) '\n'])
resp_dFoverF = squeeze(mean(tc_dfof(nOff:nOn+nOff,:,:),1));
base_dFoverF = squeeze(mean(tc_dfof(nOff/2:nOff,:,:),1));
p_ttest = zeros(nCells,nStim,nDir+1);
h_ttest = zeros(nCells,nStim,nDir+1);
h_all = zeros(1,nCells,nDir+1);
ret_run_str = run_str;

fprintf('Begin shuffling...\n')
for count_shuf = 0:Nshuf
    fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
    Im_mat_USE = zeros(nCells, nStim,nDir+1);
    for iCond = 1:nStim
        for idir = 1:nDir
            ind_all = Ind_struct(iCond).dir_trials{idir};
            if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
            else
                ind_all_1 = ind_all;
                [h_ttest(:,iCond,idir), p_ttest(:,iCond,idir)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
            end
            Im_mat_USE(:,iCond,idir) = mean(resp_dFoverF(:,ind_all_1),2);
        end
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;
            [h_ttest(:,iCond,nDir+1), p_ttest(:,iCond,nDir+1)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
        end
        Im_mat_USE(:,iCond,nDir+1) = mean(resp_dFoverF(:,ind_all_1),2);
    end
    
    ifig = 1;
    start = 1;
    for iCell = 1:nCells
        for idir = 1:nDir+1
            if count_shuf == 0
                if sum(h_ttest(iCell,:,idir),2) == 0
                    ind_p = find(p_ttest(iCell,:,idir)< 0.05./((nStim-1)/2));
                    if length(ind_p)<2
                        ind_p = find(p_ttest(iCell,:,idir)< 0.05./((nStim-1)/3));
                        if length(ind_p)<3
                            ind_p = find(p_ttest(iCell,:,idir)< 0.05./((nStim-1)/4));
                            if length(ind_p)<4
                                h_all(1,iCell,idir) = 0;
                            else
                                h_all(1,iCell,idir) = 1;
                            end
                        else
                            h_all(1,iCell,idir) = 1;
                        end
                    else
                        h_all(1,iCell,idir) = 1;
                    end
                else
                    h_all(1,iCell,idir) = 1;
                end
            end
            if count_shuf>0
                if h_all(1,iCell,idir) == 0
                    continue
                end
            end

            a = Im_mat_USE(iCell,:,idir);
            if max(a,[],2) > 0
                b = reshape(a',length(Azs),length(Els));
                data = b';
                if count_shuf == 0
                    if idir<3
                        PLOTIT_FIT = 1;
                    else
                        PLOTIT_FIT = 0;
                    end
                    SAVEALLDATA = 1;
                    Fit_2Dellipse_LG_Ret_KM % modified due to error from file saving in script, saves to kevin analysis folder
                    eval(['Fit_struct(iCell,idir).True.s_',' = s;']);
                else
                    SAVEALLDATA = 0;
                    PLOTIT_FIT = 0;
                    Fit_2Dellipse_LG_Ret_KM
                    eval(['Fit_struct(iCell,idir).Shuf(count_shuf).s_',' = s;']);
                end
            end
        end
    end
    if count_shuf == 0  
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_RFfits' num2str(ifig) '_dir' num2str(idir) '.pdf']);
        print(fn_out,'-dpdf')
    end
end
fprintf('Shuffling done, saving fit results\n')
props = whos('Fit_struct');
fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_Fit_struct.mat']);
if props.bytes>2000000000
    save(fn_out, 'Fit_struct','-v7.3')
else
    save(fn_out, 'Fit_struct')
end

Npars = 10;
fit_true_vec = zeros(nCells, 10, nDir+1);
fit_shuf_vec = NaN(nCells,10,Nshuf,nDir+1);
goodfit_ind = cell(1,nDir+1);
lbub_fits = NaN(nCells,Npars,5,nDir+1);

alpha_bound = .025;
ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
for idir = 1:nDir+1
resp_ind{idir} = find(h_all(:,:,idir)); % h_all indicates responsive cell (by t-test against baseline)

fprintf('Assessing goodness of fit\n')
if Nshuf>1
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell,idir).True)
            eval('tmp = Fit_struct(iCell,idir).True.s_.x;');
            eval('tmp = [tmp Fit_struct(iCell,idir).True.s_.Elhicut_50];');
            eval('tmp = [tmp Fit_struct(iCell,idir).True.s_.Azhicut_50];');
            eval('tmp = [tmp Fit_struct(iCell,idir).True.s_.Elhicut_10];');
            eval('tmp = [tmp Fit_struct(iCell,idir).True.s_.Azhicut_10];');
            % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
            fit_true_vec(iCell,:,idir) = tmp;
        end
    end
    
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell,idir).Shuf)
                eval('tmp = Fit_struct(iCell,idir).Shuf(count_shuf).s_.x;');
                eval('tmp = [tmp Fit_struct(iCell,idir).Shuf(count_shuf).s_.Elhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell,idir).Shuf(count_shuf).s_.Azhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell,idir).Shuf(count_shuf).s_.Elhicut_10];');
                eval('tmp = [tmp Fit_struct(iCell,idir).Shuf(count_shuf).s_.Azhicut_10];');
                % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
                fit_shuf_vec(iCell,:,count_shuf,idir) = tmp;
            end
        end
    end
    
    
    for iCell = 1:nCells
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1,idir) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2,idir) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3,idir) = mean(i); %mean
            lbub_fits(iCell,count2,5,idir) = std(i); %stdev
        end
        lbub_fits(iCell,:,4,idir) = fit_true_vec(iCell,:,idir); % true (no shuffle)
    end
end

lbub_diff = lbub_fits(:,:,2,:)-lbub_fits(:,:,1,:);

goodfit_ind{idir} = [];
for iCell = 1:nCells
    if lbub_diff(iCell,4,idir)<input.gratingAzimuthStepDeg*2
        if lbub_diff(iCell,5,idir)<input.gratingAzimuthStepDeg*2
            goodfit_ind{idir} = [goodfit_ind{idir} iCell];
        end
    end
end
fprintf(['#Good cells = ' num2str(length(goodfit_ind{idir})) ' (first pass)...\nNow removing RFs at ret perimeter\n'])

goodfit_ind2 = zeros(size(goodfit_ind{idir}));
for i=1:length(goodfit_ind{idir})
    if sum(round(lbub_fits(goodfit_ind{idir}(i),4,4))==[min(Azs) max(Azs)])
        continue
    elseif sum(round(lbub_fits(goodfit_ind{idir}(i),5,4))==[min(Els) max(Els)])
        continue
    end
    goodfit_ind2(i) = goodfit_ind{idir}(i);
end
goodfit_ind2(goodfit_ind2==0) = [];
goodfit_ind{idir} = goodfit_ind2;
fprintf(['#Good cells = ' num2str(length(goodfit_ind{idir})) ' (final)\nSaving good fits\n'])


figure;
subplot(2,2,1)
for i = goodfit_ind{idir}
    plot(lbub_fits(i,4,4,idir), lbub_fits(i,5,4,idir), 'o')
    hold on;
end
xlim([min(Azs,[],2) max(Azs,[],2)])
ylim([min(Els,[],2) max(Els,[],2)])
axis equal
title('RF center (elim)')
subplot(2,2,2)
for i = goodfit_ind{idir}
    ellipse(lbub_fits(i,2,4,idir), lbub_fits(i,3,4,idir), 0, lbub_fits(i,4,4,idir), lbub_fits(i,5,4,idir));
    hold on;
end
axis equal
title('RF- 1 sigma (elim)')
subplot(2,2,3)
erat = (lbub_fits(goodfit_ind{idir},3,4,idir).^2./lbub_fits(goodfit_ind{idir},2,4,idir).^2);
erat(erat>1) = 1./erat(erat>1);
e = sqrt(1-erat);
hist(e)
title('eccentricity')
subplot(2,2,4)
cdiff = lbub_fits(goodfit_ind{idir},3,4,idir).^2-lbub_fits(goodfit_ind{idir},2,4,idir).^2;
cdiff(cdiff<0) = -cdiff(cdiff<0);
c = sqrt(cdiff);
hist(c)
title('linear eccentricity (elim)')
set(gcf, 'Position', [0 0 800 1000]);
fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_RFs_Dir' num2str(idir) '.pdf']);
print(fn_out,'-dpdf')

figure;
subplot(2,2,1)
hist(lbub_fits(goodfit_ind{idir},2,4,idir))
title('Sigma azimuth (elim)')
subplot(2,2,2)
hist(lbub_fits(goodfit_ind{idir},3,4,idir))
title('Sigma elevation (elim)')
subplot(2,2,3)
a = lbub_fits(goodfit_ind{idir},3,4,idir).*lbub_fits(goodfit_ind{idir},2,4,idir).*pi;
hist(a)
title('Area (elim)')
subplot(2,2,4)
scatter(a, lbub_fits(goodfit_ind{idir},1,4,idir),'o')
xlabel('Area (elim)')
ylabel('Peak dF/F')
set(gcf, 'Position', [0 0 800 1000]);
fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_RFdists_Dir' num2str(idir) '.pdf']);
print(fn_out,'-dpdf')


% visualize retinotopic organization
% takes each of the goodfit_inds and colors masks by El+Az of RF center

retMap_El = NaN(size(mask_cell));
retMap_Az = retMap_El;
for i=1:length(goodfit_ind{idir})
    ind = find(mask_cell == goodfit_ind{idir}(i));
    retMap_El(ind) = lbub_fits(goodfit_ind{idir}(i),5,4,idir);
    retMap_Az(ind) = lbub_fits(goodfit_ind{idir}(i),4,4,idir);
end

imAlpha=ones(size(retMap_El));
imAlpha(isnan(retMap_El))=0; % set all unmasked pixels to alpha=0

figure;clf;
colormap default
subplot(1,2,1)
imagesc(retMap_El,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by El')
h = colorbar;
ylabel(h,'El (deg)','Rotation',270.0,'VerticalAlignment','bottom')
%set(get(h,'label'),'string','El (deg)','Rotation',270.0); 
set(gca,'color',0*[1 1 1]);
subplot(1,2,2)
imagesc(retMap_Az,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by Az')
h = colorbar;
ylabel(h,'Az (deg)','Rotation',270.0,'VerticalAlignment','bottom')
set(gca,'color',0*[1 1 1]); 
set(gcf, 'Position', [100,300,1200,400])
print(fullfile(fnout, datemouse, datemouserun, [datemouserun '_retMap_Dir' num2str(idir) '.pdf']), '-dpdf','-bestfit')

end

fn_out = fullfile(fnout, datemouse, datemouserun, [datemouserun '_lbub_fits.mat']);
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind')