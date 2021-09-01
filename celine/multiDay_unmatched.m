clear all; clear global; close all
clc

day_id = 43; %enter the experiment ID here

ds = 'DART_V1_contrast_ori'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds)
base = 'Z:\home\ACh\Analysis\2p_analysis\';


%% load data for day

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;

fn = fullfile(rc.achAnalysis,'2p_analysis',mouse,expDate); %can make this flexible if folder structure is different


run =  cell2mat(expt(day_id).contrastxori_runs);
time =  cell2mat(expt(day_id).contrastxori_time);


if strcmp(expt(day_id).data_loc,'ACh')
    root = rc.achData;
    CD = fullfile(base,mouse, expDate,run);
    dat = 'data-i''';
else
    fprintf('Cannot find data'); %in case the location is not ACh
end

cd(CD);

if username == 'celine'
    fName = ['Z:\Behavior\Data\' dat mouse '''-' expDate '-' time '.mat'];    
else
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' times{irun} '.mat'];
end
load(fName); %load the mworks behavioral file


out_path = fullfile(base,mouse, expDate, 'multiDay_analysis');
mkdir(out_path);


%create a text file to save basic results
fid = fopen('results all.txt','wt'); 

load('TCs.mat');
load('input.mat');
load('mask_cell.mat');

%% find key values 
nOn = input.nScansOn;
nOff = input.nScansOff;
nTrials = size(input.tGratingContrast,2);
nCells = size(npSub_tc,2);
red_ind = find(mask_label);

nRed = length(red_ind);
PC_ind = find(red_ind ==0);

fprintf(fid, ['\nTotal number of cells = ' num2str(nCells)]);
fprintf(fid, ['\n Red cells = ' num2str(nRed)]);


tCon = celleqel2mat_padded(input.tGratingContrast);
tCon = tCon(1:720);
cons = unique(tCon);
nCon = length(cons);
ind_con = find(tCon == max(cons(:)));
tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
tDir = tDir(1:720);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;
oris = unique(tOri);
nOri = length(oris);


%% 1 timecourses
%now we have the raw tc, the np tc, and the np-subtracted tc
% first plot timecourses
data_tc_trial = reshape(npSub_tc, [nOn+nOff,nTrials,nCells]);
data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg = mean(tc_cell_avrg,2);%average over trials and cells
tc_trial_avrg_PC = tc_trial_avrg(:,PC_ind);
tc_trial_avrg_IN = tc_trial_avrg(:,red_ind);

figure;
plot(tc_trial_avrg_PC, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_trial_avrg_IN, 'LineWidth',.005,'color',[.7 0 0]);
hold on;
plot(tc_cell_trial_avrg, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Timecourses after np subtraction');
dim = [.2 .5 .3 .3];
str = strcat(num2str(nCells),' cells averaged over', num2str(nTrials), ' trials') ;
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off
print(fullfile(out_path,['timecourses']),'-dpdf');

% timecourses without np subtracted
data_tc_trial_noSub = reshape(data_tc, [nOn+nOff,nTrials,nCells]);
data_f_trial_noSub = mean(data_tc_trial_noSub(nOff/2:nOff,:,:),1);
data_dfof_trial_noSub = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial_noSub, data_f_trial_noSub), data_f_trial);

%looking at data with np subtracted
tc_cell_avrg_noSub = mean(data_dfof_trial_noSub,3);%average pver cells, one row per trial
tc_trial_avrg_noSub  = squeeze(mean(data_dfof_trial_noSub,2));%average over trials, one row per cell
tc_cell_trial_avrg_noSub  = mean(tc_trial_avrg_noSub,2);%average over trials and cells

figure;
plot(tc_trial_avrg_noSub, 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_cell_trial_avrg_noSub, 'LineWidth',2, 'color','k');
hold on;
vline(60,'g')
title('Timecourses before np subtraction');
dim = [.2 .5 .3 .3];
str = strcat(num2str(nCells),' cells averaged over', num2str(nTrials), ' trials') ;
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off


%% 2 plotting activated cells - this code come straight from multiDay_cellSelect

%split into baseline and response windows, run paired t-test to see if
%cells have a significant response (elevation in df/f comapred to baseline)
%for any orientations/contrasts
resp_win = nOff+2:nOn+nOff;
base_win = nOff/2:nOff;
data_resp = zeros(nCells, nOri, nCon,2);
h = zeros(nCells, nOri, nCon);
p = zeros(nCells, nOri, nCon);
for iOri = 1:nOri
    ind_ori = find(tOri == oris(iOri));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        ind = intersect(ind_ori,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
        data_resp(:,iOri,iCon,1) = squeeze(mean(mean(data_dfof_trial(resp_win,ind,:),1),2));
        data_resp(:,iOri,iCon,2) = squeeze(std(mean(data_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(mean(data_dfof_trial(resp_win,ind,:),1), mean(data_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri.*3-1));
    end
end


h_all = sum(sum(h,2),3); 
nResp = length(find(h_all))
fprintf(fid, ['\nTotal number of cells with significant response to at least\n one stimulus = ' num2str(nResp)]);


if length(find(h_all))<36
    [n n2] = subplotn(length(find(h_all)));
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

% plot ori tuning at each contrast for significantly responsive cells only
figure;
movegui('center')
start = 1;
pref_ori = zeros(1,nCells);
for iCell = 1:nCells
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    if find(find(h_all)==iCell)
        for iCon = 1:nCon
            errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-o')
            hold on
        end
        if find(find(mask_label)==iCell)
            title('R')
        end
        start= start+1;
        ylim([-0.1 inf])
    end
    [max_val, pref_ori(1,iCell)] = max(mean(data_resp(iCell,:,:,1),3),[],2);
end


% plots contrast preference at preferred orientation
figure;
movegui('center')
start = 1;
for iCell = 1:nCells
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)
    if find(find(h_all)==iCell)
        errorbar(cons, squeeze(data_resp(iCell,pref_ori(iCell),:,1)), squeeze(data_resp(iCell,pref_ori(iCell),:,2)),'-o')
        if find(find(mask_label)==iCell)
            title('R')
        end
        ylim([-0.1 inf])
        start = start+1;
    end
end

