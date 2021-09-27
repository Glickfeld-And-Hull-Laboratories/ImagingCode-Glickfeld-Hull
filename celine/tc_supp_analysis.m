% to perform additional analyses on already extracted timecourses
clear all; clear global; close all; clc
ds = 'con_ori_nonDART';
%ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds)


day_id =2;
%% identifying animal and run
mouse = expt(day_id).mouse;
date = expt(day_id).date;

imgFolder = expt(day_id).contrastxori_runs{1};

frame_rate = expt(day_id).frame_rate;

%setting my paths
fn_base = 'Z:\home\ACh\Analysis\2p_Analysis\';
fn = fullfile(fn_base,mouse,date,imgFolder);
cd(fn);

load('input.mat')
load('mask_cell.mat')
load('regOuts&Img.mat')
load('TCs.mat')
% get info
nCells = size(npSub_tc,2);
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingContrast,2);
tCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tCon);
nCon = length(cons);
ind_con = find(tCon == max(cons(:)));
tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;
oris = unique(tOri);
nOri = length(oris);

%% get dfof for each trial
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;
pad = double(stimStart-1);
resp_win = stimStart:double(stimEnd-(nOn/2)); %looking at the first half of the response period
base_win = 1:(nOff/2);

npSub_tc_new=npSub_tc(stimStart:size(npSub_tc,1),:); %trimming npSub_tc to correct size( sacrifiace one trial)
npSub_tc_new=padarray(npSub_tc_new,pad,0,'post');
data_tc_trial = reshape(npSub_tc_new, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(1:stimStart,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%% get responses to stimuli
data_resp = zeros(nCells, nOri, nCon,2);
h = zeros(nCells, nOri, nCon);
p = zeros(nCells, nOri, nCon);
tCon = tCon(:,1:ntrials);
tOri = tOri(:,1:ntrials);
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

resp=logical(h_all);


%% identify cells I want to analyze
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));
RedAll =find(mask_label);
nGreen = length(green_inds)
GreenResp = intersect(green_inds,find(resp));
nGreenResp = length(GreenResp)

%% narrow down to verticle stimuli but keep all contrasts
tc_trial_conAvg = nan(nOn+nOff,nCells,length(cons));
tc_trial_avrg=nan(nOn+nOff,nCells);
for i=1:nCells
    temp_TCs=data_dfof_trial(:,:,i);
    ori_inds = find(tDir==0); %these are the trials with verticle orientation
    %all contrasts
    tc_trial_avrg(:,i)=mean(temp_TCs(:,ori_inds),2);
    %to get a mean for each contrast, if there's more than one
    %I can later add this for size too
    for iCon=1:length(cons)
        con_inds=find(tCon==cons(iCon));
        temp_trials = intersect(ori_inds, con_inds);
        tc_trial_conAvg(:,i,iCon)=mean(temp_TCs(:,temp_trials),2);
    end
end

%% plotting

tc_green = tc_trial_avrg(:,GreenResp);
tc_red = tc_trial_avrg(:,RedAll);
tc_green_avrg=cell(1,2);
tc_red_avrg=cell(1,2);
tc_green_avrg{1} = mean(tc_green,2);%average tc for responsive green cells
tc_green_avrg{2}=std(tc_green,[],2);
tc_red_avrg{1} = mean(tc_red,2); %average tc for all red cells
tc_red_avrg{2}=std(tc_red,[],2);

%convert to se 
tc_green_avrg{2} = tc_green_avrg{2}/sqrt(size(GreenResp,1));
tc_red_avrg{2} = tc_red_avrg{2}/sqrt(size(RedAll,2));
figure
t=1:(size(tc_green,1));
t=(t-pad)/frame_rate;
shadedErrorBar(t,tc_red_avrg{1},tc_red_avrg{2},'lineprops','r');
%plot(x,tc_red,'r');
%ylim([-.02 .23]);
hold on
shadedErrorBar(t,tc_green_avrg{1},tc_green_avrg{2});
title('responsive green and all red');
saveas(gcf,fullfile(fn,[mouse '-' date 'tcPlot.jpg']));


%% raw timecourses for all stimuli / all cells

tc_cell_avrg_all = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg_all = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg_all = mean(tc_cell_avrg_all,2);%average over trials and cells

figure;
ylim([min(min(tc_trial_avrg_all)) max(max(tc_trial_avrg_all))]);
patch([stimStart stimEnd stimEnd stimStart], [max(ylim) max(ylim) min(ylim) min(ylim)],'b','FaceAlpha',.25)
hold on
plot(tc_trial_avrg_all(:,green_inds), 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_trial_avrg_all(:,RedAll), 'LineWidth',.005, 'color','r');
hold on
plot(tc_cell_trial_avrg_all, 'LineWidth',2, 'color','k');
title('Timecourses all cells all stimuli');
hold off
saveas(gcf,fullfile(fn,[mouse '-' date 'raw_tc.jpg']));
%% identify cells that were visually responsive and that show a late deflection
% will look at responses to preferred stimuli, pyr cells using tc_green
deflect_ratio = tc_green(32,:)./tc_green(35,:);

figure
subplot(1,3,1)
plot(tc_green(:,find(deflect_ratio<.5)));
hold on;
plot(mean(tc_green(:,find(deflect_ratio<0.5)),2), 'LineWidth',2, 'color','k');
title('ratio < 0.5');

subplot(1,3,2)
plot(tc_green(:,find(deflect_ratio>=0.5|deflect_ratio<1)));
hold on;
plot(mean(tc_green(:,find(deflect_ratio>=0.5|deflect_ratio<1)),2), 'LineWidth',2, 'color','k');
title('ratio bet 0.5 and 1');

subplot(1,3,3)
plot(tc_green(:,find(deflect_ratio>=1)));
hold on;
plot(mean(tc_green(:,find(deflect_ratio>=1)),2), 'LineWidth',2, 'color','k');
title('ratio > 1');
%%  plotting cells individualls - I want to split this up by contrast, but keep at preferred ori
figure;
for i = 1:nGreenResp
    subplot(4,6,i);
    plot(tc_green(:,i), 'LineWidth',.005,'color',[.25 .25 .25]);
    ylim([min(min(tc_green)) max(max(tc_green))]);
    hline(0)
    vline(31)
end