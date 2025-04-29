clear all; clear global; close all;
clc
% ----- THIS SCRIPT ONLY USES ITI FRAMES -----

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
NoIntermediateFigures = 1;

fn_analysis_root = 'G:\home\ACh\Analysis\2p_analysis';
fn_data_root = 'G:\home\ACh\Data\2p_data';
fn_out_root = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask';
cd(fn_out_root)



load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask\TCs_SST.mat');
load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask\TCs_PV.mat');
load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\synchro\iti_sync_PV.mat');
load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\synchro\iti_sync_SST.mat');

%% TH
nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrial = 960;
nbins = 10;

ds = 'DART_expt_info';
eval(ds);

nMice = size(TCs_PV,1);
mice = TCs_PV(:,2);
std_cell_PV_ITI = cell(nMice,2);
std_cell_PV_ITI(:,2) = mice;
% skew_cell_PV_ITI = cell(nMice,2);
% skew_cell_PV_ITI(:,2) = mice;
frac_rms_PV_ITI = cell(nMice,2);
frac_rms_PV_ITI(:,2) = mice;
burst_PV = cell(nMice,2);
burst_PV(:,2) = mice;
ent_PV = cell(nMice,2);
ent_PV(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_PV{mouse,1};
    mouse_id = TCs_PV{mouse,2};
    nSesh = size(mouse_data,1);

    perc_mat = NaN(nSesh,1);
    std_mat = NaN(nSesh,1);
    burst_mat = NaN(nSesh,1);
    ent_mat = NaN(nSesh,1);
    % skew_mat = NaN(nSesh,1);

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sesh = 1:size(mouse_data,1)
        curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        nplots = 2; % hardcoded for 2 plots per mouse

        if sesh == 1
            this_tp_title = 'pre';
        else
            this_tp_title = 'post';
        end

        % only use ITI data
        data_trial = reshape(curr_data,nTot,[]);
        data_trial_ITI = data_trial(nOff/2+1:nOff,:);
        % ntrials = size(data_trial_ITI,2);
        % temp_ent_trial = NaN(ntrials,1);
        % for tr = 1:ntrials
        %     temp_ent_trial(tr) = perm_entropy(data_trial_ITI(:,tr),3,1);
        % end
        % ent_mat(sesh) = nanmean(temp_ent_trial);

        % get only ITI frames 
        curr_sesh = reshape(data_trial_ITI,[],1);

        ent_mat(sesh) = perm_entropy(curr_sesh,3,1);

        this_burst = abs(diff(data_trial_ITI));
        mean_burst = mean(mean(this_burst));
        burst_mat(sesh) = mean_burst;

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        % skew_curr_sesh = skewness(curr_sesh);
        % frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        sorted = sort(curr_sesh);
        bot20 = sorted(1:length(curr_sesh)*0.75);
        bot20_avg = mean(bot20);
        bot20_std = std(bot20);
        rows_vector = repelem(length(curr_sesh)/nbins,nbins);
        nelem = length(curr_sesh)/nbins;
        split_sesh = mat2cell(curr_sesh,rows_vector);
        

        % perc_mat(sesh) = frac_past;
        std_mat(sesh) = std_curr_sesh;
        % skew_mat(sesh) = skew_curr_sesh;

        subplot(nplots,1,sesh);
        hold on
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        % title([this_tp_title ' ITI std= ' num2str(round(std_curr_sesh,2)) ' skew=' num2str(round(skew_curr_sesh,4)) ' frac=' num2str(round(frac_past,4))]);
        title([this_tp_title '-DART ITI Bot75 std= ' num2str(round(bot20_std,2)) ' mean= ' num2str(round(bot20_avg,2))]);
        % yline(bot20_avg+3*bot20_std,"Color",'r');
        % yline(bot20_avg+5*bot20_std,"Color",'g');
        % yline(bot20_avg,"Color",'black');
        bin_over = NaN(nbins,1);
        for bin = 1:nbins
            this_piece = split_sesh{bin,1};
            this_sorted = sort(this_piece);
            this_bot = this_sorted(1:nelem*0.75);
            this_bot_avg = mean(this_bot);
            this_bot_std = std(this_bot);
            this_bot_thresh = this_bot_avg+5*this_bot_std;
            bin_over(bin) = sum(this_piece > this_bot_thresh);
            %plotting
            xstart = (bin-1)*nelem;
            xend = bin*nelem;
            plot([xstart xend],[this_bot_avg+3*this_bot_std this_bot_avg+3*this_bot_std],"Color",'r');
            plot([xstart xend],[this_bot_avg+5*this_bot_std this_bot_avg+5*this_bot_std],"Color",'g');
            plot([xstart xend],[this_bot_avg this_bot_avg],"Color",'black');
        end
        hold off
        frac_past = sum(bin_over)/length(curr_sesh);
        perc_mat(sesh) = frac_past;
    end
    std_cell_PV_ITI{mouse,1} = std_mat;
    frac_rms_PV_ITI{mouse,1} = perc_mat;
    burst_PV{mouse,1} = burst_mat;
    ent_PV{mouse,1} = ent_mat;
    % skew_cell_PV_ITI{mouse,1} = skew_mat;
    % saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC_ITI.pdf']));
    % saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC_ITI_bot50.pdf']));

    % PLOT DISTRIBUTION TO LOOK AT SKEWNESS
    % figure;
    % sgtitle(['Distribution ' mouse_id]);
    % for sesh = 1:size(mouse_data,1)
    %     curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
    %     nplots = 2; % hardcoded for 2 plots per mouse
    % 
    %     if sesh == 1
    %         this_tp_title = 'pre';
    %     else
    %         this_tp_title = 'post';
    %     end
    % 
    %     % only use ITI data
    %     data_trial = reshape(curr_data,nTot,[]);
    %     data_trial_ITI = data_trial(nOff/2+1:nOff,:);
    %     curr_sesh = reshape(data_trial_ITI,[],1);
    %     [N,edges] = histcounts(curr_sesh,100);
    %     skw = skew_cell_PV_ITI{mouse,1}(sesh,1);
    % 
    %     subplot(nplots,1,sesh);
    %     plot(N);
    %     title([this_tp_title '-DART ITI frames skew = ' num2str(round(skw,4))]);
    %     ylabel('counts')
    % end
    % saveas(gcf,fullfile(fn_out_root,[mouse_id '_distribution.pdf']));
end

% save(fullfile(fn_out_root,'stdPlusSkew_PV_ITI'),'std_cell_PV_ITI','skew_cell_PV_ITI','-v7.3');


%% CC
nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrial = 480;
nbins = 10;

ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

nMice = size(TCs_SST,1);
mice = TCs_SST(:,2);
std_cell_SST_ITI = cell(nMice,2);
std_cell_SST_ITI(:,2) = mice;
% skew_cell_SST_ITI = cell(nMice,2);
% skew_cell_SST_ITI(:,2) = mice;
frac_rms_SST_ITI = cell(nMice,2);
frac_rms_SST_ITI(:,2) = mice;
burst_SST = cell(nMice,2);
burst_SST(:,2) = mice;
ent_SST = cell(nMice,2);
ent_SST(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_SST{mouse,1};
    mouse_id = TCs_SST{mouse,2};
    nSesh = size(mouse_data,1);

    perc_mat = NaN(nSesh,1);
    std_mat = NaN(nSesh,1);
    burst_mat = NaN(nSesh,1);
    % ent_mat = NaN(nSesh,1);
    % skew_mat = NaN(nSesh,1);

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sesh = 1:size(mouse_data,1)
        curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        nplots = 2; % hardcoded for 2 plots per mouse

        if sesh == 1
            this_tp_title = 'pre';
        else
            this_tp_title = 'post';
        end

        % only use ITI data
        data_trial = reshape(curr_data,nTot,[]);
        data_trial_ITI = data_trial(nOff/2+1:nOff,:);

        curr_sesh = reshape(data_trial_ITI,[],1);
        % ent_mat(sesh) = SampleEn_TH(curr_sesh,0.2,3);

        this_burst = abs(diff(data_trial_ITI));
        mean_burst = mean(mean(this_burst));
        burst_mat(sesh) = mean_burst;

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        % skew_curr_sesh = skewness(curr_sesh);
        % frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        sorted = sort(curr_sesh);
        bot20 = sorted(1:length(curr_sesh)*0.75);
        bot20_avg = mean(bot20);
        bot20_std = std(bot20);
        rows_vector = repelem(length(curr_sesh)/nbins,nbins);
        nelem = length(curr_sesh)/nbins;
        split_sesh = mat2cell(curr_sesh,rows_vector);
        

        % perc_mat(sesh) = frac_past;
        std_mat(sesh) = std_curr_sesh;
        % skew_mat(sesh) = skew_curr_sesh;

        subplot(nplots,1,sesh);
        hold on
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        % title([this_tp_title ' ITI std= ' num2str(round(std_curr_sesh,2)) ' skew=' num2str(round(skew_curr_sesh,4)) ' frac=' num2str(round(frac_past,4))]);
        title([this_tp_title '-DART ITI Bot75 std= ' num2str(round(bot20_std,2)) ' mean= ' num2str(round(bot20_avg,2))]);
        % yline(bot20_avg+3*bot20_std,"Color",'r');
        % yline(bot20_avg+5*bot20_std,"Color",'g');
        % yline(bot20_avg,"Color",'black');
        bin_over = NaN(nbins,1);
        for bin = 1:nbins
            this_piece = split_sesh{bin,1};
            this_sorted = sort(this_piece);
            this_bot = this_sorted(1:nelem*0.75);
            this_bot_avg = mean(this_bot);
            this_bot_std = std(this_bot);
            this_bot_thresh = this_bot_avg+5*this_bot_std;
            bin_over(bin) = sum(this_piece > this_bot_thresh);
            %plotting
            xstart = (bin-1)*nelem;
            xend = bin*nelem;
            plot([xstart xend],[this_bot_avg+3*this_bot_std this_bot_avg+3*this_bot_std],"Color",'r');
            plot([xstart xend],[this_bot_avg+5*this_bot_std this_bot_avg+5*this_bot_std],"Color",'g');
            plot([xstart xend],[this_bot_avg this_bot_avg],"Color",'black');
        end
        hold off
        frac_past = sum(bin_over)/length(curr_sesh);
        perc_mat(sesh) = frac_past;
    end
    std_cell_SST_ITI{mouse,1} = std_mat;
    frac_rms_SST_ITI{mouse,1} = perc_mat;
    burst_SST{mouse,1} = burst_mat;
    % ent_SST{mouse,1} = ent_mat;
    % skew_cell_SST_ITI{mouse,1} = skew_mat;
    % saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC_ITI.pdf']));
    % saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC_ITI_bot50.pdf']));

    % PLOT DISTRIBUTION TO LOOK AT SKEWNESS
    % figure;
    % sgtitle(['Distribution ' mouse_id]);
    % for sesh = 1:size(mouse_data,1)
    %     curr_data = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
    %     nplots = 2; % hardcoded for 2 plots per mouse
    % 
    %     if sesh == 1
    %         this_tp_title = 'pre';
    %     else
    %         this_tp_title = 'post';
    %     end
    % 
    %     % only use ITI data
    %     data_trial = reshape(curr_data,nTot,[]);
    %     data_trial_ITI = data_trial(nOff/2+1:nOff,:);
    %     curr_sesh = reshape(data_trial_ITI,[],1);
    %     [N,edges] = histcounts(curr_sesh,100);
    %     skw = skew_cell_SST_ITI{mouse,1}(sesh,1);
    % 
    %     subplot(nplots,1,sesh);
    %     plot(N);
    %     title([this_tp_title '-DART ITI frames skew = ' num2str(round(skw,4))]);
    %     ylabel('counts')
    % end
    % saveas(gcf,fullfile(fn_out_root,[mouse_id '_distribution.pdf']));
end

% save(fullfile(fn_out_root,'stdPlusSkew_SST_ITI'),'std_cell_SST_ITI','skew_cell_SST_ITI','-v7.3');

%% plotting synchrony
if NoIntermediateFigures == 1
    close all
end

nMice_PV = size(iti_synch_PV,1);
nDay = 2;
start = 1;
SI_PV = NaN(nMice_PV,1);
for iM = 1:nMice_PV
    for iD = 1:nDay
        med_synch(iM,iD) = mean(iti_synch_PV{iM,iD});
    end
    SI_PV(iM) = (med_synch(iM,2)./med_synch(iM,1));
end

nMice_SST = size(iti_synch_SST,1);
nDay = 2;
start = 1;
SI_SST = NaN(nMice_SST,1);
for iM = 1:nMice_SST
    for iD = 1:nDay
        start = start+1;
        med_synch(iM,iD) = mean(iti_synch_SST{iM,iD});
    end
    SI_SST(iM) = (med_synch(iM,2)./med_synch(iM,1));
end

SI_PV_DART = vertcat(SI_PV(1:2),SI_PV(4));
SI_PV_PEG = vertcat(SI_PV(3),SI_PV(5));

% plotting 
f1 = figure;
f1.Name = 'Synchrony';
sgtitle('Synchrony');
hold on
plot(0.7,SI_PV_DART(:,1),'o','Color',"#0047AB")
plot(1.4,SI_PV_PEG(:,1),'o','Color',"#89CFF0")
plot(2.1,SI_SST(:,1),'o','Color',"#D95319")
ylim([0 3])
ylabel('Synchrony Index')
xlim([0 2.8])
xticks([0.7 1.4 2.1])
xticklabels({'PV+DART','PV+PEG','SST+DART'})
hold off


%% plot adjacent delta F

vector_df = cell2mat(burst_PV(:,1));
PV_df = reshape(vector_df,2,[])';
vector_df = cell2mat(burst_SST(:,1));
SST_df = reshape(vector_df,2,[])';
PV_dfIdx = PV_df(:,2) ./ PV_df(:,1);
SST_dfIdx = SST_df(:,2) ./ SST_df(:,1);

PV_dfIdx_dart = vertcat(PV_dfIdx(1:2),PV_dfIdx(4));
PV_dfIdx_peg = vertcat(PV_dfIdx(3),PV_dfIdx(5));

f2 = figure;
f2.Name = 'AdjacentDeltaF';
sgtitle('AdjacentDeltaF');
hold on
plot(0.7,PV_dfIdx_dart(:,1),'o','Color',"#0047AB")
plot(1.4,PV_dfIdx_peg(:,1),'o','Color',"#89CFF0")
plot(2.1,SST_dfIdx(:,1),'o','Color',"#D95319")
ylim([0 10])
ylabel('Post/Pre Trial Average Frame-by-Frame Delta F')
xlim([0 2.8])
xticks([0.7 1.4 2.1])
xticklabels({'PV+DART','PV+PEG','SST+DART'})
hold off


%% plotting something else
load('i3309_itiTC_pre.mat');
load('i3309_itiTC_post.mat');
ymax = max(vertcat(i3309_ITI_tc_pre,i3309_ITI_tc_post))+100;
ymin = min(vertcat(i3309_ITI_tc_pre,i3309_ITI_tc_post))-100;
f3 = figure;
f3.Name = 'i3309_ITI_TC';
sgtitle('i3309 ITI Timecourses')
hold on
subplot(2,1,1)
plot(i3309_ITI_tc_pre);
title('Control')
ylim([ymin ymax])
ylabel('F')
subplot(2,1,2)
plot(i3309_ITI_tc_post);
title('Post-DART')
ylim([ymin ymax])
ylabel('F')
hold off


i3309_iti_sync = [iti_synch_PV{1,1} iti_synch_PV{1,2}];
nDay = 2;
f4 = figure;
f4.Name = 'Synch_hist';
sgtitle('i3309 Percent Active Histogram');
hold on
subplot(2,1,1)
histogram(i3309_iti_sync(:,1),[0:0.01:1]);
xlim([0 1])
ylim([0 100])
title('Baseline')
ylabel('nFrames')
subplot(2,1,2)
histogram(i3309_iti_sync(:,2),[0:0.01:1]);
xlim([0 1])
ylim([0 100])
title('Post-DART')
ylabel('nFrames')
xlabel('Percent Active')
hold off

% save_all_open_figs('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\finalPlots');
%% test

x = curr_sesh;         % example time series
p = 2;                    % lag order
includeTrend = 0;     % no deterministic trend

[t_stat, crit_vals, result] = adftest_TH(x, p, includeTrend);

% Assume your time series is in a vector called 'y'
autocorr_TH(diff(curr_sesh), 40);
title('Autocorrelation Function (ACF)');

y = curr_sesh;
t = (1:length(y))';  % Time vector
mdl = fitlm(t, y);  % Fit a linear model to the data
trend = mdl.Fitted;  % Get the fitted trend
y_detrended = y - trend;  % Remove the trend

figure
plot(y_detrended)

figure
hold on
plot(y,"Color",'b')
plot(trend,"Color",'black')
hold off