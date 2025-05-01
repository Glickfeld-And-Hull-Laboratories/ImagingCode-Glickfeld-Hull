clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories

fn_analysis_root = 'G:\home\ACh\Analysis\2p_analysis';
fn_data_root = 'G:\home\ACh\Data\2p_data';
fn_out_root = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask';
cd(fn_out_root)

load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask\TCs_SST.mat');
load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask\TCs_PV.mat');

%% TH

ds = 'DART_expt_info';
eval(ds);

nMice = size(TCs_PV,1);
mice = TCs_PV(:,2);
std_cell_PV = cell(nMice,2);
std_cell_PV(:,2) = mice;
skew_cell_PV = cell(nMice,2);
skew_cell_PV(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_PV{mouse,1};
    mouse_id = TCs_PV{mouse,2};
    nSesh = size(mouse_data,1);

    perc_mat = NaN(nSesh,1);
    std_mat = NaN(nSesh,1);
    skew_mat = NaN(nSesh,1);

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sesh = 1:size(mouse_data,1)
        curr_sesh = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        nplots = 2; % hardcoded for 2 plots per mouse

        if sesh == 1
            this_tp_title = 'pre';
        else
            this_tp_title = 'post';
        end

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        skew_curr_sesh = skewness(curr_sesh);
        frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        perc_mat(sesh) = frac_past;
        std_mat(sesh) = std_curr_sesh;
        skew_mat(sesh) = skew_curr_sesh;

        subplot(nplots,1,sesh);
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        title([this_tp_title ' std= ' num2str(round(std_curr_sesh,2)) ' skew=' num2str(round(skew_curr_sesh,4)) ' frac=' num2str(round(frac_past,4))]);
    end
    std_cell_PV{mouse,1} = std_mat;
    skew_cell_PV{mouse,1} = skew_mat;
    saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC.pdf']));
end

% save(fullfile(fn_out_root,'stdPlusSkew_PV'),'std_cell_PV','skew_cell_PV','-v7.3');

%% CC

ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

nMice = size(TCs_SST,1);
mice = TCs_SST(:,2);
std_cell_SST = cell(nMice,2);
std_cell_SST(:,2) = mice;
skew_cell_SST = cell(nMice,2);
skew_cell_SST(:,2) = mice;

for mouse = 1:nMice %iterate through mouse
    mouse_data = TCs_SST{mouse,1};
    mouse_id = TCs_SST{mouse,2};
    nSesh = size(mouse_data,1);

    perc_mat = NaN(nSesh,1);
    std_mat = NaN(nSesh,1);
    skew_mat = NaN(nSesh,1);

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sesh = 1:size(mouse_data,1)
        curr_sesh = mouse_data{sesh,1}; % data for the current imaging session of the current mouse
        nplots = 2; % hardcoded for 2 plots per mouse

        if sesh == 1
            this_tp_title = 'pre';
        else
            this_tp_title = 'post';
        end

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        skew_curr_sesh = skewness(curr_sesh);
        frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        perc_mat(sesh) = frac_past;
        std_mat(sesh) = std_curr_sesh;
        skew_mat(sesh) = skew_curr_sesh;

        subplot(nplots,1,sesh);
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        title([this_tp_title ' std= ' num2str(round(std_curr_sesh,2)) ' skew=' num2str(round(skew_curr_sesh,4)) ' frac=' num2str(round(frac_past,4))]);
    end
    std_cell_SST{mouse,1} = std_mat;
    skew_cell_SST{mouse,1} = skew_mat;
    saveas(gcf,fullfile(fn_out_root,[mouse_id 'combinedTC.pdf']));
end

% save(fullfile(fn_out_root,'stdPlusSkew_SST'),'std_cell_SST','skew_cell_SST','-v7.3');