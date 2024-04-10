close all;
clear all;

datafile = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/MIPs.mat';
outpath = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks';
mask_fig_path = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks';
mask_fig_path_vSize = 'G:/home/ACh/Analysis/histology_quant/IHC03_confocal_TH/masks/mask_pngs/px300';
analysis_out = fullfile(mask_fig_path_vSize,'analysis');

mkdir(analysis_out);
load(fullfile(mask_fig_path_vSize,'final_masks_300.mat'));
load(fullfile(mask_fig_path_vSize,'cell_counts_300.mat'));
load(datafile);

%% extract time courses
% Two scenarios:
% tc_other - ch1 mask on ch2 & ch2 mask on ch1
% tc_self - mask applied to the channel it was segmented from
% all_ttests - ch1 mask on ch2 vs ch1 mask rotated cc 90 degrees on ch2


% tc_other = cell(23,3);
% tc_self = cell(23,3);
% all_ttests = zeros(23,2);
% 
nSer = size(final_masks,1);
nCh = size(final_masks,2)-1;
% 
% for ser = 1:nSer
%     for ch = 1:nCh
%         if ch == 1
%             extract_img = 2;
%         else
%             extract_img = 1;
%         end
%         % gets intensity of each object in THE OTHER CHANNEL using the mask
%         % obtained by segmenting the current channel 
%         otherCh_intensity = stackGetTimeCourses(img_data{ser,extract_img},final_masks{ser,ch});
%         tc_other{ser,extract_img} = otherCh_intensity;
%         maskCh_intensity = stackGetTimeCourses(img_data{ser,ch},final_masks{ser,ch});
%         tc_self{ser,ch} = maskCh_intensity;
%         % rotates the mask counterclock wise 90 degrees 3 times
%         mask_rot = rot90(final_masks{ser,ch});
%         control_intensity = stackGetTimeCourses(img_data{ser,extract_img},mask_rot);
%         all_ttests(ser,extract_img) = ttest(otherCh_intensity,control_intensity);
%     end
%     tc_other{ser,3} = img_data{ser,3}; 
% end


% save(fullfile(analysis_out,'tc_other'),'tc_other');
% save(fullfile(analysis_out,'tc_self'),'tc_self');
% save(fullfile(analysis_out,'ttests'),'all_ttests');
clear img_data;
%% prep data
load(fullfile(analysis_out,'tc_other'),'tc_other');
load(fullfile(analysis_out,'tc_self'),'tc_self');
load(fullfile(analysis_out,'ttests'),'all_ttests');

% excluding these image series due to overlapping FOV
excluded_rows = [2 9];
masks = final_masks;
masks(excluded_rows,:) = [];
tc_other(excluded_rows,:) = [];
tc_self(excluded_rows,:) = [];

nSer = size(masks,1);
nCh = size(masks,2)-1;
%% plots

% scatter plot - avg channel intensity by condition and channel

ChAvg_other = zeros(nSer,nCh);
ChAvg_self = zeros(nSer,nCh);

for ser = 1:nSer
    for ch = 1:nCh
        if ch == 1
            extract_img = 2;
        else
            extract_img = 1;
        end
        ChAvg_other(ser,ch) = mean(tc_other{ser,extract_img});
        ChAvg_self(ser,ch) = mean(tc_self{ser,ch});
    end
end

% other_table = cell2table(ChAvg_other,"VariableNames",["ch1" "ch2" "sample"]);
% self_table = cell2table(ChAvg_self);

% check if image series contains certain strings

series = string(tc_other(:,3));
mb_names = ["TH040","TH041"];
nes_names = ["TH042","TH043"];
mSc_names = ["TH044","TH045"];

mb_idx = find(contains(series,mb_names));
nes_idx = find(contains(series,nes_names));
mSc_idx = find(contains(series,mSc_names));



% self
close all;
figure;
subplot(1,3,1);
sgtitle('mask channel intensities');
hold on
title('MB HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_self(mb_idx,1),'black');
scatter(0.66,ChAvg_self(mb_idx,2),'black');
avg_mb_1 = mean(ChAvg_self(mb_idx,1));
sem_mb_1 = std(ChAvg_self(mb_idx,1))/sqrt(length(mb_idx));
avg_mb_2 = mean(ChAvg_self(mb_idx,2));
sem_mb_2 = std(ChAvg_self(mb_idx,2))/sqrt(length(mb_idx));
plot(0.33,avg_mb_1,'+r');
errorbar(0.33,avg_mb_1,sem_mb_1,"o");
plot(0.66,avg_mb_2,'+r');
errorbar(0.66,avg_mb_2,sem_mb_2,"o");

subplot(1,3,2);
hold on
title('NES HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_self(nes_idx,1),'black');
scatter(0.66,ChAvg_self(nes_idx,2),'black');
avg_nes_1 = mean(ChAvg_self(nes_idx,1));
sem_nes_1 = std(ChAvg_self(nes_idx,1))/sqrt(length(nes_idx));
avg_nes_2 = mean(ChAvg_self(nes_idx,2));
sem_nes_2 = std(ChAvg_self(nes_idx,2))/sqrt(length(nes_idx));
plot(0.33,avg_nes_1,'+r');
errorbar(0.33,avg_nes_1,sem_nes_1,"o");
plot(0.66,avg_nes_2,'+r');
errorbar(0.66,avg_nes_2,sem_nes_2,"o");

subplot(1,3,3);
hold on
title('mScarlett HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_self(mSc_idx,1),'black');
scatter(0.66,ChAvg_self(mSc_idx,2),'black');
avg_mSc_1 = mean(ChAvg_self(mSc_idx,1));
sem_mSc_1 = std(ChAvg_self(mSc_idx,1))/sqrt(length(mSc_idx));
avg_mSc_2 = mean(ChAvg_self(mSc_idx,2));
sem_mSc_2 = std(ChAvg_self(mSc_idx,2))/sqrt(length(mSc_idx));
plot(0.33,avg_mSc_1,'+r');
errorbar(0.33,avg_mSc_1,sem_mSc_1,"o");
plot(0.66,avg_mSc_2,'+r');
errorbar(0.66,avg_mSc_2,sem_mSc_2,"o");

% other
figure;
subplot(1,3,1);
sgtitle('other channel intensities');
hold on
title('MB HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_other(mb_idx,1),'black');
scatter(0.66,ChAvg_other(mb_idx,2),'black');
avg_mb_1 = mean(ChAvg_other(mb_idx,1));
sem_mb_1 = std(ChAvg_other(mb_idx,1))/sqrt(length(mb_idx));
avg_mb_2 = mean(ChAvg_other(mb_idx,2));
sem_mb_2 = std(ChAvg_other(mb_idx,1))/sqrt(length(mb_idx));
plot(0.33,avg_mb_1,'+r');
errorbar(0.33,avg_mb_1,sem_mb_1,"o");
plot(0.66,avg_mb_2,'+r');
errorbar(0.66,avg_mb_2,sem_mb_2,"o");

subplot(1,3,2);
hold on
title('NES HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_other(nes_idx,1),'black');
scatter(0.66,ChAvg_other(nes_idx,2),'black');
avg_nes_1 = mean(ChAvg_other(nes_idx,1));
sem_nes_1 = std(ChAvg_other(nes_idx,1))/sqrt(length(nes_idx));
avg_nes_2 = mean(ChAvg_other(nes_idx,2));
sem_nes_2 = std(ChAvg_other(nes_idx,1))/sqrt(length(nes_idx));
plot(0.33,avg_nes_1,'+r');
errorbar(0.33,avg_nes_1,sem_nes_1,"o");
plot(0.66,avg_nes_2,'+r');
errorbar(0.66,avg_nes_2,sem_nes_2,"o");

subplot(1,3,3);
hold on
title('mScarlett HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,ChAvg_other(mSc_idx,1),'black');
scatter(0.66,ChAvg_other(mSc_idx,2),'black');
avg_mSc_1 = mean(ChAvg_other(mSc_idx,1));
sem_mSc_1 = std(ChAvg_other(mSc_idx,1))/sqrt(length(mSc_idx));
avg_mSc_2 = mean(ChAvg_other(mSc_idx,2));
sem_mSc_2 = std(ChAvg_other(mSc_idx,1))/sqrt(length(mSc_idx));
plot(0.33,avg_mSc_1,'+r');
errorbar(0.33,avg_mSc_1,sem_mSc_1,"o");
plot(0.66,avg_mSc_2,'+r');
errorbar(0.66,avg_mSc_2,sem_mSc_2,"o");

hold off
%% normalized

% "other channel" mask intensities/75% of the minimum "mask channel" mask intensities

intensity_cell = tc_self(:,[1 2]);
baseline = cellfun(@min,intensity_cell) * 0.75;

% since tc_other & tc_self organization is based on mask origin, intensities for
% Ch2 on tc_other is stored in column 1 and vice versa. Here, invert before
% normalization 

inv_tc_other(:,1) = tc_other(:,2);
inv_tc_other(:,2) = tc_other(:,1);
inv_tc_other(:,3) = tc_other(:,3);

norm_tc_other = cell(nSer,nCh+1);

for ser = 1:size(norm_tc_other,1)
    for ch = 1:nCh
        norm_tc_other{ser,ch} = inv_tc_other{ser,ch} / baseline(ser,ch);
    end
end

norm_tc_other(:,3) = inv_tc_other(:,3);

% now normalize the self tc

norm_tc_self = cell(nSer,nCh+1);

for ser = 1:size(norm_tc_self,1)
    for ch = 1:nCh
        norm_tc_self{ser,ch} = tc_self{ser,ch} / baseline(ser,ch);
    end
end

norm_tc_self(:,3) = tc_other(:,3);

%% plot the normalized version

norm_ChAvg_other = zeros(nSer,nCh);
norm_ChAvg_self = zeros(nSer,nCh);

for ser = 1:nSer
    for ch = 1:nCh
        norm_ChAvg_other(ser,ch) = mean(norm_tc_other{ser,ch});
        norm_ChAvg_self(ser,ch) = mean(norm_tc_self{ser,ch});
    end
end


% self

figure;
subplot(1,3,1);
sgtitle('normalized mask channel intensities');
hold on
title('MB HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
ylabel('mean ch intensities / 75% of min obj intensity');
scatter(0.33,norm_ChAvg_self(mb_idx,1),'black');
scatter(0.66,norm_ChAvg_self(mb_idx,2),'black');
avg_mb_1 = mean(norm_ChAvg_self(mb_idx,1));
avg_mb_2 = mean(norm_ChAvg_self(mb_idx,2));
plot(0.33,avg_mb_1,'+r');
plot(0.66,avg_mb_2,'+r');

subplot(1,3,2);
hold on
title('NES HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,norm_ChAvg_self(nes_idx,1),'black');
scatter(0.66,norm_ChAvg_self(nes_idx,2),'black');
avg_nes_1 = mean(norm_ChAvg_self(nes_idx,1));
avg_nes_2 = mean(norm_ChAvg_self(nes_idx,2));
plot(0.33,avg_nes_1,'+r');
plot(0.66,avg_nes_2,'+r');

subplot(1,3,3);
hold on
title('mScarlett HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,norm_ChAvg_self(mSc_idx,1),'black');
scatter(0.66,norm_ChAvg_self(mSc_idx,2),'black');
avg_mSc_1 = mean(norm_ChAvg_self(mSc_idx,1));
avg_mSc_2 = mean(norm_ChAvg_self(mSc_idx,2));
plot(0.33,avg_mSc_1,'+r');
plot(0.66,avg_mSc_2,'+r');

% other
figure;
subplot(1,3,1);
sgtitle('normalized other channel intensities');
hold on
title('MB HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
ylabel('mean ch intensities / 75% of min obj intensity');
scatter(0.33,norm_ChAvg_other(mb_idx,1),'black');
scatter(0.66,norm_ChAvg_other(mb_idx,2),'black');
avg_mb_1 = mean(norm_ChAvg_other(mb_idx,1));
avg_mb_2 = mean(norm_ChAvg_other(mb_idx,2));
plot(0.33,avg_mb_1,'+r');
plot(0.66,avg_mb_2,'+r');

subplot(1,3,2);
hold on
title('NES HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,norm_ChAvg_other(nes_idx,1),'black');
scatter(0.66,norm_ChAvg_other(nes_idx,2),'black');
avg_nes_1 = mean(norm_ChAvg_other(nes_idx,1));
avg_nes_2 = mean(norm_ChAvg_other(nes_idx,2));
plot(0.33,avg_nes_1,'+r');
plot(0.66,avg_nes_2,'+r');

subplot(1,3,3);
hold on
title('mScarlett HTP');
xticks([0.33 0.66]);
xlim([0 1]);
xticklabels({'r','c'});
scatter(0.33,norm_ChAvg_other(mSc_idx,1),'black');
scatter(0.66,norm_ChAvg_other(mSc_idx,2),'black');
avg_mSc_1 = mean(norm_ChAvg_other(mSc_idx,1));
avg_mSc_2 = mean(norm_ChAvg_other(mSc_idx,2));
plot(0.33,avg_mSc_1,'+r');
plot(0.66,avg_mSc_2,'+r');