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
% data in corresponding columns indicate MASK ORIGIN


% tc_other = cell(23,3);
% tc_self = cell(23,3);
% ctrl_tc_other = cell(23,3);
% all_ttests = zeros(23,2);
% 
% nSer = size(final_masks,1);
% nCh = size(final_masks,2)-1;

% for ser = 1:nSer
%     for ch = 1:nCh
%         if ch == 1
%             extract_img = 2;
%         else
%             extract_img = 1;
%         end
%         % gets intensity of each object in THE OTHER CHANNEL using the mask
%         % obtained by segmenting the current channel 
%         otherCh_intensity = stackGetTimeCourses(img_data{ser,ch},final_masks{ser,extract_img});
%         tc_other{ser,extract_img} = otherCh_intensity;
%         maskCh_intensity = stackGetTimeCourses(img_data{ser,ch},final_masks{ser,ch});
%         tc_self{ser,ch} = maskCh_intensity;
%         % rotates the mask counterclock wise 90 degrees 3 times
%         mask_rot = rot90(final_masks{ser,extract_img});
%         control_intensity = stackGetTimeCourses(img_data{ser,ch},mask_rot);
%         ctrl_tc_other{ser,extract_img} = control_intensity;
%         all_ttests(ser,extract_img) = ttest(otherCh_intensity,control_intensity);
%     end
%     tc_self{ser,3} = img_data{ser,3}; 
%     tc_other{ser,3} = img_data{ser,3}; 
%     ctrl_tc_other{ser,3} = img_data{ser,3};
% end


% save(fullfile(analysis_out,'tc_other'),'tc_other');
% save(fullfile(analysis_out,'tc_self'),'tc_self');
% save(fullfile(analysis_out,'ctrl_tc_other'),'ctrl_tc_other');
% save(fullfile(analysis_out,'ttests'),'all_ttests');
clear img_data;
%% prep data
load(fullfile(analysis_out,'tc_other'),'tc_other');
load(fullfile(analysis_out,'tc_self'),'tc_self');
load(fullfile(analysis_out,'ctrl_tc_other'),'ctrl_tc_other');
load(fullfile(analysis_out,'ttests'),'all_ttests');

% excluding these image series due to overlapping FOV
excluded_rows = [2 9];
masks = final_masks;
masks(excluded_rows,:) = [];
tc_other(excluded_rows,:) = [];
tc_self(excluded_rows,:) = [];
ctrl_tc_other(excluded_rows,:) = [];
all_ttests(excluded_rows,:) = [];

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

virus_indices = cell(3,2);
mb_idx = find(contains(series,mb_names));
virus_indices{1,1} = mb_idx;
virus_indices{1,2} = 'mb';
nes_idx = find(contains(series,nes_names));
virus_indices{2,1} = nes_idx;
virus_indices{2,2} = 'nes';
mSc_idx = find(contains(series,mSc_names));
virus_indices{3,1} = mSc_idx;
virus_indices{3,2} = 'mSc';


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
% deprecated "normalization"
% % "other channel" mask intensities/75% of the minimum "mask channel" mask intensities
% 
% intensity_cell = tc_self(:,[1 2]);
% baseline = cellfun(@min,intensity_cell) * 0.75;
% 
% % since tc_other & tc_self organization is based on mask origin, intensities for
% % Ch2 on tc_other is stored in column 1 and vice versa. Here, invert before
% % normalization 
% 
% inv_tc_other(:,1) = tc_other(:,2);
% inv_tc_other(:,2) = tc_other(:,1);
% inv_tc_other(:,3) = tc_other(:,3);
% 
% norm_tc_other = cell(nSer,nCh+1);
% 
% for ser = 1:size(norm_tc_other,1)
%     for ch = 1:nCh
%         norm_tc_other{ser,ch} = inv_tc_other{ser,ch} / baseline(ser,ch);
%     end
% end
% 
% norm_tc_other(:,3) = inv_tc_other(:,3);
% 
% % now normalize the self tc
% 
% norm_tc_self = cell(nSer,nCh+1);
% 
% for ser = 1:size(norm_tc_self,1)
%     for ch = 1:nCh
%         norm_tc_self{ser,ch} = tc_self{ser,ch} / baseline(ser,ch);
%     end
% end
% 
% norm_tc_self(:,3) = tc_other(:,3);
% 
% %% plot the normalized version
% 
% norm_ChAvg_other = zeros(nSer,nCh);
% norm_ChAvg_self = zeros(nSer,nCh);
% 
% for ser = 1:nSer
%     for ch = 1:nCh
%         norm_ChAvg_other(ser,ch) = mean(norm_tc_other{ser,ch});
%         norm_ChAvg_self(ser,ch) = mean(norm_tc_self{ser,ch});
%     end
% end
% 
% 
% % self
% 
% figure;
% subplot(1,3,1);
% sgtitle('normalized mask channel intensities');
% hold on
% title('MB HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% ylabel('mean ch intensities / 75% of min obj intensity');
% scatter(0.33,norm_ChAvg_self(mb_idx,1),'black');
% scatter(0.66,norm_ChAvg_self(mb_idx,2),'black');
% avg_mb_1 = mean(norm_ChAvg_self(mb_idx,1));
% avg_mb_2 = mean(norm_ChAvg_self(mb_idx,2));
% plot(0.33,avg_mb_1,'+r');
% plot(0.66,avg_mb_2,'+r');
% 
% subplot(1,3,2);
% hold on
% title('NES HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% scatter(0.33,norm_ChAvg_self(nes_idx,1),'black');
% scatter(0.66,norm_ChAvg_self(nes_idx,2),'black');
% avg_nes_1 = mean(norm_ChAvg_self(nes_idx,1));
% avg_nes_2 = mean(norm_ChAvg_self(nes_idx,2));
% plot(0.33,avg_nes_1,'+r');
% plot(0.66,avg_nes_2,'+r');
% 
% subplot(1,3,3);
% hold on
% title('mScarlett HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% scatter(0.33,norm_ChAvg_self(mSc_idx,1),'black');
% scatter(0.66,norm_ChAvg_self(mSc_idx,2),'black');
% avg_mSc_1 = mean(norm_ChAvg_self(mSc_idx,1));
% avg_mSc_2 = mean(norm_ChAvg_self(mSc_idx,2));
% plot(0.33,avg_mSc_1,'+r');
% plot(0.66,avg_mSc_2,'+r');
% 
% % other
% figure;
% subplot(1,3,1);
% sgtitle('normalized other channel intensities');
% hold on
% title('MB HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% ylabel('mean ch intensities / 75% of min obj intensity');
% scatter(0.33,norm_ChAvg_other(mb_idx,1),'black');
% scatter(0.66,norm_ChAvg_other(mb_idx,2),'black');
% avg_mb_1 = mean(norm_ChAvg_other(mb_idx,1));
% avg_mb_2 = mean(norm_ChAvg_other(mb_idx,2));
% plot(0.33,avg_mb_1,'+r');
% plot(0.66,avg_mb_2,'+r');
% 
% subplot(1,3,2);
% hold on
% title('NES HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% scatter(0.33,norm_ChAvg_other(nes_idx,1),'black');
% scatter(0.66,norm_ChAvg_other(nes_idx,2),'black');
% avg_nes_1 = mean(norm_ChAvg_other(nes_idx,1));
% avg_nes_2 = mean(norm_ChAvg_other(nes_idx,2));
% plot(0.33,avg_nes_1,'+r');
% plot(0.66,avg_nes_2,'+r');
% 
% subplot(1,3,3);
% hold on
% title('mScarlett HTP');
% xticks([0.33 0.66]);
% xlim([0 1]);
% xticklabels({'r','c'});
% scatter(0.33,norm_ChAvg_other(mSc_idx,1),'black');
% scatter(0.66,norm_ChAvg_other(mSc_idx,2),'black');
% avg_mSc_1 = mean(norm_ChAvg_other(mSc_idx,1));
% avg_mSc_2 = mean(norm_ChAvg_other(mSc_idx,2));
% plot(0.33,avg_mSc_1,'+r');
% plot(0.66,avg_mSc_2,'+r');



%% distribution and inclusion criteria

% visualizing intensity distributions

[plotx, ploty] = subplotn(nSer*nCh);

figure();
hold on;
sgtitle('Distribution of self mask intensities vs other mask intensities');
for ser = 1:nSer
    for ch = 1:nCh
        if ch == 1
            extract = 2;
        else
            extract = 1;
        end
        bin_n1 = round(size(tc_self{ser,ch},2)/2);
        bin_n2 = round(size(tc_other{ser,extract},2)/2);
        current_plotn = (ser-1)*2 + ch;
        subplot(plotx,ploty,current_plotn);
        hold on;
        histogram(tc_self{ser,ch},bin_n1);
        histogram(tc_other{ser,extract},bin_n2);
        title_str = ['Series ' num2str(ser) ' CH ' num2str(ch) ' as baseline'];
        title(title_str);
        hold off;
    end
end
hold off;

% saveas(gcf,fullfile(analysis_out,'intensity_distributions'));

% calculations
coloc_count = zeros(nSer,nCh); % records # of cells above threshold for other mask
coloc_percent = zeros(nSer,nCh);

for ser = 1:nSer
    for ch = 1:nCh
        if ch == 1
            extract = 2;
        else
            extract = 1;
        end
        threshold = min(tc_self{ser,ch})*0.75;
        %tc_self and tc_other are organized by the origin of the MASK
        this_other_mask_int = tc_other{ser,extract};
        this_other_mask_int(this_other_mask_int<threshold) = 0;
        this_cell_count = length(find(this_other_mask_int));
        coloc_count(ser,ch) = this_cell_count;
        coloc_percent(ser,ch) = this_cell_count/length(tc_other{ser,extract});
    end
end

% save(fullfile(analysis_out,'coloc_count'),"coloc_count");
% save(fullfile(analysis_out,'coloc_percent'),"coloc_percent");

% tb_count = array2table(coloc_count);
% tb_percent = array2table(coloc_percent);
% 
% tb_count.Properties.VariableNames("coloc_count1") = "Red";
% tb_count.Properties.VariableNames("coloc_count2") = "Cyan";
% tb_percent.Properties.VariableNames("coloc_percent1") = "Red";
% tb_percent.Properties.VariableNames("coloc_percent2") = "Cyan";

%graph these as boxplots etc. 
figure();
sgtitle('75% as threshold')
subplot(1,2,1);
title('nCells passing threshold');
hold on;
boxplot(coloc_count(:,1));
boxplot(coloc_count(:,2),'Positions',2);
scatter(1,coloc_count(:,1),20,'black');
scatter(2,coloc_count(:,2),20,'black');
xtl1 = '\begin{tabular}{c} SST+ objects\\ also HTP+\end{tabular}';
xtl2 = '\begin{tabular}{c} HTP+ objects\\ also SST+\end{tabular}';
set(gca(),'XTick',[1 2],'XTickLabels',{xtl1,xtl2},'TickLabelInterpreter', 'latex');
xlim([0.33 2.66]);
hold off;

subplot(1,2,2);
title('% passing threshold');
hold on;
boxplot(coloc_percent(:,1));
boxplot(coloc_percent(:,2),'Positions',2);
scatter(1,coloc_percent(:,1),20,'black');
scatter(2,coloc_percent(:,2),20,'black');
xtl1 = '\begin{tabular}{c} SST+ objects\\ also HTP+\end{tabular}';
xtl2 = '\begin{tabular}{c} HTP+ objects\\ also SST+\end{tabular}';
set(gca(),'XTick',[1 2],'XTickLabels',{xtl1,xtl2},'TickLabelInterpreter', 'latex');
xlim([0.33 2.66]);
hold off;

%% rotated mask comparisons

rot_count = zeros(nSer,nCh); % records # of cells above threshold for other mask
rot_percent = zeros(nSer,nCh);
thresholds = zeros(nSer,nCh);

for ser = 1:nSer
    for ch = 1:nCh
        if ch == 1
            extract = 2;
        else
            extract = 1;
        end
        threshold = min(tc_self{ser,ch})*0.75;
        thresholds(ser,ch) = threshold;
        %tc_self and tc_other are organized by the origin of the MASK
        rot_mask_int = ctrl_tc_other{ser,extract};
        rot_mask_int(rot_mask_int<threshold) = 0;
        this_cell_count = length(find(rot_mask_int));
        rot_count(ser,ch) = this_cell_count;
        rot_percent(ser,ch) = this_cell_count/length(ctrl_tc_other{ser,extract});
    end
end

figure();
sgtitle('75% as threshold, rotated control')
hold on;
boxplot(rot_percent(:,1));
boxplot(rot_percent(:,2),'Positions',2);
scatter(1,rot_percent(:,1),15,'black');
scatter(2,rot_percent(:,2),15,'black');
set(gca(),'XTick',[1 2],'XTickLabels',{'SST+ objects also HTP+','HTP+ objects also SST+'})
ylabel('% passing threshold')
xlim([0.33 2.66]);
hold off;


%% plot within virus
% color codes r,g,b = MB, NES, mSc

colors = [1,0,0;0.4660 0.6740 0.1880;0,0,1];
figure();
sgtitle('75% as threshold')
% cell counts subplot
subplot(1,2,1);
title('nCells passing threshold');
hold on;
boxplot(coloc_count(:,1));
boxplot(coloc_count(:,2),'Positions',2);

for virus = 1:size(virus_indices,1) % order: MB, NES, mScarlett
    scatter(1,coloc_count(virus_indices{virus,1},1),20,colors(virus,:));
end
for virus = 1:size(virus_indices,1) % order: MB, NES, mScarlett
    scatter(2,coloc_count(virus_indices{virus,1},2),20,colors(virus,:));
end

% scatter(1,coloc_count(:,1),20,'black');
% scatter(2,coloc_count(:,2),20,'black');
xtl1 = '\begin{tabular}{c} SST+ objects\\ also HTP+\end{tabular}';
xtl2 = '\begin{tabular}{c} HTP+ objects\\ also SST+\end{tabular}';
set(gca(),'XTick',[1 2],'XTickLabels',{xtl1,xtl2},'TickLabelInterpreter', 'latex');
xlim([0.33 2.66]);
hold off;


% percentage subplot

subplot(1,2,2);
title('% passing threshold');
hold on;
boxplot(coloc_percent(:,1));
boxplot(coloc_percent(:,2),'Positions',2);

for virus = 1:size(virus_indices,1) % order: MB, NES, mScarlett
    scatter(1,coloc_percent(virus_indices{virus,1},1),20,colors(virus,:));
end
for virus = 1:size(virus_indices,1) % order: MB, NES, mScarlett
    scatter(2,coloc_percent(virus_indices{virus,1},2),20,colors(virus,:));
end

avgs = mean(coloc_percent);
plot(1,avgs(1),'+r','Color','black','Markersize',20);
plot(2,avgs(2),'+r','Color','black','Markersize',20);
% scatter(1,coloc_percent(:,1),20,'black');
% scatter(2,coloc_percent(:,2),20,'black');

xtl1 = '\begin{tabular}{c} SST+ objects\\ also HTP+\end{tabular}';
xtl2 = '\begin{tabular}{c} HTP+ objects\\ also SST+\end{tabular}';
set(gca(),'XTick',[1 2],'XTickLabels',{xtl1,xtl2},'TickLabelInterpreter', 'latex');
xlim([0.33 2.66]);
hold off;

%% mask size

m_sizes = cell(nSer,nCh);
avg_m_sizes = zeros(nSer,nCh);

for ser = 1:nSer
    for ch = 1:nCh
        edges = unique(masks{ser,ch});
        real_edges = [edges' max(edges+1)];
        counts_w_zero = histcounts(masks{ser,ch}(:),real_edges);
        counts = counts_w_zero(2:end);
        m_sizes{ser,ch} = counts;
        avg_m_sizes(ser,ch) = mean(counts);
    end
end

m_sizes_mean = mean(avg_m_sizes);

figure;
% sgtitle('mask sizes');
hold on
title('average mask object sizes');
xlim([0.5 2.5]);
ylabel('pixels');
scatter(1,avg_m_sizes(:,1),'black');
scatter(2,avg_m_sizes(:,2),'black');
plot(1,m_sizes_mean(1),'+r','MarkerSize',20);
% errorbar(0.33,avg_mb_1,sem_mb_1,"o");
plot(2,m_sizes_mean(2),'+r','MarkerSize',20);
% errorbar(0.66,avg_mb_2,sem_mb_2,"o");
set(gca(),'XTick',[1 2],'XTickLabels',{'red','cyan'});
hold off

%% testing for colored scatterplot
% data1 = rand(100, 2); % Assuming you have 100 data points with 2 features
% categories = randi([1, 5], 100, 1); % Assuming you have 5 categories
% 
% % Define colors for each category
% colors = hsv(max(categories)); % Using HSV colormap for variety
% 
% % Plot each category separately
% figure();
% hold on;
% for i = 1:max(categories)
%     indices = categories == i;
%     scatter(data1(indices, 1), data1(indices, 2), [], colors(i,:), 'filled');
% end
% hold off;
% 
% % Add legend
% legend('Category 1', 'Category 2', 'Category 3', 'Category 4', 'Category 5');
% 
% % Add labels
% xlabel('Feature 1');
% ylabel('Feature 2');
% title('Scatter Plot by Category');


%% testing for mask boundaries
load(datafile);
img_data(excluded_rows,:) = [];

overlay_mask_path = fullfile(analysis_out,'mask_overlays');

% close all
% figure;
% imshow(img_data{1,1});
% hold on
% bound = cell2mat(bwboundaries(masks{1,1}));
% plot(bound(:,2),bound(:,1),'.','Color','red','MarkerSize',1);


% red masks overlay both channels
% for ser = 1:nSer
%     this_bound = cell2mat(bwboundaries(masks{ser,1})); %only using red masks
%     for ch = 1:nCh
%         imshow(img_data{ser,ch});
%         hold on
%         plot(this_bound(:,2),this_bound(:,1),'.','Color','red','MarkerSize',1);
%         saveas(gcf,fullfile(overlay_mask_path,['RedMaskOverlay_' 'ser' num2str(ser) '_ch' num2str(ch) '.png']));
%         close;
%     end
% end
% 
% % cyan masks overlay both channels
% for ser = 1:nSer
%     this_bound = cell2mat(bwboundaries(masks{ser,2})); %only using red masks
%     for ch = 1:nCh
%         imshow(img_data{ser,ch});
%         hold on
%         plot(this_bound(:,2),this_bound(:,1),'.','Color','cyan','MarkerSize',1);
%         saveas(gcf,fullfile(overlay_mask_path,['CyanMaskOverlay_' 'ser' num2str(ser) '_ch' num2str(ch) '.png']));
%         close;
%     end
% end
