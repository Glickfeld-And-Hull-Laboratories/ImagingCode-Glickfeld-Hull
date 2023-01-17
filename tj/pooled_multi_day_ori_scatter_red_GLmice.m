%clear everything
clear all
clear all global
clc
%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\2P\Arc_Pooled_Multi_Day'; %folder to save files to
dataset = 'exp_list_arc_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
ref_str = 'runs-003';
lacz_d1_ori_all = [];
lacz_d2_ori_all = [];
lacz_d3_ori_all = [];
lacz_d1_matches_all = [];
lacz_d2_matches_all = [];
lacz_d3_matches_all = [];
lacz_d1_tc_all = [];
lacz_d2_tc_all = [];
lacz_d3_tc_all = [];
lacz_d1_k_all = [];
lacz_d1_k_max_all = [];
lacz_d1_max_all = [];
lacz_d1_k_red_all = [];
lacz_d1_k_tuned_red_all = [];
lacz_d1_max_red_all = [];
lacz_d1_max_tuned_red_all = [];

enh_d1_ori_all = [];
enh_d2_ori_all = [];
enh_d3_ori_all = [];
enh_d1_matches_all = [];
enh_d2_matches_all = [];
enh_d3_matches_all = [];
enh_d1_tc_all = [];
enh_d2_tc_all = [];
enh_d3_tc_all = [];
enh_d1_k_all = [];
enh_d1_max_all = [];
enh_d1_k_max_all = [];
enh_d1_k_red_all = [];
enh_d1_k_tuned_red_all = [];
enh_d1_max_red_all = [];
enh_d1_max_tuned_red_all = [];

lacz_d1 = [49 52 55 58 61 64];
lacz_d2 = lacz_d1+1;
lacz_d3 = lacz_d2+1;

enh_d1 = [67 70 73 76 79];
enh_d2 = enh_d1+1;
enh_d3 = enh_d2+1;
%%
%lacz d1
for isess = lacz_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    lacz_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    lacz_d1_ori_all = [lacz_d1_ori_all lacz_d1_ori];
    lacz_d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    lacz_d1_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    lacz_d1_k_max_all = [lacz_d1_k_max_all lacz_d1_k_max];
    lacz_d1_k_red = lacz_d1_k_max.k1(lacz_d1_matches.redCells);
    lacz_d1_k_red_all = [lacz_d1_k_red_all lacz_d1_k_red];
    lacz_d1_max_red = lacz_d1_k_max.max_dfof(lacz_d1_matches.redCells);
    lacz_d1_max_red_all = [lacz_d1_max_red_all lacz_d1_max_red];
    lacz_d1_k_red_tuned_index = intersect(lacz_d1_matches.redCells, lacz_d1_ori.ind_theta90);
    lacz_d1_k_tuned_red = lacz_d1_k_max.k1(lacz_d1_k_red_tuned_index);
    lacz_d1_k_tuned_red_all = [lacz_d1_k_tuned_red_all lacz_d1_k_tuned_red];
    lacz_d1_max_tuned_red = lacz_d1_k_max.max_dfof(lacz_d1_k_red_tuned_index);
    lacz_d1_max_tuned_red_all = [lacz_d1_max_tuned_red_all lacz_d1_max_tuned_red];
    lacz_d1_matches_all = [lacz_d1_matches_all lacz_d1_matches];

end

%lacz d2
for isess = lacz_d2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    lacz_d2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    lacz_d2_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    lacz_d2_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    lacz_d2_ori_all = [lacz_d2_ori_all lacz_d2_ori];
    lacz_d2_matches_all = [lacz_d2_matches_all lacz_d2_matches];
    lacz_d2_tc_all = [lacz_d2_tc_all lacz_d2_tc];
end

%lacz d3
for isess = lacz_d3
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    lacz_d3_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    lacz_d3_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    lacz_d3_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    lacz_d3_ori_all = [lacz_d3_ori_all lacz_d3_ori];
    lacz_d3_matches_all = [lacz_d3_matches_all lacz_d3_matches];
    lacz_d3_tc_all = [lacz_d3_tc_all lacz_d3_tc];
end

%%
lacz_match_d2_all = [];
lacz_match_d3_all = [];
lacz_prefori_d1_d2_match_all = [];
lacz_prefori_d2_match_all = [];
lacz_prefori_d1_d3_match_all = [];
lacz_prefori_d3_match_all = [];
lacz_prefori_d1_d2_match_tune_all = [];
lacz_prefori_d2_match_tune_all = [];
lacz_prefori_d1_d3_match_tune_all = [];
lacz_prefori_d3_match_tune_all = [];
lacz_tuned_d1_all = [];
lacz_tuned_d2_all = [];
lacz_tuned_d3_all = [];
lacz_tuned_matched_d2_all = [];
lacz_tuned_matched_d3_all = [];

lacz_d1_k_tuned_all = [];
lacz_d1_max_tuned_all = [];
lacz_d1_k_any_all = [];
lacz_d1_max_any_all = [];

for idata = 1:length(lacz_d1_ori_all)
    lacz_match_d2 = lacz_d2_tc_all(idata).match_ind;
    lacz_match_d3 = lacz_d3_tc_all(idata).match_ind;
    lacz_prefori_d1_d2_match = lacz_d1_ori_all(idata).prefOri(1,lacz_match_d2);
    lacz_prefori_d2_match = lacz_d2_ori_all(idata).prefOri(1,lacz_match_d2);
    lacz_prefori_d1_d3_match = lacz_d1_ori_all(idata).prefOri(1,lacz_match_d3);
    lacz_prefori_d3_match = lacz_d3_ori_all(idata).prefOri(1,lacz_match_d3);
    lacz_match_d2_all = [lacz_match_d2_all lacz_match_d2];
    lacz_match_d3_all = [lacz_match_d3_all lacz_match_d3];
    lacz_prefori_d1_d2_match_all = [lacz_prefori_d1_d2_match_all lacz_prefori_d1_d2_match];
    lacz_prefori_d2_match_all = [lacz_prefori_d2_match_all lacz_prefori_d2_match];
    lacz_prefori_d1_d3_match_all = [lacz_prefori_d1_d3_match_all lacz_prefori_d1_d3_match];
    lacz_prefori_d3_match_all = [lacz_prefori_d3_match_all lacz_prefori_d3_match];
    
    lacz_tuned_d1 = lacz_d1_ori_all(idata).ind_theta90;
    lacz_tuned_d2 = lacz_d2_ori_all(idata).ind_theta90;
    lacz_tuned_d3 = lacz_d3_ori_all(idata).ind_theta90;
    lacz_tuned_matched_d2 = find(ismember(lacz_match_d2, lacz_tuned_d1) | ismember(lacz_match_d2, lacz_tuned_d2));
    lacz_tuned_matched_d3 = find(ismember(lacz_match_d3, lacz_tuned_d1) | ismember(lacz_match_d3, lacz_tuned_d3));
    lacz_tuned_d1_all = [lacz_tuned_d1_all, lacz_tuned_d1];
    lacz_tuned_d2_all = [lacz_tuned_d2_all, lacz_tuned_d2];
    lacz_tuned_d3_all = [lacz_tuned_d3_all, lacz_tuned_d3];
    lacz_tuned_matched_d2_all = [lacz_tuned_matched_d2_all, lacz_tuned_matched_d2];
    lacz_tuned_matched_d3_all = [lacz_tuned_matched_d3_all, lacz_tuned_matched_d3];

    lacz_prefori_d1_d2_match_tune = lacz_d1_ori_all(idata).prefOri(1, lacz_tuned_matched_d2);
    lacz_prefori_d1_d2_match_tune_all = [lacz_prefori_d1_d2_match_tune_all lacz_prefori_d1_d2_match_tune];
    lacz_prefori_d2_match_tune = lacz_d2_ori_all(idata).prefOri(1, lacz_tuned_matched_d2);
    lacz_prefori_d2_match_tune_all = [lacz_prefori_d2_match_tune_all lacz_prefori_d2_match_tune];    

    lacz_prefori_d1_d3_match_tune = lacz_d1_ori_all(idata).prefOri(1, lacz_tuned_matched_d3);
    lacz_prefori_d1_d3_match_tune_all = [lacz_prefori_d1_d3_match_tune_all lacz_prefori_d1_d3_match_tune];
    lacz_prefori_d3_match_tune = lacz_d3_ori_all(idata).prefOri(1, lacz_tuned_matched_d3);
    lacz_prefori_d3_match_tune_all = [lacz_prefori_d3_match_tune_all lacz_prefori_d3_match_tune];
    
    lacz_d1_k_tuned = lacz_d1_k_max_all(idata).k1(lacz_tuned_d1);
    lacz_d1_k_tuned_all = [lacz_d1_k_tuned_all lacz_d1_k_tuned];
    lacz_d1_max_tuned = lacz_d1_k_max_all(idata).max_dfof(lacz_tuned_d1);
    lacz_d1_max_tuned_all = [lacz_d1_max_tuned_all lacz_d1_max_tuned];

    lacz_d1_k_any = lacz_d1_k_max_all(idata).k1;
    lacz_d1_k_any_all = [lacz_d1_k_any_all lacz_d1_k_any];
    lacz_d1_max_any = lacz_d1_k_max_all(idata).max_dfof;
    lacz_d1_max_any_all = [lacz_d1_max_any_all lacz_d1_max_any];
end

%%
figure('Position', [400 20 650 700]);
sgtitle(['All LacZ Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(lacz_prefori_d1_d2_match_all, lacz_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
%hold on
%scatter(enh_prefori_d1_d2_match_all(enh_tuned_matched_d2_all), enh_prefori_d2_match_all(enh_tuned_matched_d2_all),'r');
%hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Not Tuned (n = ', num2str(length(enh_match_d2_all)-length(enh_tuned_matched_d2_all)), ')'], ['Tuned (n = ', num2str(length(enh_tuned_matched_d2_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(lacz_prefori_d1_d3_match_all, lacz_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
%hold on
%scatter(enh_prefori_d1_d3_match_all(enh_tuned_matched_d3_all), enh_prefori_d3_match_all(enh_tuned_matched_d3_all),'r');
%hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
print(fullfile(newfnout, ['Grace LacZ mice', '_ori_scatter.pdf']), '-dpdf', '-bestfit')
%legend(['Not Tuned (n = ', num2str(length(enh_match_d3_all)-length(enh_tuned_matched_d3_all)), ')'], ['Tuned (n = ', num2str(length(enh_tuned_matched_d3_all)), ')'], 'Location', 'northwest')

%%
%enh d1
for isess = enh_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    enh_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    enh_d1_ori_all = [enh_d1_ori_all enh_d1_ori];
    enh_d1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    enh_d1_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    enh_d1_k_max_all = [enh_d1_k_max_all enh_d1_k_max];
    enh_d1_k_red = enh_d1_k_max.k1(enh_d1_matches.redCells);
    enh_d1_k_red_all = [enh_d1_k_red_all enh_d1_k_red];
    enh_d1_max_red = enh_d1_k_max.max_dfof(enh_d1_matches.redCells);
    enh_d1_max_red_all = [enh_d1_max_red_all enh_d1_max_red];
    enh_d1_k_red_tuned_index = intersect(enh_d1_matches.redCells, enh_d1_ori.ind_theta90);
    enh_d1_k_tuned_red = enh_d1_k_max.k1(enh_d1_k_red_tuned_index);
    enh_d1_k_tuned_red_all = [enh_d1_k_tuned_red_all enh_d1_k_tuned_red];
    enh_d1_max_tuned_red = enh_d1_k_max.max_dfof(enh_d1_k_red_tuned_index);
    enh_d1_max_tuned_red_all = [enh_d1_max_tuned_red_all enh_d1_max_tuned_red];
    enh_d1_matches_all = [enh_d1_matches_all enh_d1_matches];
end

%enh d2
for isess = enh_d2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    enh_d2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    enh_d2_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    enh_d2_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    enh_d2_ori_all = [enh_d2_ori_all enh_d2_ori];
    enh_d2_matches_all = [enh_d2_matches_all enh_d2_matches];
    enh_d2_tc_all = [enh_d2_tc_all enh_d2_tc];
end

%enh d3
for isess = enh_d3
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    enh_d3_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    enh_d3_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    enh_d3_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    enh_d3_ori_all = [enh_d3_ori_all enh_d3_ori];
    enh_d3_matches_all = [enh_d3_matches_all enh_d3_matches];
    enh_d3_tc_all = [enh_d3_tc_all enh_d3_tc];
end
%%
enh_match_d2_all = [];
enh_match_d3_all = [];
enh_prefori_d1_d2_match_all = [];
enh_prefori_d2_match_all = [];
enh_prefori_d1_d3_match_all = [];
enh_prefori_d3_match_all = [];
enh_prefori_d1_d2_match_tune_all = [];
enh_prefori_d2_match_tune_all = [];
enh_prefori_d1_d3_match_tune_all = [];
enh_prefori_d3_match_tune_all = [];
enh_tuned_d1_all = [];
enh_tuned_d2_all = [];
enh_tuned_d3_all = [];
enh_tuned_matched_d2_all = [];
enh_tuned_matched_d3_all = [];

enh_d1_k_tuned_all = [];
enh_d1_max_tuned_all = [];
enh_d1_k_any_all = [];
enh_d1_max_any_all = [];

for idata = 1:length(enh_d1_ori_all)
    enh_match_d2 = enh_d2_tc_all(idata).match_ind;
    enh_match_d3 = enh_d3_tc_all(idata).match_ind;
    enh_prefori_d1_d2_match = enh_d1_ori_all(idata).prefOri(1,enh_match_d2);
    enh_prefori_d2_match = enh_d2_ori_all(idata).prefOri(1,enh_match_d2);
    enh_prefori_d1_d3_match = enh_d1_ori_all(idata).prefOri(1,enh_match_d3);
    enh_prefori_d3_match = enh_d3_ori_all(idata).prefOri(1,enh_match_d3);
    enh_match_d2_all = [enh_match_d2_all enh_match_d2];
    enh_match_d3_all = [enh_match_d3_all enh_match_d3];
    enh_prefori_d1_d2_match_all = [enh_prefori_d1_d2_match_all enh_prefori_d1_d2_match];
    enh_prefori_d2_match_all = [enh_prefori_d2_match_all enh_prefori_d2_match];
    enh_prefori_d1_d3_match_all = [enh_prefori_d1_d3_match_all enh_prefori_d1_d3_match];
    enh_prefori_d3_match_all = [enh_prefori_d3_match_all enh_prefori_d3_match];
    
    enh_tuned_d1 = enh_d1_ori_all(idata).ind_theta90;
    enh_tuned_d2 = enh_d2_ori_all(idata).ind_theta90;
    enh_tuned_d3 = enh_d3_ori_all(idata).ind_theta90;
    enh_tuned_matched_d2 = find(ismember(enh_match_d2, enh_tuned_d1) | ismember(enh_match_d2, enh_tuned_d2));
    enh_tuned_matched_d3 = find(ismember(enh_match_d3, enh_tuned_d1) | ismember(enh_match_d3, enh_tuned_d3));
    enh_tuned_d1_all = [enh_tuned_d1_all, enh_tuned_d1];
    enh_tuned_d2_all = [enh_tuned_d2_all, enh_tuned_d2];
    enh_tuned_d3_all = [enh_tuned_d3_all, enh_tuned_d3];
    enh_tuned_matched_d2_all = [enh_tuned_matched_d2_all, enh_tuned_matched_d2];
    enh_tuned_matched_d3_all = [enh_tuned_matched_d3_all, enh_tuned_matched_d3];

    enh_prefori_d1_d2_match_tune = enh_d1_ori_all(idata).prefOri(1, enh_tuned_matched_d2);
    enh_prefori_d1_d2_match_tune_all = [enh_prefori_d1_d2_match_tune_all enh_prefori_d1_d2_match_tune];
    enh_prefori_d2_match_tune = enh_d2_ori_all(idata).prefOri(1, enh_tuned_matched_d2);
    enh_prefori_d2_match_tune_all = [enh_prefori_d2_match_tune_all enh_prefori_d2_match_tune];    

    enh_prefori_d1_d3_match_tune = enh_d1_ori_all(idata).prefOri(1, enh_tuned_matched_d3);
    enh_prefori_d1_d3_match_tune_all = [enh_prefori_d1_d3_match_tune_all enh_prefori_d1_d3_match_tune];
    enh_prefori_d3_match_tune = enh_d3_ori_all(idata).prefOri(1, enh_tuned_matched_d3);
    enh_prefori_d3_match_tune_all = [enh_prefori_d3_match_tune_all enh_prefori_d3_match_tune];    

    enh_d1_k_tuned = enh_d1_k_max_all(idata).k1(enh_tuned_d1);
    enh_d1_k_tuned_all = [enh_d1_k_tuned_all enh_d1_k_tuned];
    enh_d1_max_tuned = enh_d1_k_max_all(idata).max_dfof(enh_tuned_d1);
    enh_d1_max_tuned_all = [enh_d1_max_tuned_all enh_d1_max_tuned];
end
%%
%plot - something is wrong with the n sizes*** - it is because the index
%numbers repeat for the scatter
figure('Position', [400 20 650 700]);
sgtitle(['All enh-Enh Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(enh_prefori_d1_d2_match_all, enh_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
%hold on
%scatter(enh_prefori_d1_d2_match_all(enh_tuned_matched_d2_all), enh_prefori_d2_match_all(enh_tuned_matched_d2_all),'r');
%hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
%legend(['Not Tuned (n = ', num2str(length(enh_match_d2_all)-length(enh_tuned_matched_d2_all)), ')'], ['Tuned (n = ', num2str(length(enh_tuned_matched_d2_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(enh_prefori_d1_d3_match_all, enh_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
%hold on
%scatter(enh_prefori_d1_d3_match_all(enh_tuned_matched_d3_all), enh_prefori_d3_match_all(enh_tuned_matched_d3_all),'r');
%hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
print(fullfile(newfnout, ['Grace enh mice', '_ori_scatter.pdf']), '-dpdf', '-bestfit')
%%
%plots of d1 not matched cells (red only!)
figure;
cdfplot(enh_d1_k_red_all);
hold on
cdfplot(lacz_d1_k_red_all);
legend(['enh-Enh (n = ', num2str(length(enh_d1_k_red_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_k_red_all)), ')']', 'Location', 'Best')
title('Day 1 K Value Distribution')
hold off
print(fullfile(newfnout, ['Grace all mice red cells', '_k_dist.pdf']), '-dpdf', '-bestfit')

figure;
cdfplot(enh_d1_max_red_all);
hold on
cdfplot(lacz_d1_max_red_all);
legend(['enh-Enh (n = ', num2str(length(enh_d1_k_red_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_k_red_all)), ')']', 'Location', 'Best')
title('Day 1 Max dF/F Distribution')
hold off
print(fullfile(newfnout, ['Grace all mice red cells', 'max_dist.pdf']), '-dpdf', '-bestfit')
%%
%days 1 and 2 ori scatter same graph
%enh
figure;
sgtitle(['Grace Enhancer Mice'], 'Interpreter', 'None');
scatter(enh_prefori_d1_d2_match_all, enh_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8);
hold on
scatter(enh_prefori_d1_d3_match_all, enh_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day x Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend('Day 2', 'Day 3', 'Location', 'best')
print(fullfile(newfnout, ['Grace enh mice', 'ori_scatter_both_days.pdf']), '-dpdf', '-bestfit')

%lacz
figure;
sgtitle(['Grace LacZ Mice'], 'Interpreter', 'None');
scatter(lacz_prefori_d1_d2_match_all, lacz_prefori_d2_match_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8);
hold on
scatter(lacz_prefori_d1_d3_match_all, lacz_prefori_d3_match_all, 'filled', 'MarkerEdgeColor', 'y', 'MarkerFaceColor', 'k', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day x Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend('Day 2', 'Day 3', 'Location', 'best')
print(fullfile(newfnout, ['Grace lacz mice', 'ori_scatter_both_days.pdf']), '-dpdf', '-bestfit')

%%
%enh-enh and lacz scatter pref ori day 1 same plot
figure;
%sgtitle(['All Mice Day 1 vs. Day 2 Pref Ori'], 'Interpreter', 'None');
scatter(enh_prefori_d1_d2_match_tune_all, enh_prefori_d2_match_tune_all, 'filled', 'MarkerEdgeColor', '#00841a', 'MarkerFaceColor', '#00841a', 'LineWidth', 0.8);
hold on
scatter(lacz_prefori_d1_d2_match_tune_all, lacz_prefori_d2_match_tune_all, 'MarkerEdgeColor', 'r', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.8);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend('Enhancer', 'LacZ', 'Location', 'best')
print(fullfile(newfnout, ['Grace mice tuned matched red', 'ori_scatter_d1_d2.pdf']), '-dpdf', '-bestfit')

figure;
%sgtitle(['All Mice Day 1 vs. Day 3 Pref Ori'], 'Interpreter', 'None');
scatter(enh_prefori_d1_d3_match_tune_all, enh_prefori_d3_match_tune_all, 'filled', 'MarkerEdgeColor', '#00841a', 'MarkerFaceColor', '#00841a', 'LineWidth', 0.8);
hold on
scatter(lacz_prefori_d1_d3_match_tune_all, lacz_prefori_d3_match_tune_all, 'MarkerEdgeColor', 'r', 'LineWidth', 0.8, 'MarkerFaceAlpha', 0.8);
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend('Enhancer', 'LacZ', 'Location', 'best')
print(fullfile(newfnout, ['Grace mice tuned matched red', 'ori_scatter_d1_d3.pdf']), '-dpdf', '-bestfit')

%%
%only tuned k and max dist
figure;
h = cdfplot(enh_d1_k_tuned_red_all);
hold on
j = cdfplot(lacz_d1_k_tuned_red_all);
%cdfplot(enh_d1_k_all);
%legend('enh', 'LacZ', 'Enhancer', 'Location', 'Best')
%legend(['enh (n = ', num2str(length(enh_d1_k_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_k_all)), ')']', 'Location', 'Best')
%title('Day 1 K Value Distribution')
title('')
set(h, 'Color', '#00841a', 'LineWidth', 1.0)
set(j, 'Color', 'r', 'LineWidth', 1.0)
legend(['Enhancer (n = ', num2str(length(enh_d1_k_tuned_red_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_k_tuned_red_all)), ')']', 'Location', 'Best')
xlabel('k Value')
print(fullfile(newfnout, ['Grace mice red and tuned cells', '_k_dist.pdf']), '-dpdf', '-bestfit')
hold off

figure;
m = cdfplot(enh_d1_max_tuned_red_all);
hold on
n = cdfplot(lacz_d1_max_tuned_red_all);
%cdfplot(enh_d1_k_all);
%legend('enh', 'LacZ', 'Enhancer', 'Location', 'Best')
%legend(['enh (n = ', num2str(length(enh_d1_k_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_k_all)), ')']', 'Location', 'Best')
%title('Day 1 Max dF/F Value Distribution')
title('')
set(m, 'Color', '#00841a', 'LineWidth', 1.0)
set(n, 'Color', 'r', 'LineWidth', 1.0)
legend(['Enhancer (n = ', num2str(length(enh_d1_max_tuned_red_all)), ')'], ['LacZ (n = ', num2str(length(lacz_d1_max_tuned_red_all)), ')']', 'Location', 'Best')
xlabel('Max dF/F')
hold off
print(fullfile(newfnout, ['Grace mice red and tuned cells', 'max_dist.pdf']), '-dpdf', '-bestfit')
hold off

%%
%individual mice d1 v d2 enh and lacz
mouse1_enh_d1_k_red_tuned_index_ses1 = intersect(enh_d1_matches_all(1).redCells, enh_d1_ori_all(1).ind_theta90);
mouse1_enh_k_d1_ses1 = enh_d1_k_max_all(1).k1(mouse1_enh_d1_k_red_tuned_index_ses1);
mouse1_enh_max_d1_ses1 = enh_d1_k_max_all(1).max_dfof(mouse1_enh_d1_k_red_tuned_index_ses1);
mouse1_enh_d1_k_red_tuned_index_ses2 = intersect(enh_d1_matches_all(2).redCells, enh_d1_ori_all(2).ind_theta90);
mouse1_enh_k_d1_ses2 = enh_d1_k_max_all(2).k1(mouse1_enh_d1_k_red_tuned_index_ses2);
mouse1_enh_max_d1_ses2 = enh_d1_k_max_all(2).max_dfof(mouse1_enh_d1_k_red_tuned_index_ses2);
mouse1_enh_k_d1_all = [mouse1_enh_k_d1_ses1 mouse1_enh_k_d1_ses2];
mouse1_enh_max_d1_all = [mouse1_enh_max_d1_ses1 mouse1_enh_max_d1_ses2];

mouse2_enh_d1_k_red_tuned_index_ses1 = intersect(enh_d1_matches_all(3).redCells, enh_d1_ori_all(3).ind_theta90);
mouse2_enh_k_d1_ses1 = enh_d1_k_max_all(3).k1(mouse2_enh_d1_k_red_tuned_index_ses1);
mouse2_enh_max_d1_ses1 = enh_d1_k_max_all(3).max_dfof(mouse2_enh_d1_k_red_tuned_index_ses1);
%mouse2_enh_d1_k_red_tuned_index_ses2 = intersect(enh_d1_matches_all(4).redCells, enh_d1_ori_all(4).ind_theta90);
%mouse2_enh_k_d1_ses2 = enh_d1_k_max_all(4).k1(mouse2_enh_d1_k_red_tuned_index_ses2);
%mouse2_enh_max_d1_ses2 = enh_d1_k_max_all(4).max_dfof(mouse2_enh_d1_k_red_tuned_index_ses2);
%mouse2_enh_k_d1_all = [mouse2_enh_k_d1_ses1 mouse2_enh_k_d1_ses2];
mouse2_enh_k_d1_all = [mouse2_enh_k_d1_ses1];
%mouse2_enh_max_d1_all = [mouse2_enh_max_d1_ses1 mouse2_enh_max_d1_ses2];
mouse2_enh_max_d1_all = [mouse2_enh_max_d1_ses1];

mouse3_enh_d1_k_red_tuned_index_ses1 = intersect(enh_d1_matches_all(4).redCells, enh_d1_ori_all(4).ind_theta90);
mouse3_enh_k_d1_ses1 = enh_d1_k_max_all(4).k1(mouse3_enh_d1_k_red_tuned_index_ses1);
mouse3_enh_max_d1_ses1 = enh_d1_k_max_all(4).max_dfof(mouse3_enh_d1_k_red_tuned_index_ses1);
mouse3_enh_d1_k_red_tuned_index_ses2 = intersect(enh_d1_matches_all(5).redCells, enh_d1_ori_all(5).ind_theta90);
mouse3_enh_max_d1_ses2 = enh_d1_k_max_all(5).max_dfof(mouse3_enh_d1_k_red_tuned_index_ses2);
mouse3_enh_k_d1_ses2 = enh_d1_k_max_all(5).k1(mouse3_enh_d1_k_red_tuned_index_ses2);
mouse3_enh_k_d1_all = [mouse3_enh_k_d1_ses1 mouse3_enh_k_d1_ses2];
mouse3_enh_max_d1_all = [mouse3_enh_max_d1_ses1 mouse3_enh_max_d1_ses2];

mouse1_lacz_d1_k_red_tuned_index_ses1 = intersect(lacz_d1_matches_all(1).redCells, lacz_d1_ori_all(1).ind_theta90);
mouse1_lacz_k_d1_ses1 = lacz_d1_k_max_all(1).k1(mouse1_lacz_d1_k_red_tuned_index_ses1);
mouse1_lacz_max_d1_ses1 = lacz_d1_k_max_all(1).max_dfof(mouse1_lacz_d1_k_red_tuned_index_ses1);
mouse1_lacz_d1_k_red_tuned_index_ses2 = intersect(lacz_d1_matches_all(2).redCells, lacz_d1_ori_all(2).ind_theta90);
mouse1_lacz_k_d1_ses2 = lacz_d1_k_max_all(2).k1(mouse1_lacz_d1_k_red_tuned_index_ses2);
mouse1_lacz_max_d1_ses2 = lacz_d1_k_max_all(2).max_dfof(mouse1_lacz_d1_k_red_tuned_index_ses2);
mouse1_lacz_k_d1_all = [mouse1_lacz_k_d1_ses1 mouse1_lacz_k_d1_ses2];
mouse1_lacz_max_d1_all = [mouse1_lacz_max_d1_ses1 mouse1_lacz_max_d1_ses2];

mouse2_lacz_d1_k_red_tuned_index_ses1 = intersect(lacz_d1_matches_all(3).redCells, lacz_d1_ori_all(3).ind_theta90);
mouse2_lacz_k_d1_ses1 = lacz_d1_k_max_all(3).k1(mouse2_lacz_d1_k_red_tuned_index_ses1);
mouse2_lacz_max_d1_ses1 = lacz_d1_k_max_all(3).max_dfof(mouse2_lacz_d1_k_red_tuned_index_ses1);
mouse2_lacz_d1_k_red_tuned_index_ses2 = intersect(lacz_d1_matches_all(4).redCells, lacz_d1_ori_all(4).ind_theta90);
mouse2_lacz_k_d1_ses2 = lacz_d1_k_max_all(4).k1(mouse2_lacz_d1_k_red_tuned_index_ses2);
mouse2_lacz_max_d1_ses2 = lacz_d1_k_max_all(4).max_dfof(mouse2_lacz_d1_k_red_tuned_index_ses2);
mouse2_lacz_k_d1_all = [mouse2_lacz_k_d1_ses1 mouse2_lacz_k_d1_ses2];
mouse2_lacz_max_d1_all = [mouse2_lacz_max_d1_ses1 mouse2_lacz_max_d1_ses2];

mouse3_lacz_d1_k_red_tuned_index_ses1 = intersect(lacz_d1_matches_all(5).redCells, lacz_d1_ori_all(5).ind_theta90);
mouse3_lacz_k_d1_ses1 = lacz_d1_k_max_all(5).k1(mouse3_lacz_d1_k_red_tuned_index_ses1);
mouse3_lacz_max_d1_ses1 = lacz_d1_k_max_all(5).max_dfof(mouse3_lacz_d1_k_red_tuned_index_ses1);
mouse3_lacz_d1_k_red_tuned_index_ses2 = intersect(lacz_d1_matches_all(6).redCells, lacz_d1_ori_all(6).ind_theta90);
mouse3_lacz_k_d1_ses2 = lacz_d1_k_max_all(6).k1(mouse3_lacz_d1_k_red_tuned_index_ses2);
mouse3_lacz_max_d1_ses2 = lacz_d1_k_max_all(6).max_dfof(mouse3_lacz_d1_k_red_tuned_index_ses2);
mouse3_lacz_k_d1_all = [mouse3_lacz_k_d1_ses1 mouse3_lacz_k_d1_ses2];
mouse3_lacz_max_d1_all = [mouse3_lacz_max_d1_ses1 mouse3_lacz_max_d1_ses2];

figure; 
m1_enh_k_d1 = cdfplot(mouse1_enh_k_d1_all);
hold on
m2_enh_k_d1 = cdfplot(mouse2_enh_k_d1_all);
m3_enh_k_d1 = cdfplot(mouse3_enh_k_d1_all);
set(m1_enh_k_d1, 'Color', '#00841a', 'LineWidth', 1.0);
set(m2_enh_k_d1, 'Color', '#00841a', 'LineWidth', 1.0);
set(m3_enh_k_d1, 'Color', '#00841a', 'LineWidth', 1.0);

m1_lacz_k_d1 = cdfplot(mouse1_lacz_k_d1_all);
m2_lacz_k_d1 = cdfplot(mouse2_lacz_k_d1_all);
m3_lacz_k_d1 = cdfplot(mouse3_lacz_k_d1_all);
set(m1_lacz_k_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m2_lacz_k_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m3_lacz_k_d1, 'Color', 'r', 'LineWidth', 1.0);

title('')
xlabel('k Value')
xlim([0 30])
legend([m1_enh_k_d1(1), m1_lacz_k_d1(1)], 'Arc Enhancer', ' LacZ')
hold off

print(fullfile(newfnout, ['individual_Grace_mice', '_d1_k_cdf.pdf']), '-dpdf', '-bestfit')

figure; 
m1_enh_max_d1 = cdfplot(mouse1_enh_max_d1_all);
hold on
m2_enh_max_d1 = cdfplot(mouse2_enh_max_d1_all);
m3_enh_max_d1 = cdfplot(mouse3_enh_max_d1_all);
set(m1_enh_max_d1, 'Color', '#00841a', 'LineWidth', 1.0);
set(m2_enh_max_d1, 'Color', '#00841a', 'LineWidth', 1.0);
set(m3_enh_max_d1, 'Color', '#00841a', 'LineWidth', 1.0);

m1_lacz_max_d1 = cdfplot(mouse1_lacz_max_d1_all);
m2_lacz_max_d1 = cdfplot(mouse2_lacz_max_d1_all);
m3_lacz_max_d1 = cdfplot(mouse3_lacz_max_d1_all);
set(m1_lacz_max_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m2_lacz_max_d1, 'Color', 'r', 'LineWidth', 1.0);
set(m3_lacz_max_d1, 'Color', 'r', 'LineWidth', 1.0);

title('')
xlabel('Max dF/F Value')
%xlim([0 30])
legend([m1_enh_max_d1(1), m1_lacz_max_d1(1)], 'Arc Enhancer', ' LacZ')
hold off

print(fullfile(newfnout, ['individual_Grace_mice', '_d1_max_cdf.pdf']), '-dpdf', '-bestfit')

%now do this for enhancer vs. lacz!