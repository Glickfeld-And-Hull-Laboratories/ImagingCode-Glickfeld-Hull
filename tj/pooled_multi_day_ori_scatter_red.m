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
arc_list_mouse = [];
arc_list_img_area = [];
arc_list_img_layer = [];
arc_list_img_day = [];
arc_list_date = [];
lacz_list_mouse = [];
lacz_list_img_area = [];
lacz_list_img_layer = [];
lacz_list_img_day = [];
lacz_list_date = [];
arc_data = [];
lacz_data = [];
ref_str = 'runs-003';

%%trying to make a loop that knows the groups and such

for isess = 1:length(expt)
    
    if strcmp(expt(isess).red_indicator{1}, 'Arc-gRNA')
        arc_mouse = expt(isess).mouse;
        arc_img_area = expt(isess).img_loc{1};
        arc_img_layer = expt(isess).img_loc{2};
        arc_img_day = expt(isess).img_day;
            if strcmp(arc_img_day, '4')
                arc_img_day = '1';
            elseif strcmp(arc_img_day, '5')
                arc_img_day = '2';
            elseif strcmp(arc_img_day, '6')
                arc_img_day = '3';
            end
        arc_date = expt(isess).date;
        arc_list_mouse = [arc_list_mouse; arc_mouse];
        arc_list_img_area = [arc_list_img_area; arc_img_area];
        arc_list_img_layer = [arc_list_img_layer; arc_img_layer];
        arc_list_img_day = [arc_list_img_day; arc_img_day];
        arc_list_date = [arc_list_date; arc_date];
        
    elseif strcmp(expt(isess).red_indicator{1}, 'LacZ-gRNA')
        lacz_mouse = expt(isess).mouse;
        lacz_img_area = expt(isess).img_loc{1};
        lacz_img_layer = expt(isess).img_loc{2};
        lacz_img_day = expt(isess).img_day;
            if strcmp(lacz_img_day, '4')
                lacz_img_day = '1';
            elseif strcmp(lacz_img_day, '5')
                lacz_img_day = '2';
            elseif strcmp(lacz_img_day, '6')
                lacz_img_day = '3';
            end
        lacz_date = expt(isess).date;
        lacz_list_mouse = [lacz_list_mouse; lacz_mouse];
        lacz_list_img_area = [lacz_list_img_area; lacz_img_area];
        lacz_list_img_layer = [lacz_list_img_layer; lacz_img_layer];
        lacz_list_img_day = [lacz_list_img_day; lacz_img_day];
        lacz_list_date = [lacz_list_date; lacz_date];
        
    end
    
end    


%%
ref_str = 'runs-003';
arc_d1_ori_all = [];
arc_d2_ori_all = [];
arc_d3_ori_all = [];
arc_d1_matches_all = [];
arc_d2_matches_all = [];
arc_d3_matches_all = [];
arc_d1_tc_all = [];
arc_d2_tc_all = [];
arc_d3_tc_all = [];

arc_d1 = [7 10 19 22];
arc_d2 = arc_d1+1;
arc_d3 = arc_d2+1;

lacz_d1 = [1 4 13 16];
lacz_d2 = lacz_d1+1;
lacz_d3 = lacz_d2+1;

%%

%arc d1
for isess = arc_d1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    arc_d1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    arc_d1_ori_all = [arc_d1_ori_all arc_d1_ori];
end

%arc d2
for isess = arc_d2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    arc_d2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    arc_d2_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    arc_d2_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    arc_d2_ori_all = [arc_d2_ori_all arc_d2_ori];
    arc_d2_matches_all = [arc_d2_matches_all arc_d2_matches];
    arc_d2_tc_all = [arc_d2_tc_all arc_d2_tc];
end

%arc d3
for isess = arc_d3
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    arc_d3_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    arc_d3_matches = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'multiday_alignment.mat']));
    arc_d3_tc = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'TCs.mat']));
    arc_d3_ori_all = [arc_d3_ori_all arc_d3_ori];
    arc_d3_matches_all = [arc_d3_matches_all arc_d3_matches];
    arc_d3_tc_all = [arc_d3_tc_all arc_d3_tc];
end

%%
arc_match_d2_all = [];
arc_match_d3_all = [];
arc_prefori_d1_d2_match_all = [];
arc_prefori_d2_match_all = [];
arc_prefori_d1_d3_match_all = [];
arc_prefori_d3_match_all = [];
arc_tuned_d1_all = [];
arc_tuned_d2_all = [];
arc_tuned_d3_all = [];
arc_tuned_matched_d2_all = [];
arc_tuned_matched_d3_all = [];

for idata = 1:length(arc_d1_ori_all)
    arc_match_d2 = arc_d2_tc_all(idata).match_ind;
    arc_match_d3 = arc_d3_tc_all(idata).match_ind;
    arc_prefori_d1_d2_match = arc_d1_ori_all(idata).prefOri(1,arc_match_d2);
    arc_prefori_d2_match = arc_d2_ori_all(idata).prefOri(1,arc_match_d2);
    arc_prefori_d1_d3_match = arc_d1_ori_all(idata).prefOri(1,arc_match_d3);
    arc_prefori_d3_match = arc_d3_ori_all(idata).prefOri(1,arc_match_d3);
    arc_match_d2_all = [arc_match_d2_all arc_match_d2];
    arc_match_d3_all = [arc_match_d3_all arc_match_d3];
    arc_prefori_d1_d2_match_all = [arc_prefori_d1_d2_match_all arc_prefori_d1_d2_match];
    arc_prefori_d2_match_all = [arc_prefori_d2_match_all arc_prefori_d2_match];
    arc_prefori_d1_d3_match_all = [arc_prefori_d1_d3_match_all arc_prefori_d1_d3_match];
    arc_prefori_d3_match_all = [arc_prefori_d3_match_all arc_prefori_d3_match];
    
    arc_tuned_d1 = arc_d1_ori_all(idata).ind_theta90;
    arc_tuned_d2 = arc_d2_ori_all(idata).ind_theta90;
    arc_tuned_d3 = arc_d3_ori_all(idata).ind_theta90;
    arc_tuned_matched_d2 = find(ismember(arc_match_d2, arc_tuned_d1) & ismember(arc_match_d2, arc_tuned_d2));
    arc_tuned_matched_d3 = find(ismember(arc_match_d3, arc_tuned_d1) & ismember(arc_match_d3, arc_tuned_d3));
    arc_tuned_d1_all = [arc_tuned_d1_all, arc_tuned_d1];
    arc_tuned_d2_all = [arc_tuned_d2_all, arc_tuned_d2];
    arc_tuned_d3_all = [arc_tuned_d3_all, arc_tuned_d3];
    arc_tuned_matched_d2_all = [arc_tuned_matched_d2_all, arc_tuned_matched_d2];
    arc_tuned_matched_d3_all = [arc_tuned_matched_d3_all, arc_tuned_matched_d3];
end


%%
%plot - something is wrong with the n sizes***
figure('Position', [400 40 650 800]);
sgtitle(['All Arc Mice'], 'Interpreter', 'None');
subplot(2,1,1);
scatter(arc_prefori_d1_d2_match_all, arc_prefori_d2_match_all);
hold on
scatter(arc_prefori_d1_d2_match_all(arc_tuned_matched_d2_all), arc_prefori_d2_match_all(arc_tuned_matched_d2_all),'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 2 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
line = refline(1,0);
line.Color = 'k';
legend(['Not Tuned (n = ', num2str(length(arc_match_d2_all)-length(arc_tuned_matched_d2_all)), ')'], ['Tuned (n = ', num2str(length(arc_tuned_matched_d2_all)), ')'], 'Location', 'northwest')
subplot(2,1,2);
scatter(arc_prefori_d1_d3_match_all, arc_prefori_d3_match_all);
hold on
scatter(arc_prefori_d1_d3_match_all(arc_tuned_matched_d3_all), arc_prefori_d3_match_all(arc_tuned_matched_d3_all),'r');
hold off
xlabel('Day 1 Pref Ori');
ylabel('Day 3 Pref Ori');
xlim([0,180]);
ylim([0,180]);
xticks(0:20:180);
yticks(0:20:180);
refline(1,0);
line = refline(1,0);
line.Color = 'k';
legend(['Not Tuned (n = ', num2str(length(arc_match_d3_all)-length(arc_tuned_matched_d3_all)), ')'], ['Tuned (n = ', num2str(length(arc_tuned_matched_d3_all)), ')'], 'Location', 'northwest')

%now do the same as above but with LacZ***