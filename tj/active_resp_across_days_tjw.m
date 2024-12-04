%good_ind is index of active cells
clear all;
clc;
close all;

%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\new_virusTG\active_across_days'; %folder to save files to
dataset = 'exp_list_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%we start with loading the proper information based on the mouse's info in
%the experiment list

d1 = 40; 
d2 = 41; 
d3 = 42;

mouse = expt(d1).mouse; 
group = expt(d1).indicator{1};
ref_str_d1 = ['runs-',expt(d1).runs];
ref_str_d2 = ['runs-',expt(d2).runs]; 
ref_str_d3 = ['runs-',expt(d3).runs]; 
img_area = expt(d1).img_loc{1};
img_layer = expt(d1).img_loc{2};
date_d1 = expt(d1).date; 
date_d2 = expt(d2).date; 
date_d3 = expt(d3).date; 
img_folder_d1 = expt(d1).runs; 
img_folder_d2 = expt(d2).runs; 
img_folder_d3 = expt(d3).runs; 

%load ori info for each day
d1_ori = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningInfo.mat']));
d2_ori = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningInfo.mat']));
d3_ori = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningInfo.mat']));

%load multiday data for day 2
d2_matches = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'multiday_alignment.mat']));
d3_matches = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'multiday_alignment.mat']));

%load k and max vals
d1_k_max = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'k_and_max_vals.mat']));
d2_k_max = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'k_and_max_vals.mat']));
d3_k_max = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'k_and_max_vals.mat']));

%load resp info for each day
d1_resp = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'trialData.mat']));
d2_resp = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'trialData.mat']));
d3_resp = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'trialData.mat']));

% %load fits info for each day
% d1_fits = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'oriTuningAndFits.mat']));
% d2_fits = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'oriTuningAndFits.mat']));
% d3_fits = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'oriTuningAndFits.mat']));
% 
% %load input info for each day
% d1_input = load(fullfile(fnout, [date_d1 '_' mouse], [date_d1 '_' mouse '_' ref_str_d1], [date_d1 '_' mouse '_' ref_str_d1 '_' 'input.mat']));
% d2_input = load(fullfile(fnout, [date_d2 '_' mouse], [date_d2 '_' mouse '_' ref_str_d2], [date_d2 '_' mouse '_' ref_str_d2 '_' 'input.mat']));
% d3_input = load(fullfile(fnout, [date_d3 '_' mouse], [date_d3 '_' mouse '_' ref_str_d3], [date_d3 '_' mouse '_' ref_str_d3 '_' 'input.mat']));



%% matching indices

active_d1 = d1_resp.good_ind';
active_d2 = d2_resp.good_ind';
active_d3 = d3_resp.good_ind';

tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;
tuned_d3 = d3_ori.ind_theta90;

tuned_all = intersect(intersect(tuned_d1, tuned_d2), tuned_d3);

%find green cells that match from d1 to d2 and d1 to d3 
match_d2 = find([d2_matches.cellImageAlign.pass]); 
match_d3 = find([d3_matches.cellImageAlign.pass]); 

match_all = intersect(match_d2, match_d3);

match_tune_all = intersect(match_all, tuned_all);

%%
active_d1_match = intersect(match_all,active_d1);
active_d2_match = intersect(match_all,active_d2);
active_d3_match = intersect(match_all,active_d3);

active_d1_match_tune = intersect(match_tune_all,active_d1);
active_d2_match_tune = intersect(match_tune_all,active_d2);
active_d3_match_tune = intersect(match_tune_all,active_d3);

%%

% resp_d1_match = intersect(match_all,tuned_d1);
% resp_d2_match = intersect(match_all,tuned_d2);
% resp_d3_match = intersect(match_all,tuned_d3);

%%
%total active neurons that are matched
n_active_d1 = length(active_d1_match);
n_active_d2 = length(active_d2_match);
n_active_d3 = length(active_d3_match);
n_active_all = [n_active_d1 n_active_d2 n_active_d3];

figure;
plot(n_active_all)
ylim([0 100])

%fraction of active neurons that are resp
% n_resp_d1 = length(resp_d1_match);
% n_resp_d2 = length(resp_d2_match);
% n_resp_d3 = length(resp_d3_match);
% n_resp_all = [n_resp_d1; n_resp_d2; n_resp_d3];

%%
n_active_tune_d1 = length(active_d1_match_tune);
n_active_tune_d2 = length(active_d2_match_tune);
n_active_tune_d3 = length(active_d3_match_tune);
n_active_tune_all = [n_active_tune_d1 n_active_tune_d2 n_active_tune_d3];

figure;
plot(n_active_tune_all)
ylim([0 100])


%%
%try just d1-d7
tuned_d1_d3 = intersect(tuned_d1, tuned_d3);
match_d1_d3 = match_d3;

tune_match_d1_d3 = intersect(tuned_d1_d3, match_d1_d3);

active_d1_d3 = intersect(active_d1, tune_match_d1_d3);
active_d3_d1 = intersect(active_d3, tune_match_d1_d3);


%%



%%


%%
%if i need to actually do the anova from the paper


for i = 1:nCells
        if ~sum(no_match==i) 
        subplot(n, n2, i)
            for idir = 1:ndir %for each dir
                if nOn>29
                    ind = find(Dir == Dirs(idir)); %find trials of that dir
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1)); %ttest of each cell at each dir for base vs. resp
                if h_dir(i,idir) %if this cell/dir is sig make red, if not make black and plot
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2); %avg response at each dir for each cell
                hold on
            end
            if sum(h_dir(i,:),2)>0 %if cell has one sig dir
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2); %find max index and value
                max_dir(i,:) = max_ind; %which index/dir was highest?
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end




