%this script is the extended dark/light analysis paradigm, with 7 different imaging timepoints -
%this was replaced by a more efficient 4 timepoint version, so although these data were analyzed I
%will not comment this code

%%
%clear everything
clear all
clear all global
clc
close all


%%
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day'; %folder to load files 
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\update'; %folder to save files 
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%
%rad2deg(acos((1./k) .* log((exp(k)./2) + (exp(-k)./2)))) this is the formula
%to convert k value into full width at half max (FWHM)
%%
%baseline1
base1_ori_all = [];
base1_k_max_all = [];
base1_k_tuned_green_all = [];
base1_max_tuned_green_all = [];
base1_mouse_id_all = [];
base1_avgtuningcurve_all = [];


base1 = [1 1+7 1+14];

for isess = base1
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base1_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base1_ori_all = [base1_ori_all base1_ori];
    base1_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base1_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base1_tuned = base1_ori.ind_theta90; %only taking tuned cells
    base1_k_max.k1 = base1_k_max.k1(base1_tuned);
    base1_k_max.max_dfof = base1_k_max.max_dfof(base1_tuned);
    base1_k_max_all = [base1_k_max_all base1_k_max];
    base1_k_tuned_green = base1_k_max.k1;
    base1_k_tuned_green_all = [base1_k_tuned_green_all base1_k_tuned_green];
    base1_max_tuned_green = base1_k_max.max_dfof;
    base1_max_tuned_green_all = [base1_max_tuned_green_all base1_max_tuned_green];
    base1_ncells = length(base1_k_tuned_green);
    base1_mouse_id = string(repmat(mouse,base1_ncells,1))';
    base1_mouse_id_all = [base1_mouse_id_all base1_mouse_id];
    base1_avgtuningcurve = base1_newavg.newAvg(base1_tuned,:)';
    base1_avgtuningcurve_all = [base1_avgtuningcurve_all base1_avgtuningcurve];
end

base1_avgtuningcurve_all = [base1_avgtuningcurve_all]';

base1_tuning_width_all = rad2deg(acos((1./base1_k_tuned_green_all) .* log((exp(base1_k_tuned_green_all)./2) + (exp(-base1_k_tuned_green_all)./2))));

%%
sz = size(base1_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if base1_avgtuningcurve_all(i,j) > 0.5
            base1_avgtuningcurve_all(i,:) = [];
        elseif base1_avgtuningcurve_all(i,j) < -0.5
            base1_avgtuningcurve_all(i,:) = [];
        end
    end
end 


%%
%baseline1 + 4d
base1_4d_ori_all = [];
base1_4d_k_max_all = [];
base1_4d_k_tuned_green_all = [];
base1_4d_max_tuned_green_all = [];
base1_4d_mouse_id_all = [];
base1_4d_avgtuningcurve_all = [];


base1_4d = [2 2+7 2+14];

for isess = base1_4d
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base1_4d_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base1_4d_ori_all = [base1_4d_ori_all base1_4d_ori];
    base1_4d_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base1_4d_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base1_4d_tuned = base1_4d_ori.ind_theta90; %only taking tuned cells
    base1_4d_k_max.k1 = base1_4d_k_max.k1(base1_4d_tuned);
    base1_4d_k_max.max_dfof = base1_4d_k_max.max_dfof(base1_4d_tuned);
    base1_4d_k_max_all = [base1_4d_k_max_all base1_4d_k_max];
    base1_4d_k_tuned_green = base1_4d_k_max.k1;
    base1_4d_k_tuned_green_all = [base1_4d_k_tuned_green_all base1_4d_k_tuned_green];
    base1_4d_max_tuned_green = base1_4d_k_max.max_dfof;
    base1_4d_max_tuned_green_all = [base1_4d_max_tuned_green_all base1_4d_max_tuned_green];
    base1_4d_ncells = length(base1_4d_k_tuned_green);
    base1_4d_mouse_id = string(repmat(mouse,base1_4d_ncells,1))';
    base1_4d_mouse_id_all = [base1_4d_mouse_id_all base1_4d_mouse_id];
    base1_4d_avgtuningcurve = base1_4d_newavg.newAvg(base1_4d_tuned,:)';
    base1_4d_avgtuningcurve_all = [base1_4d_avgtuningcurve_all base1_4d_avgtuningcurve];
end

base1_4d_avgtuningcurve_all = [base1_4d_avgtuningcurve_all]';

base1_4d_tuning_width_all = rad2deg(acos((1./base1_4d_k_tuned_green_all) .* log((exp(base1_4d_k_tuned_green_all)./2) + (exp(-base1_4d_k_tuned_green_all)./2))));

%%
sz = size(base1_4d_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if base1_4d_avgtuningcurve_all(i,j) > 0.5
            base1_4d_avgtuningcurve_all(i,:) = [];
        elseif base1_4d_avgtuningcurve_all(i,j) < -0.5
            base1_4d_avgtuningcurve_all(i,:) = [];
        end
    end
end    
%%
%baseline2
base2_ori_all = [];
base2_k_max_all = [];
base2_k_tuned_green_all = [];
base2_max_tuned_green_all = [];
base2_mouse_id_all = [];
base2_avgtuningcurve_all = [];


base2 = [3 3+7 3+14];

for isess = base2
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base2_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base2_ori_all = [base2_ori_all base2_ori];
    base2_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base2_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base2_tuned = base2_ori.ind_theta90; %only taking tuned cells
    base2_k_max.k1 = base2_k_max.k1(base2_tuned);
    base2_k_max.max_dfof = base2_k_max.max_dfof(base2_tuned);
    base2_k_max_all = [base2_k_max_all base2_k_max];
    base2_k_tuned_green = base2_k_max.k1;
    base2_k_tuned_green_all = [base2_k_tuned_green_all base2_k_tuned_green];
    base2_max_tuned_green = base2_k_max.max_dfof;
    base2_max_tuned_green_all = [base2_max_tuned_green_all base2_max_tuned_green];
    base2_ncells = length(base2_k_tuned_green);
    base2_mouse_id = string(repmat(mouse,base2_ncells,1))';
    base2_mouse_id_all = [base2_mouse_id_all base2_mouse_id];
    base2_avgtuningcurve = base2_newavg.newAvg(base2_tuned,:)';
    base2_avgtuningcurve_all = [base2_avgtuningcurve_all base2_avgtuningcurve];
end

base2_avgtuningcurve_all = [base2_avgtuningcurve_all]';

base2_tuning_width_all = rad2deg(acos((1./base2_k_tuned_green_all) .* log((exp(base2_k_tuned_green_all)./2) + (exp(-base2_k_tuned_green_all)./2))));

%%
sz = size(base2_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if base2_avgtuningcurve_all(i,j) > 0.5
            base2_avgtuningcurve_all(i,:) = [];
        elseif base2_avgtuningcurve_all(i,j) < -0.5
            base2_avgtuningcurve_all(i,:) = [];
        end
    end
end 


%%
%baseline2 + 4hr
base2_4hr_ori_all = [];
base2_4hr_k_max_all = [];
base2_4hr_k_tuned_green_all = [];
base2_4hr_max_tuned_green_all = [];
base2_4hr_mouse_id_all = [];
base2_4hr_avgtuningcurve_all = [];


base2_4hr = [4 4+7 4+14];

for isess = base2_4hr
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    base2_4hr_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    base2_4hr_ori_all = [base2_4hr_ori_all base2_4hr_ori];
    base2_4hr_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    base2_4hr_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    base2_4hr_tuned = base2_4hr_ori.ind_theta90; %only taking tuned cells
    base2_4hr_k_max.k1 = base2_4hr_k_max.k1(base2_4hr_tuned);
    base2_4hr_k_max.max_dfof = base2_4hr_k_max.max_dfof(base2_4hr_tuned);
    base2_4hr_k_max_all = [base2_4hr_k_max_all base2_4hr_k_max];
    base2_4hr_k_tuned_green = base2_4hr_k_max.k1;
    base2_4hr_k_tuned_green_all = [base2_4hr_k_tuned_green_all base2_4hr_k_tuned_green];
    base2_4hr_max_tuned_green = base2_4hr_k_max.max_dfof;
    base2_4hr_max_tuned_green_all = [base2_4hr_max_tuned_green_all base2_4hr_max_tuned_green];
    base2_4hr_ncells = length(base2_4hr_k_tuned_green);
    base2_4hr_mouse_id = string(repmat(mouse,base2_4hr_ncells,1))';
    base2_4hr_mouse_id_all = [base2_4hr_mouse_id_all base2_4hr_mouse_id];
    base2_4hr_avgtuningcurve = base2_4hr_newavg.newAvg(base2_4hr_tuned,:)';
    base2_4hr_avgtuningcurve_all = [base2_4hr_avgtuningcurve_all base2_4hr_avgtuningcurve];
end

base2_4hr_avgtuningcurve_all = [base2_4hr_avgtuningcurve_all]';

base2_4hr_tuning_width_all = rad2deg(acos((1./base2_4hr_k_tuned_green_all) .* log((exp(base2_4hr_k_tuned_green_all)./2) + (exp(-base2_4hr_k_tuned_green_all)./2))));

%%
sz = size(base2_4hr_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if base2_4hr_avgtuningcurve_all(i,j) > 0.5
            base2_4hr_avgtuningcurve_all(i,:) = [];
        elseif base2_4hr_avgtuningcurve_all(i,j) < -0.5
            base2_4hr_avgtuningcurve_all(i,:) = [];
        end
    end
end 




%%
%postdark + 15min
post_dark_ori_all = [];
post_dark_k_max_all = [];
post_dark_k_tuned_green_all = [];
post_dark_max_tuned_green_all = [];
post_dark_mouse_id_all = [];
post_dark_avgtuningcurve_all = [];


post_dark = [5 5+7 5+14];

for isess = post_dark
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_dark_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_dark_ori_all = [post_dark_ori_all post_dark_ori];
    post_dark_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_dark_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_dark_tuned = post_dark_ori.ind_theta90; %only taking tuned cells
    post_dark_k_max.k1 = post_dark_k_max.k1(post_dark_tuned);
    post_dark_k_max.max_dfof = post_dark_k_max.max_dfof(post_dark_tuned);
    post_dark_k_max_all = [post_dark_k_max_all post_dark_k_max];
    post_dark_k_tuned_green = post_dark_k_max.k1;
    post_dark_k_tuned_green_all = [post_dark_k_tuned_green_all post_dark_k_tuned_green];
    post_dark_max_tuned_green = post_dark_k_max.max_dfof;
    post_dark_max_tuned_green_all = [post_dark_max_tuned_green_all post_dark_max_tuned_green];
    post_dark_ncells = length(post_dark_k_tuned_green);
    post_dark_mouse_id = string(repmat(mouse,post_dark_ncells,1))';
    post_dark_mouse_id_all = [post_dark_mouse_id_all post_dark_mouse_id];
    post_dark_avgtuningcurve = post_dark_newavg.newAvg(post_dark_tuned,:)';
    post_dark_avgtuningcurve_all = [post_dark_avgtuningcurve_all post_dark_avgtuningcurve];
end

post_dark_avgtuningcurve_all = [post_dark_avgtuningcurve_all]';

post_dark_tuning_width_all = rad2deg(acos((1./post_dark_k_tuned_green_all) .* log((exp(post_dark_k_tuned_green_all)./2) + (exp(-post_dark_k_tuned_green_all)./2))));

%%
sz = size(post_dark_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if post_dark_avgtuningcurve_all(i,j) > 0.5
            post_dark_avgtuningcurve_all(i,:) = [];
        elseif post_dark_avgtuningcurve_all(i,j) < -0.5
            post_dark_avgtuningcurve_all(i,:) = [];
        end
    end
end 



%%
%post 4hr light
post_light_ori_all = [];
post_light_k_max_all = [];
post_light_k_tuned_green_all = [];
post_light_max_tuned_green_all = [];
post_light_mouse_id_all = [];
post_light_avgtuningcurve_all = [];


post_light = [6 6+7 6+14];

for isess = post_light
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_light_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_light_ori_all = [post_light_ori_all post_light_ori];
    post_light_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_light_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_light_tuned = post_light_ori.ind_theta90; %only taking tuned cells
    post_light_k_max.k1 = post_light_k_max.k1(post_light_tuned);
    post_light_k_max.max_dfof = post_light_k_max.max_dfof(post_light_tuned);
    post_light_k_max_all = [post_light_k_max_all post_light_k_max];
    post_light_k_tuned_green = post_light_k_max.k1;
    post_light_k_tuned_green_all = [post_light_k_tuned_green_all post_light_k_tuned_green];
    post_light_max_tuned_green = post_light_k_max.max_dfof;
    post_light_max_tuned_green_all = [post_light_max_tuned_green_all post_light_max_tuned_green];
    post_light_ncells = length(post_light_k_tuned_green);
    post_light_mouse_id = string(repmat(mouse,post_light_ncells,1))';
    post_light_mouse_id_all = [post_light_mouse_id_all post_light_mouse_id];
    post_light_avgtuningcurve = post_light_newavg.newAvg(post_light_tuned,:)';
    post_light_avgtuningcurve_all = [post_light_avgtuningcurve_all post_light_avgtuningcurve];
end

post_light_avgtuningcurve_all = [post_light_avgtuningcurve_all]';

post_light_tuning_width_all = rad2deg(acos((1./post_light_k_tuned_green_all) .* log((exp(post_light_k_tuned_green_all)./2) + (exp(-post_light_k_tuned_green_all)./2))));

%%
sz = size(post_light_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if post_light_avgtuningcurve_all(i,j) > 0.5
            post_light_avgtuningcurve_all(i,:) = [];
        elseif post_light_avgtuningcurve_all(i,j) < -0.5
            post_light_avgtuningcurve_all(i,:) = [];
        end
    end
end 



%%
%post 4hr light + 7d
post_light_7d_ori_all = [];
post_light_7d_k_max_all = [];
post_light_7d_k_tuned_green_all = [];
post_light_7d_max_tuned_green_all = [];
post_light_7d_mouse_id_all = [];
post_light_7d_avgtuningcurve_all = [];


post_light_7d = [7 7+7 7+14];

for isess = post_light_7d
    mouse = expt(isess).mouse;
    date = expt(isess).date;
    img_area = expt(isess).img_loc{1};
    img_layer = expt(isess).img_loc{2};
    ref_str = ['runs-' expt(isess).runs];
    post_light_7d_ori = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'oriTuningInfo.mat']));
    post_light_7d_ori_all = [post_light_7d_ori_all post_light_7d_ori];
    post_light_7d_k_max = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_' 'k_and_max_vals.mat']));
    post_light_7d_newavg = load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_newAvg.mat']));
    post_light_7d_tuned = post_light_7d_ori.ind_theta90; %only taking tuned cells
    post_light_7d_k_max.k1 = post_light_7d_k_max.k1(post_light_7d_tuned);
    post_light_7d_k_max.max_dfof = post_light_7d_k_max.max_dfof(post_light_7d_tuned);
    post_light_7d_k_max_all = [post_light_7d_k_max_all post_light_7d_k_max];
    post_light_7d_k_tuned_green = post_light_7d_k_max.k1;
    post_light_7d_k_tuned_green_all = [post_light_7d_k_tuned_green_all post_light_7d_k_tuned_green];
    post_light_7d_max_tuned_green = post_light_7d_k_max.max_dfof;
    post_light_7d_max_tuned_green_all = [post_light_7d_max_tuned_green_all post_light_7d_max_tuned_green];
    post_light_7d_ncells = length(post_light_7d_k_tuned_green);
    post_light_7d_mouse_id = string(repmat(mouse,post_light_7d_ncells,1))';
    post_light_7d_mouse_id_all = [post_light_7d_mouse_id_all post_light_7d_mouse_id];
    post_light_7d_avgtuningcurve = post_light_7d_newavg.newAvg(post_light_7d_tuned,:)';
    post_light_7d_avgtuningcurve_all = [post_light_7d_avgtuningcurve_all post_light_7d_avgtuningcurve];
end

post_light_7d_avgtuningcurve_all = [post_light_7d_avgtuningcurve_all]';

post_light_7d_tuning_width_all = rad2deg(acos((1./post_light_7d_k_tuned_green_all) .* log((exp(post_light_7d_k_tuned_green_all)./2) + (exp(-post_light_7d_k_tuned_green_all)./2))));

%%
sz = size(post_light_7d_avgtuningcurve_all);
sz1 = sz(1);
sz2 = sz(2);

for i = 1:sz(1)
    for j = 1:sz(2)
        if post_light_7d_avgtuningcurve_all(i,j) > 0.5
            post_light_7d_avgtuningcurve_all(i,:) = [];
        elseif post_light_7d_avgtuningcurve_all(i,j) < -0.5
            post_light_7d_avgtuningcurve_all(i,:) = [];
        end
    end
end 



%%
figure; 
h=cdfplot(base1_k_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_k_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_k_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_k_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_k_tuned_green_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_k_tuned_green_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_k_tuned_green_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('k Values')
xlim([0 30]);
title('K Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_k.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(base1_max_tuned_green_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_max_tuned_green_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_max_tuned_green_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_max_tuned_green_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_max_tuned_green_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_max_tuned_green_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_max_tuned_green_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values')
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_max.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(base1_tuning_width_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all);
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_tuning_width.pdf']), '-dpdf', '-bestfit')


%%
%avg tuning curves
figure;
fast_errbar(1:8,base1_avgtuningcurve_all,1,'color',[0 0 0]);hold on;fast_errbar(1:8,base1_4d_avgtuningcurve_all,1,'color',[.15 .15 .15]);fast_errbar(1:8,base2_avgtuningcurve_all,1,'color',[.3 .3 .3]); fast_errbar(1:8,base2_4hr_avgtuningcurve_all,1,'color',[.45 .45 .45]); fast_errbar(1:8,post_dark_avgtuningcurve_all,1,'color',[.6 .6 .6]); fast_errbar(1:8,post_light_avgtuningcurve_all,1,'color',[.75 .75 .75]); fast_errbar(1:8,post_light_7d_avgtuningcurve_all,1,'color',[.9 .9 .9]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'])
hold off;
print(fullfile(realfnout, ['pooldarklight_allpoints_avgtuning.pdf']), '-dpdf', '-bestfit')

%%
prefori_scores = [];
pref_dscores_all = [];
k_scores = [];
max_scores = [];

list = [1 1+7 1+14];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    prefori_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'ori_changes']));
    prefori_scores = [prefori_scores prefori_ses_all];
    pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'd_scores']));
    pref_dscores_all = [pref_dscores_all pref_dscores_ses_all];
%     k_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
%     k_scores = [k_scores k_ses_all];
%     max_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
%     max_scores = [max_scores max_ses_all];
end

%%
pref_d_1_2_all = [];
pref_d_1_3_all = [];
pref_d_3_4_all = [];
pref_d_4_5_all = [];
pref_d_5_6_all = [];
pref_d_5_7_all = [];
mouse_id_1_2_all = [];


for idata = 1:length(pref_dscores_all)
    pref_d_1_2 = pref_dscores_all(idata).d_score_prefori_1_2_match;
    pref_d_1_2_all = [pref_d_1_2_all pref_d_1_2];
    pref_d_1_3 = pref_dscores_all(idata).d_score_prefori_1_3_match;
    pref_d_1_3_all = [pref_d_1_3_all pref_d_1_3];
    pref_d_3_4 = pref_dscores_all(idata).d_score_prefori_3_4_match;
    pref_d_3_4_all = [pref_d_3_4_all pref_d_3_4];
    pref_d_4_5 = pref_dscores_all(idata).d_score_prefori_4_5_match;
    pref_d_4_5_all = [pref_d_4_5_all pref_d_4_5];
    pref_d_5_6 = pref_dscores_all(idata).d_score_prefori_5_6_match;
    pref_d_5_6_all = [pref_d_5_6_all pref_d_5_6];
    pref_d_5_7 = pref_dscores_all(idata).d_score_prefori_5_7_match;
    pref_d_5_7_all = [pref_d_5_7_all pref_d_5_7];
end


%%
figure; 
h=cdfplot(pref_d_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(pref_d_1_3_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(pref_d_3_4_all);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(pref_d_4_5_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(pref_d_5_6_all);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(pref_d_5_7_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
xlabel('Change in Pref Ori')
title('Change in Pref Ori');
legend(['Baseline1 to Baseline1 + 4d'], ['Baseline1 to Baseline2'], ['Baseline2 to Baseline2 + 4hr'], ['Baseline2 + 4hr to Post Dark'], ['Post Dark to Post Light 4hr'], ['Post Dark to Post Light 7d'])
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%Jeon paper comparison

%img_point: 1 = baseline1, 2 = baseline1 + 4d, 3 = baseline2, 4 = baseline2
%+ 4hrs, 5 = post 4d dark, 6 = post dark + 4hrs, 7 = post dark + 7d

figure; 
% h=cdfplot(pref_d_1_2_all);
% set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
j=cdfplot(pref_d_1_3_all);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
hold on
% k=cdfplot(pref_d_3_4_all);
% set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(pref_d_4_5_all);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
% m=cdfplot(pref_d_5_6_all);
% set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(pref_d_5_7_all);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
xlabel('Change in Pref Ori')
title('Change in Pref Ori');
legend(['Baseline1 to Baseline2'], ['Baseline2 to Post Dark'], ['Post Dark to Post Light 7d'], 'Location', 'best')
hold off
print(fullfile(realfnout, ['pooldarklight_allpoints_Jeon_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

%%
%2 session comparisons using only cells matched to those days

comparison_ids_all = [];
comparison_pref_dscores_all = [];
comparison_max_all = [];

list = [1 1+7 1+14];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    comparison_ids = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'id_matches']));
    %comparison_ids_all = [comparison_ids_all comparison_ids];
    comparison_pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'same_cell_dscores']));
    comparison_pref_dscores_all = [comparison_pref_dscores_all comparison_pref_dscores_ses_all];
    comparison_max_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'same_cell_max_vals']));
    comparison_max_all = [comparison_max_all comparison_max_ses_all];


%     k_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_k_changes']));
%     k_scores = [k_scores k_ses_all];
%     max_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'RW_max_changes']));
%     max_scores = [max_scores max_ses_all];
end

%%

comparison_max_1_2_all = [];
comparison_max_2_1_all = [];
comparison_max_1_3_all = [];
comparison_max_3_1_all = [];
comparison_max_3_4_all = [];
comparison_max_4_3_all = [];
comparison_max_4_5_all = [];
comparison_max_5_4_all = [];
comparison_max_5_6_all = [];
comparison_max_6_5_all = [];
comparison_max_5_7_all = [];
comparison_max_7_5_all = [];


for idata = 1:length(comparison_max_all)
    comparison_max_1_2 = comparison_max_all(idata).comparison_max_1_2_match_tune;
    comparison_max_1_2_all = [comparison_max_1_2_all comparison_max_1_2];
    comparison_max_2_1 = comparison_max_all(idata).comparison_max_2_1_match_tune;
    comparison_max_2_1_all = [comparison_max_2_1_all comparison_max_2_1];

    comparison_max_1_3 = comparison_max_all(idata).comparison_max_1_3_match_tune;
    comparison_max_1_3_all = [comparison_max_1_3_all comparison_max_1_3];
    comparison_max_3_1 = comparison_max_all(idata).comparison_max_3_1_match_tune;
    comparison_max_3_1_all = [comparison_max_3_1_all comparison_max_3_1];

    comparison_max_3_4 = comparison_max_all(idata).comparison_max_3_4_match_tune;
    comparison_max_3_4_all = [comparison_max_3_4_all comparison_max_3_4];
    comparison_max_4_3 = comparison_max_all(idata).comparison_max_4_3_match_tune;
    comparison_max_4_3_all = [comparison_max_4_3_all comparison_max_4_3];

    comparison_max_4_5 = comparison_max_all(idata).comparison_max_4_5_match_tune;
    comparison_max_4_5_all = [comparison_max_4_5_all comparison_max_4_5];
    comparison_max_5_4 = comparison_max_all(idata).comparison_max_5_4_match_tune;
    comparison_max_5_4_all = [comparison_max_5_4_all comparison_max_5_4];

    comparison_max_5_6 = comparison_max_all(idata).comparison_max_5_6_match_tune;
    comparison_max_5_6_all = [comparison_max_5_6_all comparison_max_5_6];
    comparison_max_6_5 = comparison_max_all(idata).comparison_max_6_5_match_tune;
    comparison_max_6_5_all = [comparison_max_6_5_all comparison_max_6_5];

    comparison_max_5_7 = comparison_max_all(idata).comparison_max_5_7_match_tune;
    comparison_max_5_7_all = [comparison_max_5_7_all comparison_max_5_7];
    comparison_max_7_5 = comparison_max_all(idata).comparison_max_7_5_match_tune;
    comparison_max_7_5_all = [comparison_max_7_5_all comparison_max_7_5];
end

%%

fig = figure;
%sgtitle('Max dF/F')
%a = subplot(2,2,1)
h=cdfplot(comparison_max_1_3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_3_1_all);
set(j, 'LineStyle', '--', 'Color', [.0 .0 .0], 'LineWidth',2);
m=cdfplot(comparison_max_5_7_all);
set(m, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
n=cdfplot(comparison_max_7_5_all);
set(n, 'LineStyle', '--', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base1'], ['Base2'], ['Post-dark'], ['Post-dark + 7d'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['7-day comparisons'])
hold off
print(fullfile(realfnout, ['pooled_comparison_max_7d.pdf']), '-dpdf', '-bestfit')


fig = figure;
%sgtitle('Max dF/F')
%a = subplot(2,2,1)
h=cdfplot(comparison_max_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_2_1_all);
set(j, 'LineStyle', '--', 'Color', [.0 .0 .0], 'LineWidth',2);
m=cdfplot(comparison_max_4_5_all);
set(m, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
n=cdfplot(comparison_max_5_4_all);
set(n, 'LineStyle', '--', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base1'], ['Base1 + 4d'], ['Base2 + 4hr'], ['Post-dark'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['4-day comparisons'])
hold off
print(fullfile(realfnout, ['pooled_comparison_max_4d.pdf']), '-dpdf', '-bestfit')


fig = figure;
%sgtitle('Max dF/F')
%a = subplot(2,2,1)
h=cdfplot(comparison_max_3_4_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_4_3_all);
set(j, 'LineStyle', '--', 'Color', [.0 .0 .0], 'LineWidth',2);
m=cdfplot(comparison_max_5_6_all);
set(m, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
n=cdfplot(comparison_max_6_5_all);
set(n, 'LineStyle', '--', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base2'], ['Base2 + 4hr'], ['Post-dark'], ['Post-dark + 4hr'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['4-hour comparisons'])
hold off
print(fullfile(realfnout, ['pooled_comparison_max_4hr.pdf']), '-dpdf', '-bestfit')

%%
%ks tests for above
[ks_1_3 ks_p_1_3] = kstest2(comparison_max_1_3_all, comparison_max_3_1_all)
[ks_5_7 ks_p_5_7] = kstest2(comparison_max_5_7_all, comparison_max_7_5_all)
[ks_1_2 ks_p_1_2] = kstest2(comparison_max_1_2_all, comparison_max_2_1_all)
[ks_4_5 ks_p_4_5] = kstest2(comparison_max_4_5_all, comparison_max_5_4_all)
[ks_3_4 ks_p_3_4] = kstest2(comparison_max_3_4_all, comparison_max_4_3_all)
[ks_5_6 ks_p_5_6] = kstest2(comparison_max_5_6_all, comparison_max_6_5_all)


%%

comparison_pref_d_1_2_all = [];
comparison_pref_d_1_3_all = [];
comparison_pref_d_3_4_all = [];
comparison_pref_d_4_5_all = [];
comparison_pref_d_5_6_all = [];
comparison_pref_d_5_7_all = [];

comparison_k_d_1_2_all = [];
comparison_k_d_1_3_all = [];
comparison_k_d_3_4_all = [];
comparison_k_d_4_5_all = [];
comparison_k_d_5_6_all = [];
comparison_k_d_5_7_all = [];

comparison_max_d_1_2_all = [];
comparison_max_d_1_3_all = [];
comparison_max_d_3_4_all = [];
comparison_max_d_4_5_all = [];
comparison_max_d_5_6_all = [];
comparison_max_d_5_7_all = [];




for idata = 1:length(comparison_pref_dscores_all)
    comparison_pref_d_1_2 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_1_2;
    comparison_pref_d_1_2_all = [comparison_pref_d_1_2_all comparison_pref_d_1_2];
    comparison_pref_d_1_3 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_1_3;
    comparison_pref_d_1_3_all = [comparison_pref_d_1_3_all comparison_pref_d_1_3];
    comparison_pref_d_3_4 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_3_4;
    comparison_pref_d_3_4_all = [comparison_pref_d_3_4_all comparison_pref_d_3_4];
    comparison_pref_d_4_5 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_4_5;
    comparison_pref_d_4_5_all = [comparison_pref_d_4_5_all comparison_pref_d_4_5];
    comparison_pref_d_5_6 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_5_6;
    comparison_pref_d_5_6_all = [comparison_pref_d_5_6_all comparison_pref_d_5_6];
    comparison_pref_d_5_7 = comparison_pref_dscores_all(idata).comparison_d_score_prefori_5_7;
    comparison_pref_d_5_7_all = [comparison_pref_d_5_7_all comparison_pref_d_5_7];

    comparison_k_d_1_2 = comparison_pref_dscores_all(idata).comparison_d_score_k_1_2;
    comparison_k_d_1_2_all = [comparison_k_d_1_2_all comparison_k_d_1_2];
    comparison_k_d_1_3 = comparison_pref_dscores_all(idata).comparison_d_score_k_1_3;
    comparison_k_d_1_3_all = [comparison_k_d_1_3_all comparison_k_d_1_3];
    comparison_k_d_3_4 = comparison_pref_dscores_all(idata).comparison_d_score_k_3_4;
    comparison_k_d_3_4_all = [comparison_k_d_3_4_all comparison_k_d_3_4];
    comparison_k_d_4_5 = comparison_pref_dscores_all(idata).comparison_d_score_k_4_5;
    comparison_k_d_4_5_all = [comparison_k_d_4_5_all comparison_k_d_4_5];
    comparison_k_d_5_6 = comparison_pref_dscores_all(idata).comparison_d_score_k_5_6;
    comparison_k_d_5_6_all = [comparison_k_d_5_6_all comparison_k_d_5_6];
    comparison_k_d_5_7 = comparison_pref_dscores_all(idata).comparison_d_score_k_5_7;
    comparison_k_d_5_7_all = [comparison_k_d_5_7_all comparison_k_d_5_7];

    comparison_max_d_1_2 = comparison_pref_dscores_all(idata).comparison_d_score_max_1_2;
    comparison_max_d_1_2_all = [comparison_max_d_1_2_all comparison_max_d_1_2];
    comparison_max_d_1_3 = comparison_pref_dscores_all(idata).comparison_d_score_max_1_3;
    comparison_max_d_1_3_all = [comparison_max_d_1_3_all comparison_max_d_1_3];
    comparison_max_d_3_4 = comparison_pref_dscores_all(idata).comparison_d_score_max_3_4;
    comparison_max_d_3_4_all = [comparison_max_d_3_4_all comparison_max_d_3_4];
    comparison_max_d_4_5 = comparison_pref_dscores_all(idata).comparison_d_score_max_4_5;
    comparison_max_d_4_5_all = [comparison_max_d_4_5_all comparison_max_d_4_5];
    comparison_max_d_5_6 = comparison_pref_dscores_all(idata).comparison_d_score_max_5_6;
    comparison_max_d_5_6_all = [comparison_max_d_5_6_all comparison_max_d_5_6];
    comparison_max_d_5_7 = comparison_pref_dscores_all(idata).comparison_d_score_max_5_7;
    comparison_max_d_5_7_all = [comparison_max_d_5_7_all comparison_max_d_5_7];



    
end


%%
%subplot comparisons of two sessions - need to add legend, titles, etc.
%also work on adding k and max changes

fig=figure;
sgtitle('Pref Ori Changes')
a = subplot(2,2,1)
h=cdfplot(comparison_pref_d_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_pref_d_4_5_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_pref_d_1_2_all)), ')'], ['Base2 + 4hr - Post Dark (n= ', num2str(length(comparison_pref_d_4_5_all)), ')'], 'Location', 'best')
xlim([0 90])
xlabel([])
ylabel([])
title([])
hold off


b = subplot(2,2,2)
h=cdfplot(comparison_pref_d_3_4_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_pref_d_5_6_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_pref_d_3_4_all)), ')'], ['Post Dark - Post Dark + 4hr Light (n= ', num2str(length(comparison_pref_d_5_6_all)), ')'], 'Location', 'best')
xlabel([])
ylabel([])
title([])
hold off

c = subplot(2,2,3)
h=cdfplot(comparison_pref_d_1_3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_pref_d_5_7_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base1 - Base2 (n= ', num2str(length(comparison_pref_d_1_3_all)), ')'], ['Post Dark - Post Dark + 7d Light (n= ', num2str(length(comparison_pref_d_5_7_all)), ')'], 'Location', 'best')
xlabel([])
ylabel([])
title([])
hold off

% d = subplot(2,2,4)
% h=cdfplot(extra_comparison_pref_d_1_3_all);
% set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
% hold on
% j=cdfplot(extra_comparison_pref_d_6_7_all);
% set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
% legend(['Base1 - Base2 (n= ', num2str(length(extra_comparison_pref_d_1_3_all)), ')'], ['Post Dark + 4hr - Post Dark + 7d Light (n= ', num2str(length(extra_comparison_pref_d_6_7_all)), ')'], 'Location', 'best')
% xlabel([])
% ylabel([])
% title([])
% hold off

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c d],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in Pref Ori');


print(fullfile(realfnout, ['pooled_comparison_prefori_change.pdf']), '-dpdf', '-fillpage')


%%
fig = figure;
sgtitle('k Changes')
a = subplot(2,2,1)
h=cdfplot(comparison_k_d_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_k_d_4_5_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_k_d_1_2_all)), ')'], ['Base2 - Post Dark (n= ', num2str(length(comparison_k_d_4_5_all)), ')'], 'Location', 'best')
%xlim([-30 30])
xlabel([])
ylabel([])
title([])
hold off


b = subplot(2,2,2)
h=cdfplot(comparison_k_d_3_4_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_k_d_5_6_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_k_d_3_4_all)), ')'], ['Post Dark - Post Light 4hr (n= ', num2str(length(comparison_k_d_5_6_all)), ')'], 'Location', 'best')
xlabel([])
ylabel([])
title([])
hold off

c = subplot(2,2,3)
h=cdfplot(comparison_k_d_1_3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_k_d_5_7_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base1 - Base2 (n= ', num2str(length(comparison_k_d_1_3_all)), ')'], ['Post Light 4hr - Post Light 7d (n= ', num2str(length(comparison_k_d_5_7_all)), ')'], 'Location', 'best')
hold off
xlabel([])
ylabel([])
title([])

% d = subplot(2,2,4)
% h=cdfplot(extra_comparison_k_d_1_3_all);
% set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
% hold on
% j=cdfplot(extra_comparison_k_d_6_7_all);
% set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
% %legend(['Base1 - Base2 (n= ', num2str(length(comparison_k_d_1_3_all)), ')'], ['Post Light 4hr - Post Light 7d (n= ', num2str(length(comparison_k_d_5_7_all)), ')'], 'Location', 'best')
% hold off
% xlabel([])
% ylabel([])
% title([])

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% han.XLim([-30 30])
linkaxes([a b c d],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in k Value');
xlim(han, [-30 30]);

print(fullfile(realfnout, ['pooled_comparison_k_change.pdf']), '-dpdf', '-fillpage')


%%


fig = figure;
sgtitle('Max dF/F Changes')
a = subplot(2,2,1)
h=cdfplot(comparison_max_d_1_2_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_d_4_5_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base1 - Base2 + 4d (n= ', num2str(length(comparison_max_d_1_2_all)), ')'], ['Base2 - Post Dark (n= ', num2str(length(comparison_max_d_4_5_all)), ')'], 'Location', 'southeast')
xlim([-1 1])
xlabel([])
ylabel([])
title([])
hold off


b = subplot(2,2,2)
h=cdfplot(comparison_max_d_3_4_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_d_5_6_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base2 - Base2 + 4hr (n= ', num2str(length(comparison_max_d_3_4_all)), ')'], ['Post Dark - Post Light 4hr (n= ', num2str(length(comparison_max_d_5_6_all)), ')'], 'Location', 'southeast')
xlabel([])
ylabel([])
title([])
hold off

c = subplot(2,2,3)
h=cdfplot(comparison_max_d_1_3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(comparison_max_d_5_7_all);
set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
%legend(['Base1 - Base2 (n= ', num2str(length(comparison_max_d_1_3_all)), ')'], ['Post Light 4hr - Post Light 7d (n= ', num2str(length(comparison_max_d_5_7_all)), ')'], 'Location', 'southeast')
hold off
xlabel([])
ylabel([])
title([])

% d = subplot(2,2,4)
% h=cdfplot(extra_comparison_max_d_1_3_all);
% set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
% hold on
% j=cdfplot(extra_comparison_max_d_6_7_all);
% set(j, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
% %legend(['Base1 - Base2 (n= ', num2str(length(comparison_max_d_1_3_all)), ')'], ['Post Light 4hr - Post Light 7d (n= ', num2str(length(comparison_max_d_5_7_all)), ')'], 'Location', 'southeast')
% hold off
% xlabel([])
% ylabel([])
% title([])

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
linkaxes([a b c d],'xy')
ylabel(han,'% of Cells');
xlabel(han,'Change in Max dF/F Value');

print(fullfile(realfnout, ['pooled_comparison_max_change.pdf']), '-dpdf', '-fillpage')


%%
%could make scatters of the above info too

%%
k_s_pref_1 = kstest2(comparison_pref_d_1_2_all, comparison_pref_d_4_5_all)
k_s_pref_2 = kstest2(comparison_pref_d_3_4_all, comparison_pref_d_5_6_all)
k_s_pref_3 = kstest2(comparison_pref_d_1_3_all, comparison_pref_d_5_7_all)

k_s_k_1 = kstest2(comparison_k_d_1_2_all, comparison_k_d_4_5_all)
k_s_k_2 = kstest2(comparison_k_d_3_4_all, comparison_k_d_5_6_all)
k_s_k_3 = kstest2(comparison_k_d_1_3_all, comparison_k_d_5_7_all)

k_s_max_1 = kstest2(comparison_max_d_1_2_all, comparison_max_d_4_5_all)
k_s_max_2 = kstest2(comparison_max_d_3_4_all, comparison_max_d_5_6_all)
k_s_max_3 = kstest2(comparison_max_d_1_3_all, comparison_max_d_5_7_all)

%max 1 and 3 are significant -> dark housing increases max; light
%reintroduction for 1 week decreases it

%%
%extra comparison
extra_comparison_pref_dscores_all = [];

list = [1 1+7 1+14];
for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    extra_comparison_pref_dscores_ses_all = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_' 'extra_comp']));
    extra_comparison_pref_dscores_all = [extra_comparison_pref_dscores_all extra_comparison_pref_dscores_ses_all];
end

extra_comparison_pref_d_1_3_all = [];
extra_comparison_pref_d_6_7_all = [];

extra_comparison_k_d_1_3_all = [];
extra_comparison_k_d_6_7_all = [];

extra_comparison_max_d_1_3_all = [];
extra_comparison_max_d_6_7_all = [];

%%

extra_comparison_max_1_3_all = [];
extra_comparison_max_3_1_all = [];
extra_comparison_max_6_7_all = [];
extra_comparison_max_7_6_all = [];


for idata = 1:length(extra_comparison_pref_dscores_all)
 
    %d scores
    extra_comparison_pref_d_1_3 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_prefori_1_3;
    extra_comparison_pref_d_1_3_all = [extra_comparison_pref_d_1_3_all extra_comparison_pref_d_1_3];
    extra_comparison_pref_d_6_7 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_prefori_6_7;
    extra_comparison_pref_d_6_7_all = [extra_comparison_pref_d_6_7_all extra_comparison_pref_d_6_7];

    
    extra_comparison_k_d_1_3 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_k_1_3;
    extra_comparison_k_d_1_3_all = [extra_comparison_k_d_1_3_all extra_comparison_k_d_1_3];
    extra_comparison_k_d_6_7 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_k_6_7;
    extra_comparison_k_d_6_7_all = [extra_comparison_k_d_6_7_all extra_comparison_k_d_6_7];

   
    extra_comparison_max_d_1_3 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_max_1_3;
    extra_comparison_max_d_1_3_all = [extra_comparison_max_d_1_3_all extra_comparison_max_d_1_3];
    extra_comparison_max_d_6_7 = extra_comparison_pref_dscores_all(idata).extra_comparison_d_score_max_6_7;
    extra_comparison_max_d_6_7_all = [extra_comparison_max_d_6_7_all extra_comparison_max_d_6_7];
    
    %regular max scores
    extra_comparison_max_1_3 = extra_comparison_pref_dscores_all(idata).extra_comparison_max_1_3_match_tune;
    extra_comparison_max_1_3_all = [extra_comparison_max_1_3_all extra_comparison_max_1_3];
    extra_comparison_max_3_1 = extra_comparison_pref_dscores_all(idata).extra_comparison_max_3_1_match_tune;
    extra_comparison_max_3_1_all = [extra_comparison_max_3_1_all extra_comparison_max_3_1];

    extra_comparison_max_6_7 = extra_comparison_pref_dscores_all(idata).extra_comparison_max_6_7_match_tune;
    extra_comparison_max_6_7_all = [extra_comparison_max_6_7_all extra_comparison_max_6_7];
    extra_comparison_max_7_6 = extra_comparison_pref_dscores_all(idata).extra_comparison_max_7_6_match_tune;
    extra_comparison_max_7_6_all = [extra_comparison_max_7_6_all extra_comparison_max_7_6];


end


%%
fig = figure;
%sgtitle('Max dF/F')
%a = subplot(2,2,1)
h=cdfplot(extra_comparison_max_1_3_all);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(extra_comparison_max_3_1_all);
set(j, 'LineStyle', '--', 'Color', [.0 .0 .0], 'LineWidth',2);
m=cdfplot(extra_comparison_max_6_7_all);
set(m, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
n=cdfplot(extra_comparison_max_7_6_all);
set(n, 'LineStyle', '--', 'Color', [.75 .75 .75], 'LineWidth',2);
legend(['Base1'], ['Base2'], ['Post-dark + 4hr'], ['Post-dark + 7d'], 'Location', 'southeast')
xlim([0 1])
xlabel(['Max dF/F'])
ylabel(['% of cells'])
title(['7-day (extra) comparisons'])
hold off
print(fullfile(realfnout, ['pooled_extra_comparison_max_7d.pdf']), '-dpdf', '-bestfit')

%%


%%
figure;
cdfplot(extra_comparison_max_d_1_3_all)
hold on
cdfplot(extra_comparison_max_d_6_7_all)

figure;
cdfplot(extra_comparison_k_d_1_3_all)
hold on
cdfplot(extra_comparison_k_d_6_7_all)

figure;
cdfplot(extra_comparison_pref_d_1_3_all)
hold on
cdfplot(extra_comparison_pref_d_6_7_all)
%%
%%
%%
figure;
subplot(3,2,1)
scatter(comparison_pref_d_1_2_all, comparison_max_d_1_2_all)
subplot(3,2,2)
scatter(comparison_pref_d_4_5_all, comparison_max_d_4_5_all)
subplot(3,2,3)
scatter(comparison_pref_d_3_4_all, comparison_max_d_3_4_all)
subplot(3,2,4)
scatter(comparison_pref_d_5_6_all, comparison_max_d_5_6_all)
subplot(3,2,5)
scatter(comparison_pref_d_1_3_all, comparison_max_d_1_3_all)
subplot(3,2,6)
scatter(comparison_pref_d_5_7_all, comparison_max_d_5_7_all)
%%
%indexing - use find
i2537_index_base1 = find(base1_mouse_id_all == 'i2537');
i2538_index_base1 = find(base1_mouse_id_all == 'i2538');
i2537_index_base1_4d = find(base1_4d_mouse_id_all == 'i2537');
i2538_index_base1_4d = find(base1_4d_mouse_id_all == 'i2538');
i2537_index_base2 = find(base2_mouse_id_all == 'i2537');
i2538_index_base2 = find(base2_mouse_id_all == 'i2538');
i2537_index_base2_4hr= find(base2_4hr_mouse_id_all == 'i2537');
i2538_index_base2_4hr = find(base2_4hr_mouse_id_all == 'i2538');
i2537_index_post_dark= find(post_dark_mouse_id_all == 'i2537');
i2538_index_post_dark = find(post_dark_mouse_id_all == 'i2538');
i2537_index_post_light= find(post_light_mouse_id_all == 'i2537');
i2538_index_post_light = find(post_light_mouse_id_all == 'i2538');
i2537_index_post_light_7d= find(post_light_7d_mouse_id_all == 'i2537');
i2538_index_post_light_7d = find(post_light_7d_mouse_id_all == 'i2538');

%%
%individual plots

figure; 
subplot(1,2,1)
h=cdfplot(base1_tuning_width_all(i2537_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all(i2537_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all(i2537_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all(i2537_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all(i2537_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all(i2537_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all(i2537_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')

 
subplot(1,2,2)
h=cdfplot(base1_tuning_width_all(i2538_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_tuning_width_all(i2538_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_tuning_width_all(i2538_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_tuning_width_all(i2538_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_tuning_width_all(i2538_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_tuning_width_all(i2538_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_tuning_width_all(i2538_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Tuning Width Values')
title('Tuning Width Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')
hold off
print(fullfile(realfnout, ['pooldarklight_individual_tuningwidth_cdf.pdf']), '-dpdf', '-bestfit')


figure; 
subplot(1,2,1)
h=cdfplot(base1_max_tuned_green_all(i2537_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_max_tuned_green_all(i2537_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_max_tuned_green_all(i2537_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_max_tuned_green_all(i2537_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_max_tuned_green_all(i2537_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_max_tuned_green_all(i2537_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_max_tuned_green_all(i2537_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values');
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')

 
subplot(1,2,2)
h=cdfplot(base1_max_tuned_green_all(i2538_index_base1));
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(base1_4d_max_tuned_green_all(i2538_index_base1_4d));
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(base2_max_tuned_green_all(i2538_index_base2));
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(base2_4hr_max_tuned_green_all(i2538_index_base2_4hr));
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(post_dark_max_tuned_green_all(i2538_index_post_dark));
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(post_light_max_tuned_green_all(i2538_index_post_light));
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
o=cdfplot(post_light_7d_max_tuned_green_all(i2538_index_post_light_7d));
set(o, 'LineStyle', '-', 'Color', [.9 .9 .9], 'LineWidth',2);
xlabel('Max dF/F Values');
xlim([0 1]);
title('Max dF/F Values');
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')
hold off
print(fullfile(realfnout, ['pooldarklight_individual_maxdfof_cdf.pdf']), '-dpdf', '-bestfit')



%%
%2537
figure;
subplot(1,2,1)
fast_errbar(1:8,base1_avgtuningcurve_all(i2537_index_base1,:),1,'color',[0 0 0]);
hold on;
fast_errbar(1:8,base1_4d_avgtuningcurve_all(i2537_index_base1_4d,:),1,'color',[.15 .15 .15]);
fast_errbar(1:8,base2_avgtuningcurve_all(i2537_index_base2,:),1,'color',[.3 .3 .3]); 
fast_errbar(1:8,base2_4hr_avgtuningcurve_all(i2537_index_base2_4hr,:),1,'color',[.45 .45 .45]); 
fast_errbar(1:8,post_dark_avgtuningcurve_all(i2537_index_post_dark,:),1,'color',[.6 .6 .6]); 
fast_errbar(1:8,post_light_avgtuningcurve_all(i2537_index_post_light,:),1,'color',[.75 .75 .75]); 
fast_errbar(1:8,post_light_7d_avgtuningcurve_all(i2537_index_post_light_7d,:),1,'color',[.9 .9 .9]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')
hold off

%2538
subplot(1,2,2)
fast_errbar(1:8,base1_avgtuningcurve_all(i2538_index_base1,:),1,'color',[0 0 0]);
hold on;
fast_errbar(1:8,base1_4d_avgtuningcurve_all(i2538_index_base1_4d,:),1,'color',[.15 .15 .15]);
fast_errbar(1:8,base2_avgtuningcurve_all(i2538_index_base2,:),1,'color',[.3 .3 .3]); 
fast_errbar(1:8,base2_4hr_avgtuningcurve_all(i2538_index_base2_4hr,:),1,'color',[.45 .45 .45]); 
fast_errbar(1:8,post_dark_avgtuningcurve_all(i2538_index_post_dark,:),1,'color',[.6 .6 .6]); 
fast_errbar(1:8,post_light_avgtuningcurve_all(i2538_index_post_light,:),1,'color',[.75 .75 .75]); 
fast_errbar(1:8,post_light_7d_avgtuningcurve_all(i2538_index_post_light_7d,:),1,'color',[.9 .9 .9]);
xlabel('Orientation')
ylabel('dF/F')
legend(['Baseline1'], ['Baseline1 + 4d'], ['Baseline2'], ['Baseline2 + 4hr'], ['Post Dark 15min'], ['Post Light 4hr'], ['Post Light + 7d'], 'Location','best')
hold off

print(fullfile(realfnout, ['pooldarklight_individual_avgtuning_cdf.pdf']), '-dpdf', '-bestfit')





%%
%cell ids for each matched condition and each mouse
i2537_match_cell_ids = [42 39 18 19 15 23];
i2538_match_cell_ids = [32 34 14 17 17 22];

pref_d_i2537_1_2 = pref_d_1_2_all(1:i2537_match_cell_ids(1));
pref_d_i2537_1_3 = pref_d_1_3_all(1:i2537_match_cell_ids(2));
pref_d_i2537_3_4 = pref_d_3_4_all(1:i2537_match_cell_ids(3));
pref_d_i2537_4_5 = pref_d_4_5_all(1:i2537_match_cell_ids(4));
pref_d_i2537_5_6 = pref_d_5_6_all(1:i2537_match_cell_ids(5));
pref_d_i2537_5_7 = pref_d_5_7_all(1:i2537_match_cell_ids(6));

pref_d_i2538_1_2 = pref_d_1_2_all(1:i2538_match_cell_ids(1));
pref_d_i2538_1_3 = pref_d_1_3_all(1:i2538_match_cell_ids(2));
pref_d_i2538_3_4 = pref_d_3_4_all(1:i2538_match_cell_ids(3));
pref_d_i2538_4_5 = pref_d_4_5_all(1:i2538_match_cell_ids(4));
pref_d_i2538_5_6 = pref_d_5_6_all(1:i2538_match_cell_ids(5));
pref_d_i2538_5_7 = pref_d_5_7_all(1:i2538_match_cell_ids(6));

%%
%these are wrong - it is basically the same cells for both mice
figure; 
h=cdfplot(pref_d_i2537_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(pref_d_i2537_1_3);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(pref_d_i2537_3_4);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(pref_d_i2537_4_5);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(pref_d_i2537_5_6);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(pref_d_i2537_5_7);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
xlabel('Change in Pref Ori')
title('Change in Pref Ori');
legend(['Baseline1 to Baseline1 + 4d'], ['Baseline1 to Baseline2'], ['Baseline2 to Baseline2 + 4hr'], ['Baseline2 + 4hr to Post Dark'], ['Post Dark to Post Light 4hr'], ['Post Dark to Post Light 7d'])
hold off
%print(fullfile(realfnout, ['pooldarklight_allpoints_change_pref_cdf.pdf']), '-dpdf', '-bestfit')

figure; 
h=cdfplot(pref_d_i2538_1_2);
set(h, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth',2);
hold on
j=cdfplot(pref_d_i2538_1_3);
set(j, 'LineStyle', '-', 'Color', [.15 .15 .15], 'LineWidth',2);
k=cdfplot(pref_d_i2538_3_4);
set(k, 'LineStyle', '-', 'Color', [.3 .3 .3], 'LineWidth',2);
l=cdfplot(pref_d_i2538_4_5);
set(l, 'LineStyle', '-', 'Color', [.45 .45 .45], 'LineWidth',2);
m=cdfplot(pref_d_i2538_5_6);
set(m, 'LineStyle', '-', 'Color', [.6 .6 .6], 'LineWidth',2);
n=cdfplot(pref_d_i2538_5_7);
set(n, 'LineStyle', '-', 'Color', [.75 .75 .75], 'LineWidth',2);
xlabel('Change in Pref Ori')
title('Change in Pref Ori');
legend(['Baseline1 to Baseline1 + 4d'], ['Baseline1 to Baseline2'], ['Baseline2 to Baseline2 + 4hr'], ['Baseline2 + 4hr to Post Dark'], ['Post Dark to Post Light 4hr'], ['Post Dark to Post Light 7d'])
hold off
%print(fullfile(realfnout, ['pooldarklight_allpoints_change_pref_cdf.pdf']), '-dpdf', '-bestfit')


