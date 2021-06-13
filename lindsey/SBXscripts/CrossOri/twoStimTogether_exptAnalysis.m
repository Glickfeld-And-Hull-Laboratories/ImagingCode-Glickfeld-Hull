clc; clear all; close all;
doRedChannel = 0;
ds = 'TwoStimTogether_ExptList';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);

iexp =14;

frame_rate = params.frameRate;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc;
coFolder = expt(iexp).coFolder;
nrun = length(coFolder);
run_str = catRunName(cell2mat(coFolder), nrun);
retFolder = expt(iexp).retFolder;
nrun = length(retFolder);
ret_run_str = catRunName(cell2mat(retFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fn_out = fullfile(LG_base, 'Analysis\2P');
fprintf([mouse ' ' date '\n'])

datemouse = [date '_' mouse];
datemouserun = [date '_' mouse '_' run_str];
datemouseretrun = [date '_' mouse '_' ret_run_str];

%% Pref direction analysis
load(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_TCs.mat']))
load(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_dataStim.mat']))
load(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_input.mat']))

%% cell dists
load(fullfile(fn_out, datemouse, datemouseretrun, [datemouseretrun '_lbub_fits.mat']));
fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')

% input stimulus location based on experimental choice
stimEl = (double(input.stimOneGratingElevationDeg)+double(input.stimTwoGratingElevationDeg))/2;
stimAz = (double(input.stimOneGratingAzimuthDeg)+double(input.stimTwoGratingAzimuthDeg))/2;
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])
cellAz = cell(1,3);
cellEl = cell(1,3);
cellDists = zeros(nCells,2,3);
dir_group = {'0','90','both'};
goodfit_dist_ind = [];
min_el_dist = 5;
min_az_dist = 2;
for i = 1:3
    fprintf([dir_group{i} '\n'])
    cellAz{i} = lbub_fits(:,4,4,i);
    cellEl{i} = lbub_fits(:,5,4,i);
    nC = length(goodfit_ind{i});
    fprintf(['# goodfit cells = ' num2str(nC) '\n'])
    
    % calculate cell distances
    fprintf('Calculating cell RF distances to stimulus...\n')
    cellDists(:,1,i) = abs(cellAz{i}-stimAz);
    cellDists(:,2,i) = abs(cellEl{i}-stimEl);
    
    goodDistIx = (cellDists(:,1,i)<=min_az_dist & cellDists(:,2,i)<=min_el_dist);
    n = length(intersect(find(goodDistIx),goodfit_ind{i}));
    fprintf(['# goodfit cells within ' num2str(min_az_dist) '/' num2str(min_el_dist) ' Az/El deg = ' num2str(n) '\n'])
    goodfit_dist_ind = [goodfit_dist_ind; intersect(find(goodDistIx),goodfit_ind{i})];
end
goodfit_dist_ind = unique(goodfit_dist_ind);

%%
if doRedChannel == 0
    red_cells = [];
end

prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
nFramesOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials-1);
data_f = nan(1,nCells,nTrials-1);

for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;

stimDir_all = celleqel2mat_padded(input.tStimOneGratingDirectionDeg);
maskDir_all = celleqel2mat_padded(input.tStimTwoGratingDirectionDeg);

ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));
nStim = 9;
ind_stim = cell(1,9);
stim_mat = cell(1,9);
stim_mat{1} = [{[1]},{[1]},{[1 2]},{[1 2]}];
stim_mat{2} = [{[1]},{[2]},{[1 2]},{[1]}];
stim_mat{3} = [{[1]},{[2]},{[1 2]},{[2]}];
stim_mat{4} = [{[2]},{[1]},{[1]},{[1 2]}];
stim_mat{5} = [{[2]},{[2]},{[1]},{[1]}];
stim_mat{6} = [{[2]},{[2]},{[1]},{[2]}];
stim_mat{7} = [{[2]},{[1]},{[2]},{[1 2]}];
stim_mat{8} = [{[2]},{[2]},{[2]},{[1]}];
stim_mat{9} = [{[2]},{[2]},{[2]},{[2]}];
stim_mat{10} = [{['stimCon']},{['maskCon']},{['stimDir']},{['maskDir']}];
    
resp_cell = cell(1,9);
resp_cell_tc = cell(1,9);
trialsperstim = zeros(1,9);
h_resp = zeros(nCells,9);
p_resp = zeros(nCells,9);

resp_win = prewin_frames+5:prewin_frames+nFramesOn;
base_win = 1+5:prewin_frames;
resp_blank = squeeze(mean(data_dfof_tc(prewin_frames+frame_rate:prewin_frames+nFramesOn,:,ind_blank),1));
avg_resp_stim = zeros(nCells,9,2);
avg_resp_tc = zeros(length(tt),nCells,9,2);

test_mat = cell(1,9);
for i=1:nStim
    test_mat{i} = cell(1,4);
    ind1 = [];
    ind2 = [];
    ind3 = [];
    ind4 = [];
    for ii = stim_mat{i}{1}
        ind1 = [ind1 find(stimCon_all == stimCons(ii))];
        test_mat{i}{1} = [test_mat{i}{1} stimCons(ii)];
    end
    for ii = stim_mat{i}{2}
        ind2 = [ind2 find(maskCon_all == maskCons(ii))];
        test_mat{i}{2} = [test_mat{i}{2} maskCons(ii)];
    end
    for ii = stim_mat{i}{3}
        ind3 = [ind3 find(stimDir_all == stimDirs(ii))];
        test_mat{i}{3} = [test_mat{i}{3} stimDirs(ii)];
    end
    for ii = stim_mat{i}{4}
        ind4 = [ind4 find(maskDir_all == maskDirs(ii))];
        test_mat{i}{4} = [test_mat{i}{4} maskDirs(ii)];
    end
    ind_stim{i} = intersect(intersect(intersect(ind1,ind2),ind3),ind4);
    trialsperstim(1,i) = length(ind_stim{i});
    resp_cell{i} = squeeze(mean(data_dfof_tc(resp_win,:,ind_stim{i}),1));
    base_cell{i} = squeeze(mean(data_dfof_tc(base_win,:,ind_stim{i}),1));
    resp_cell_tc{i} = squeeze(data_dfof_tc(:,:,ind_stim{i}));
    for iCell=1:nCells
        [h_resp(iCell,i), p_resp(iCell,i)] = ttest2(resp_cell{i}(iCell,:),base_cell{i}(iCell,:),'tail','right','alpha', 0.05/nStim);
    end
    avg_resp_stim(:,i,1) = squeeze(mean(mean(data_dfof_tc(resp_win,:,ind_stim{i}),1),3));
    avg_resp_stim(:,i,2) = squeeze(std(mean(data_dfof_tc(resp_win,:,ind_stim{i}),1),[],3)./sqrt(length(ind_stim{i})));
    avg_resp_tc(:,:,i,1) = mean(data_dfof_tc(:,:,ind_stim{i}),3);
    avg_resp_tc(:,:,i,2) = std(data_dfof_tc(:,:,ind_stim{i}),[],3)./sqrt(length(ind_stim{i}));
end

stim_alone = [2 3 4 7];
stim_same = [5 9];
stim_same_comp{1} = [2 4];
stim_same_comp{2} = [3 7];
stim_opp = [6 8];
stim_opp_comp{1} = [3 4];
stim_opp_comp{2} = [2 7];
resp_ind_alone = cell(1,length(stim_alone));
for i = 1:length(stim_alone)
    resp_ind_alone{i} = find(h_resp(:,stim_alone(i)));
end
resp_ind_dir = cell(1,2);
for i = 1:2
    resp_ind_dir{i} = find(sum(h_resp(:,stim_same_comp{i}),2));
end
resp_ind_all = find(sum(h_resp(:,stim_alone),2));

[h_pref_dir p_pref_dir] = ttest2([resp_cell{2} resp_cell{4} resp_cell{5}],[resp_cell{3} resp_cell{7} resp_cell{9}],'tail','both','dim',2);
pref_dir = mean([resp_cell{2} resp_cell{4} resp_cell{5}],2) > mean([resp_cell{3} resp_cell{7} resp_cell{9}],2);
ind_pref_dir{1} = intersect(resp_ind_dir{1},intersect(find(pref_dir==0),find(h_pref_dir)));
ind_pref_dir{2} = intersect(resp_ind_dir{2},intersect(find(pref_dir),find(h_pref_dir)));
ind_pref_dir{3} = intersect(resp_ind_all,find(h_pref_dir==0));

eq_resp{1} = intersect(ind_pref_dir{1},find(sum(h_resp(:,[2 4]),2)==2));
eq_resp{2} = intersect(ind_pref_dir{2},find(sum(h_resp(:,[3 7]),2)==2));
eq_resp{3} = unique([eq_resp{1}; eq_resp{2}]);


resp_tc_bystim{1} = mean(cat(3,resp_cell_tc{2}, resp_cell_tc{4}, resp_cell_tc{5}),3);
resp_tc_bystim{2} = mean(cat(3,resp_cell_tc{3}, resp_cell_tc{7}, resp_cell_tc{9}),3);

save(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_respData.mat']), 'resp_cell', 'data_dfof_tc', 'resp_ind_alone', 'resp_ind_all', 'resp_ind_dir', 'frame_rate', 'h_resp', 'avg_resp_stim', 'avg_resp_tc', 'h_pref_dir','pref_dir','ind_pref_dir', 'eq_resp','goodfit_dist_ind', 'cellDists');
save(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_stimData.mat']), 'prewin_frames', 'postwin_frames', 'resp_win', 'ind_stim', 'ind_blank','stim_mat');


%% plot all responsive cells

close all
for iCell = intersect(resp_ind_all,goodfit_dist_ind)'
    figure;
    subplot(3,3,1)
    plot(tt,resp_tc_bystim{1}(:,iCell))
    hold on
    plot(tt,resp_tc_bystim{2}(:,iCell))
    ylim([-0.2 0.7])
    for i = 1:9
        subplot(3,3,i)
        shadedErrorBar(tt,avg_resp_tc(:,iCell,i,1),avg_resp_tc(:,iCell,i,2));
        ylim([-0.2 0.7])
        if i==5
            hold on
            plot(tt,avg_resp_tc(:,iCell,2,1)+avg_resp_tc(:,iCell,4,1))
        end
        if i==6
            hold on
            plot(tt,avg_resp_tc(:,iCell,3,1)+avg_resp_tc(:,iCell,4,1))
        end
        if i==8
            hold on
            plot(tt,avg_resp_tc(:,iCell,2,1)+avg_resp_tc(:,iCell,7,1))
        end
        if i==9
            hold on
            plot(tt,avg_resp_tc(:,iCell,3,1)+avg_resp_tc(:,iCell,7,1))
        end
        if h_resp(iCell,i)
            title('**')
        end
    end
    suptitle(['Cell ' num2str(iCell) ' - Az: ' num2str(cellDists(iCell,1,3)) '; El: ' num2str(cellDists(iCell,2,3))])
    movegui('center')
    print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_respMatrix_cell' num2str(iCell) '.pdf']),'-dpdf','-bestfit')
end

%%
figure;
for i = 1:length(stim_same)
    ind = [];
    for ii = 1:length(stim_same_comp{i})
        iii = find(stim_alone==stim_same_comp{i}(ii));
        ind = [ind; resp_ind_alone{iii}];
    end
    subplot(2,2,i)
    scatter(avg_resp_stim(ind,stim_same_comp{i}(1),1)+avg_resp_stim(ind,stim_same_comp{i}(2),1),avg_resp_stim(ind,stim_same(i),1))
    title('Same')
    xlim([0 0.7])
    ylim([0 0.7])
    refline(1);
    xlabel('Left+Right')
    ylabel('Together')
    ind = [];
    for ii = 1:length(stim_opp_comp{i})
        iii = find(stim_alone==stim_opp_comp{i}(ii));
        ind = [ind; resp_ind_alone{iii}];
    end
    subplot(2,2,i+2)
    scatter(avg_resp_stim(ind,stim_opp_comp{i}(1),1)+avg_resp_stim(ind,stim_opp_comp{i}(2),1),avg_resp_stim(ind,stim_opp(i),1))
    title('Opp')
    xlim([0 0.7])
    ylim([0 0.7])
    refline(1);
    xlabel('Left+Right')
    ylabel('Together')
end
suptitle([date ' ' mouse])
print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_scatterAloneVsTogether.pdf']),'-dpdf','-bestfit')
    
avg_resp_stim_rect = avg_resp_stim;
avg_resp_stim_rect(find(avg_resp_stim<0)) = 0;

oppStim_SI{1} = (avg_resp_stim_rect(:,6,1)-(avg_resp_stim_rect(:,3,1)+avg_resp_stim_rect(:,4,1)))./(avg_resp_stim_rect(:,6,1)+(avg_resp_stim_rect(:,3,1)+avg_resp_stim_rect(:,4,1)));
oppStim_SI{2} = (avg_resp_stim_rect(:,8,1)-(avg_resp_stim_rect(:,2,1)+avg_resp_stim_rect(:,7,1)))./(avg_resp_stim_rect(:,8,1)+(avg_resp_stim_rect(:,2,1)+avg_resp_stim_rect(:,7,1)));
sameStim_SI{1} = (avg_resp_stim_rect(:,5,1)-(avg_resp_stim_rect(:,2,1)+avg_resp_stim_rect(:,4,1)))./(avg_resp_stim_rect(:,5,1)+(avg_resp_stim_rect(:,2,1)+avg_resp_stim_rect(:,4,1)));
sameStim_SI{2} = (avg_resp_stim_rect(:,9,1)-(avg_resp_stim_rect(:,3,1)+avg_resp_stim_rect(:,7,1)))./(avg_resp_stim_rect(:,9,1)+(avg_resp_stim_rect(:,3,1)+avg_resp_stim_rect(:,7,1)));
for iCell = 1:nCells
hopp_SI{1}(iCell) = ttest(resp_cell{6}(iCell,:), avg_resp_stim_rect(iCell,3,1)+avg_resp_stim_rect(iCell,4,1),'tail','both');
hopp_SI{2}(iCell) = ttest(resp_cell{8}(iCell,:), avg_resp_stim_rect(iCell,2,1)+avg_resp_stim_rect(iCell,7,1),'tail','both');
hsame_SI{1}(iCell) = ttest(resp_cell{5}(iCell,:), avg_resp_stim_rect(iCell,2,1)+avg_resp_stim_rect(iCell,4,1),'tail','both');
hsame_SI{2}(iCell) = ttest(resp_cell{9}(iCell,:), avg_resp_stim_rect(iCell,3,1)+avg_resp_stim_rect(iCell,7,1),'tail','both');
end

save(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_SI.mat']), 'oppStim_SI', 'sameStim_SI', 'hopp_SI', 'hsame_SI');

figure;
for i = 1:2
    subplot(2,2,i)
    histogram(sameStim_SI{i}(resp_ind_dir{i}),[-1:0.2:1])
    if intersect(resp_ind_dir{i},find(hsame_SI{i}))
        hold on
        histogram(sameStim_SI{i}(intersect(resp_ind_dir{i},find(hsame_SI{i}))),[-1:0.2:1])
    end
    title('Same')
    subplot(2,2,i+2)
     histogram(oppStim_SI{i}(resp_ind_dir{i}),[-1:0.2:1])
    if intersect(resp_ind_dir{i},find(hopp_SI{i}))
        hold on
        histogram(oppStim_SI{i}(intersect(resp_ind_dir{i},find(hopp_SI{i}))),[-1:0.2:1])
    end
    title('Opp')
end
suptitle([date ' ' mouse])
print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_SIhists.pdf']),'-dpdf','-bestfit')


%% by cell dists
figure;
for i = 1:length(stim_same)
    ind = [];
    for ii = 1:length(stim_same_comp{i})
        iii = find(stim_alone==stim_same_comp{i}(ii));
        ind = [ind; intersect(goodfit_dist_ind,resp_ind_alone{iii})];
    end
    subplot(2,2,i)
    scatter(avg_resp_stim(ind,stim_same_comp{i}(1),1)+avg_resp_stim(ind,stim_same_comp{i}(2),1),avg_resp_stim(ind,stim_same(i),1))
    title('Same')
    xlim([0 0.7])
    ylim([0 0.7])
    refline(1);
    xlabel('Left+Right')
    ylabel('Together')
    ind = [];
    for ii = 1:length(stim_opp_comp{i})
        iii = find(stim_alone==stim_opp_comp{i}(ii));
        ind = [ind; intersect(goodfit_dist_ind,resp_ind_alone{iii})];
    end
    subplot(2,2,i+2)
    scatter(avg_resp_stim(ind,stim_opp_comp{i}(1),1)+avg_resp_stim(ind,stim_opp_comp{i}(2),1),avg_resp_stim(ind,stim_opp(i),1))
    title('Opp')
    xlim([0 0.7])
    ylim([0 0.7])
    refline(1);
    xlabel('Left+Right')
    ylabel('Together')
end
suptitle([date ' ' mouse ' cells w/in' num2str(min_az_dist) '/' num2str(min_el_dist) ' Az/El deg'])
print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_scatterAloneVsTogether_dist.pdf']),'-dpdf','-bestfit')


figure;
sing_same = [];
sing_opp = [];
tog_same = [];
tog_opp = [];
for i = 1:length(stim_same)
    ind = [];
    for ii = 1:length(stim_same_comp{i})
        iii = find(stim_alone==stim_same_comp{i}(ii));
        ind = [ind; intersect(goodfit_dist_ind,resp_ind_alone{iii})];
    end
    subplot(1,2,1)
    scatter(avg_resp_stim(ind,stim_same_comp{i}(1),1)+avg_resp_stim(ind,stim_same_comp{i}(2),1),avg_resp_stim(ind,stim_same(i),1),'ok')
    sing_same = [sing_same; avg_resp_stim(ind,stim_same_comp{i}(1),1)+avg_resp_stim(ind,stim_same_comp{i}(2),1)]; 
    tog_same = [tog_same; avg_resp_stim(ind,stim_same(i),1)];
    hold on
    ind = [];
    for ii = 1:length(stim_opp_comp{i})
        iii = find(stim_alone==stim_opp_comp{i}(ii));
        ind = [ind; intersect(goodfit_dist_ind,resp_ind_alone{iii})];
    end
    subplot(1,2,2)
    scatter(avg_resp_stim(ind,stim_opp_comp{i}(1),1)+avg_resp_stim(ind,stim_opp_comp{i}(2),1),avg_resp_stim(ind,stim_opp(i),1),'ok')
    hold on
    sing_opp = [sing_opp; avg_resp_stim(ind,stim_opp_comp{i}(1),1)+avg_resp_stim(ind,stim_opp_comp{i}(2),1)];
    tog_opp = [tog_opp; avg_resp_stim(ind,stim_opp(i),1)];
end
subplot(1,2,1)
errorbar(mean(sing_same,1),mean(tog_same,1),std(tog_same,[],1)./sqrt(size(tog_same,1)),std(tog_same,[],1)./sqrt(size(tog_same,1)),std(sing_same,[],1)./sqrt(size(sing_same,1)),std(sing_same,[],1)./sqrt(size(sing_same,1)),'or')
title('Same')
xlim([0 2])
ylim([0 2])
refline(1);
xlabel('Left+Right')
ylabel('Together')
axis square
subplot(1,2,2)
errorbar(mean(sing_opp,1),mean(tog_opp,1),std(tog_opp,[],1)./sqrt(size(tog_opp,1)),std(tog_opp,[],1)./sqrt(size(tog_opp,1)),std(sing_opp,[],1)./sqrt(size(sing_opp,1)),std(sing_opp,[],1)./sqrt(size(sing_opp,1)),'or')
title('Opp')
xlim([0 2])
ylim([0 2])
refline(1);
xlabel('Left+Right')
ylabel('Together')
axis square
suptitle([date ' ' mouse ' cells w/in' num2str(min_az_dist) '/' num2str(min_el_dist) ' Az/El deg'])
print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_scatterAloneVsTogetherAll_dist.pdf']),'-dpdf','-bestfit')


figure;
for i = 1:2
    subplot(2,2,i)
    histogram(sameStim_SI{i}(intersect(goodfit_dist_ind,resp_ind_dir{i})),[-1:0.2:1])
    if intersect(resp_ind_dir{i},find(hsame_SI{i}))
        hold on
        histogram(sameStim_SI{i}(intersect(intersect(goodfit_dist_ind,resp_ind_dir{i}),find(hsame_SI{i}))),[-1:0.2:1])
    end
    title('Same')
    subplot(2,2,i+2)
     histogram(oppStim_SI{i}(intersect(goodfit_dist_ind,resp_ind_dir{i})),[-1:0.2:1])
    if intersect(intersect(goodfit_dist_ind,resp_ind_dir{i}),find(hopp_SI{i}))
        hold on
        histogram(oppStim_SI{i}(intersect(intersect(goodfit_dist_ind,resp_ind_dir{i}),find(hopp_SI{i}))),[-1:0.2:1])
    end
    title('Opp')
end
suptitle([date ' ' mouse ' cells w/in' num2str(min_az_dist) '/' num2str(min_el_dist) ' Az/El deg'])
print(fullfile(fn_out, datemouse, datemouserun, [datemouserun '_SIhists_dist.pdf']),'-dpdf','-bestfit')

sameStim_SI_all = [];
oppStim_SI_all= [];
for i = 1:2
    sameStim_SI_all = [sameStim_SI_all; sameStim_SI{i}(intersect(goodfit_dist_ind,resp_ind_dir{i}))];
    oppStim_SI_all = [oppStim_SI_all; oppStim_SI{i}(intersect(goodfit_dist_ind,resp_ind_dir{i}))];
end
figure;
subplot(2,2,1)
histogram(sameStim_SI_all,[-1:0.2:1])
xlabel('Selectivity Index')
title('Same')
subplot(2,2,2)
histogram(oppStim_SI_all,[-1:0.2:1])
xlabel('Selectivity Index')
title('Opp')
subplot(2,2,3)
cdfplot(sameStim_SI_all)
hold on
cdfplot(oppStim_SI_all)
legend('Same','Opp','location','southeast')
title('')
xlabel('Selectivity Index')

sameStim_SI_all = [];
oppStim_SI_all= [];
for i = 1:2
    sameStim_SI_all = [sameStim_SI_all; sameStim_SI{i}(resp_ind_dir{i})];
    oppStim_SI_all = [oppStim_SI_all; oppStim_SI{i}(resp_ind_dir{i})];
end
figure;
subplot(2,2,1)
histogram(sameStim_SI_all,[-1:0.2:1])
xlabel('Selectivity Index')
title('Same')
subplot(2,2,2)
histogram(oppStim_SI_all,[-1:0.2:1])
xlabel('Selectivity Index')
title('Opp')
subplot(2,2,3)
cdfplot(sameStim_SI_all)
hold on
cdfplot(oppStim_SI_all)
legend('Same','Opp','location','southeast')
title('')
xlabel('Selectivity Index')
[h p] = kstest2(sameStim_SI_all,oppStim_SI_all);