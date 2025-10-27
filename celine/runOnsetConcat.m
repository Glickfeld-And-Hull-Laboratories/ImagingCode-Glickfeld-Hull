clear all; %close all; clc

prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
eval(instr);

ds = instructions.ds;
eval(ds);

rc = behavConstsDART;
experimentFolder = instructions.experimentFolder;

sess_list = instructions.runOnset_sess_list;
nSess = length(sess_list);

x = instructions.refDay;
switch x
    case '1'
        pre = 1;
        post = 2;
        fprintf('Baseline used as reference\n');
    case '2'
        pre = 2;
        post = 1;
        fprintf('Post-DART used as reference\n');
end
clear x 

% Create output directory
sess_title = num2str(sess_list(1));
for iSess = 2:nSess
    sess_title = [sess_title '_' num2str(sess_list(iSess))];
end
d = char(string(datetime('today')));

if nSess == 1
    if expt(sess_list(1)).multiday_timesincedrug_hours > 0
        dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fnout = fullfile(rc.analysis, experimentFolder, expt(sess_list(1)).mouse, ['multiday_' dart_str], d);
else
    fnout = fullfile(rc.analysis, experimentFolder, ['concat' sess_title], d);
end
mkdir(fnout); cd(fnout)

if isfield(instructions, 'stillTimeList')
    stillTimeList = instructions.stillTimeList;
else
    stillTimeList = [5];
    fprintf('Using default stillTimeList: [5]\n');
end
nStillTime = length(stillTimeList);

nd = 2;

tc_data = cell(nd, nStillTime);
mean_resp_data = cell(nd, nStillTime);
tc_data_bySess = cell(nd, nStillTime, nSess);
red_ind_concat = [];
red_ind_bySess = cell(1, nSess);
nOnsets = zeros(2, nSess);

for iExp = 1:nSess
    day_id = sess_list(iExp);
    fprintf('Processing session %d of %d (day_id %d)...\n', iExp, nSess, day_id);
    
    pre_day = expt(day_id).multiday_matchdays;
    
    mouse = expt(day_id).mouse;
    if expt(day_id).multiday_timesincedrug_hours > 0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.analysis, experimentFolder, mouse, ['multiday_' dart_str]);
    
    for iStillTime = 1:nStillTime
        thisStillTime = stillTimeList(iStillTime);
        fn_still = fullfile(fn_multi, [num2str(thisStillTime), 'sec_noITI']);
        
        if ~exist(fn_still, 'dir')
            warning('Directory not found: %s', fn_still);
            continue;
        end
        
        load(fullfile(fn_still, "run_onset_analysis.mat"));
        
        for id = 1:nd
            tc_data_bySess{id, iStillTime, iExp} = data_dfof_runOnset_match{id};
            if iExp == 1
                tc_data{id, iStillTime} = data_dfof_runOnset_match{id};
                mean_resp_data{id, iStillTime} = mean_resp_runOnset_match{id};
            else
                tc_data{id, iStillTime} = [tc_data{id, iStillTime}, data_dfof_runOnset_match{id}];
                mean_resp_data{id, iStillTime} = [mean_resp_data{id, iStillTime}, mean_resp_runOnset_match{id}];
            end
        end
        
        nOnsets(pre, iExp) = sum(strcmp(onsetITIStatus{pre}, 'ITI')) + sum(strcmp(onsetITIStatus{pre}, 'stimulus'));
        nOnsets(post, iExp) = sum(strcmp(onsetITIStatus{post}, 'ITI')) + sum(strcmp(onsetITIStatus{post}, 'stimulus'));
    end
    
    red_ind_concat = [red_ind_concat, red_ind_match];
    red_ind_bySess{iExp} = logical(red_ind_match);
    length(red_ind_match)
    length(red_ind_concat)

end

red_ind_concat = logical(red_ind_concat);
frame_rate_double = double(instructions.frame_rate);

clear data_dfof_runOnset_match mean_resp_runOnset_match red_ind_match

%% visualize pre-DART running response by stillness criterion, for each cell type
figure;
subplot(2, 1, 1);
histogram(mean_resp_data{pre, 1}(~red_ind_concat));
xlim([-.1 .5]);
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(~red_ind_concat)), ' -HTP cells']);
xline(nanmean(mean_resp_data{pre, 1}(~red_ind_concat)));

subplot(2, 1, 2);
histogram(mean_resp_data{pre, 1}(red_ind_concat));
xlim([-.1 .5]);
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(red_ind_concat)), ' +HTP cells']);
xline(nanmean(mean_resp_data{pre, 1}(red_ind_concat)));

sgtitle(sprintf('%d mice with %s onsets on the bsln day', ...
    nSess, strjoin(arrayfun(@num2str, nOnsets(pre, :), 'UniformOutput', false), ', ')));
print(fullfile(fnout, 'histogram_run_onset_response'), '-dpng', '-r300');

%% barchart of means and SEM for all cells
greenMeans = cell(1, 2);
redMeans = cell(1, 2);
greenSEM = cell(1, 2);
redSEM = cell(1, 2);
for id = 1:2
    for iStillTime = 1:nStillTime
        greenData = mean_resp_data{id, iStillTime}(~red_ind_concat);
        greenMeans{id}(iStillTime) = nanmean(greenData);
        greenSEM{id}(iStillTime) = nanstd(greenData) / sqrt(sum(~isnan(greenData)));
        redData = mean_resp_data{id, iStillTime}(red_ind_concat);
        redMeans{id}(iStillTime) = nanmean(redData);
        redSEM{id}(iStillTime) = nanstd(redData) / sqrt(sum(~isnan(redData)));
    end
end

figure;
green_data = [greenMeans{pre}; greenMeans{post}]';
red_data = [redMeans{pre}; redMeans{post}]';
all_data = [green_data; red_data];

h = bar(all_data);
h(1).FaceColor = [0.3, 0.3, 0.3];
h(2).FaceColor = [0, 114, 178]/255;
hold on;

green_sem_data = [greenSEM{1}; greenSEM{2}]';
red_sem_data = [redSEM{1}; redSEM{2}]';
all_sem = [green_sem_data; red_sem_data];

nGroups = size(all_data, 1);
nBars = size(all_data, 2);
x = nan(nBars, nGroups);
for i = 1:nBars
    x(i,:) = h(i).XEndPoints;
end
errorbar(x', all_data, all_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1);

ylim_vals = [0 max(all_data(:) + all_sem(:)) * 1.15];
ylim(ylim_vals);

xlabel('Seconds of stillness required');
ylabel('Run onset response');
xticks(1:nGroups);
xticklabels([string(stillTimeList), string(stillTimeList)]);
text(mean(1:nStillTime), ylim_vals(2)*0.95, '-HTP', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(mean((nStillTime+1):nGroups), ylim_vals(2)*0.95, '+HTP', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
legend({'pre-DART', 'post-DART'}, 'Location', 'best');
set(gca, 'TickDir', 'out');
grid off;
box off;
hold off;
print(fullfile(fnout, 'barchart_means_SEM'), '-dpng', '-r300');

%% plot timecourses
redMeansTC = cell(1, 2);
greenMeansTC = cell(1, 2);
greenSEM_TC = cell(1, 2);
redSEM_TC = cell(1, 2);

for id = 1:2
    for iStillTime = 1:nStillTime
        greenData = tc_data{id, iStillTime}(:, ~red_ind_concat);
        greenMeansTC{id}(iStillTime, :) = nanmean(greenData, 2);
        greenSEM_TC{id}(iStillTime, :) = nanstd(greenData, [], 2) / sqrt(sum(~red_ind_concat));
        
        redData = tc_data{id, iStillTime}(:, red_ind_concat);
        redMeansTC{id}(iStillTime, :) = nanmean(redData, 2);
        redSEM_TC{id}(iStillTime, :) = nanstd(redData, [], 2) / sqrt(sum(red_ind_concat));
        clear greenData redData
    end
end

timeVec = 1:length(greenMeansTC{id}(iStillTime, :));
timeVec = (timeVec ./ frame_rate_double) - 5;

figure;
subplot(2,1,1)
shadedErrorBar(timeVec, greenMeansTC{pre}(1, :), greenSEM_TC{pre}(1, :), '--k');
hold on;
shadedErrorBar(timeVec, greenMeansTC{post}(1, :), greenSEM_TC{post}(1, :), '--b');
set(gca, 'TickDir', 'out');
grid off;
box off;
xlabel('Time from run onset (s)');
ylabel('dF/F');
title("-HTP")
xline(0, '--', 'Color', [0.5 0.5 0.5]);

subplot(2,1,2)
shadedErrorBar(timeVec, redMeansTC{pre}(1, :), redSEM_TC{pre}(1, :), 'k');
hold on;
shadedErrorBar(timeVec, redMeansTC{post}(1, :), redSEM_TC{post}(1, :), 'b');
set(gca, 'TickDir', 'out');
grid off;
box off;
xlabel('Time from run onset (s)');
ylabel('dF/F');
title("+HTP")
xline(0, '--', 'Color', [0.5 0.5 0.5]);
print(fullfile(fnout, 'timecourse_plot'), '-dpng', '-r300');

%% plot pre vs post timecourses for each session - HTP- cells
figure;
nRows = ceil(nSess/2);
nCols = min(nSess, 2);
for iExp = 1:nSess
    subplot(nRows, nCols, iExp)
    
    greenDataPre = tc_data_bySess{pre, 1, iExp}(:, ~red_ind_bySess{iExp});
    greenMeanPre = nanmean(greenDataPre, 2);
    greenSEMPre = nanstd(greenDataPre, [], 2) / sqrt(size(greenDataPre, 2));
    
    greenDataPost = tc_data_bySess{post, 1, iExp}(:, ~red_ind_bySess{iExp});
    greenMeanPost = nanmean(greenDataPost, 2);
    greenSEMPost = nanstd(greenDataPost, [], 2) / sqrt(size(greenDataPost, 2));
    
    shadedErrorBar(timeVec, greenMeanPre, greenSEMPre, '--k');
    hold on;
    shadedErrorBar(timeVec, greenMeanPost, greenSEMPost, '--b');
    xline(0, '--', 'Color', [0.5 0.5 0.5]);
    
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    xlabel('Time from run onset (s)');
    ylabel('dF/F');
    title(sprintf('Session %d (pre: %d, post: %d onsets)', sess_list(iExp), nOnsets(pre, iExp), nOnsets(post, iExp)));

end
sgtitle('-HTP cells');
print(fullfile(fnout, 'timecourse_by_session_HTPneg'), '-dpng', '-r300');

% plot pre vs post timecourses for each session - HTP+ cells
figure;
for iExp = 1:nSess
    subplot(nRows, nCols, iExp)
    
    redDataPre = tc_data_bySess{pre, 1, iExp}(:, red_ind_bySess{iExp});
    redMeanPre = nanmean(redDataPre, 2);
    redSEMPre = nanstd(redDataPre, [], 2) / sqrt(size(redDataPre, 2));
    
    redDataPost = tc_data_bySess{post, 1, iExp}(:, red_ind_bySess{iExp});
    redMeanPost = nanmean(redDataPost, 2);
    redSEMPost = nanstd(redDataPost, [], 2) / sqrt(size(redDataPost, 2));
    
    shadedErrorBar(timeVec, redMeanPre, redSEMPre, 'k');
    hold on;
    shadedErrorBar(timeVec, redMeanPost, redSEMPost, 'b');
    xline(0, '--', 'Color', [0.5 0.5 0.5]);
    
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    xlabel('Time from run onset (s)');
    ylabel('dF/F');
    title(sprintf('Session %d (pre: %d, post: %d onsets)', sess_list(iExp), nOnsets(pre, iExp), nOnsets(post, iExp)));

end
sgtitle('+HTP cells');
print(fullfile(fnout, 'timecourse_by_session_HTPpos'), '-dpng', '-r300');