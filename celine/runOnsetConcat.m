clear all; close all; clc

prompt = 'Enter name of instructions file: ';
instr = input(prompt, 's');
clear prompt
eval(instr);

ds = instructions.ds;
eval(ds);

rc = behavConstsDART;
experimentFolder = instructions.experimentFolder;

sess_list = instructions.sess_list;
nExp = length(sess_list);

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
clear x instr

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
red_ind_concat = [];
nOnsets = zeros(2, nExp);

for iExp = 1:nExp
    day_id = sess_list(iExp);
    fprintf('Processing session %d of %d (day_id %d)...\n', iExp, nExp, day_id);
    
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
end

red_ind_concat = logical(red_ind_concat);
frame_rate_double = double(instructions.frame_rate);

clear data_dfof_runOnset_match mean_resp_runOnset_match red_ind_match

%% visualize pre-DART running response by stillness criterion, for each cell type
figure;
subplot(2, 1, 1);
histogram(mean_resp_data{pre, 1}(~red_ind_concat));
xlim([-.2 .2]);
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(~red_ind_concat)), ' -HTP cells']);
xline(0.036);

subplot(2, 1, 2);
histogram(mean_resp_data{pre, 1}(red_ind_concat));
xlim([-.2 .2]);
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(red_ind_concat)), ' +HTP cells']);
xline(0.036);

sgtitle(sprintf('%d mice with %s onsets', ...
    nExp, strjoin(arrayfun(@num2str, nOnsets(pre, :), 'UniformOutput', false), ', ')));

%% barchart of means and SEM for all cells
greenMeans = cell(1, 2);
redMeans = cell(1, 2);
greenSEM = cell(1, 2);
redSEM = cell(1, 2);

for id = 1:2
    for iStillTime = 1:nStillTime
        greenData = mean_resp_data{id, iStillTime}(~red_ind_concat);
        greenMeans{id}(iStillTime) = mean(greenData);
        greenSEM{id}(iStillTime) = std(greenData) / sqrt(length(greenData));
        
        redData = mean_resp_data{id, iStillTime}(red_ind_concat);
        redMeans{id}(iStillTime) = mean(redData);
        redSEM{id}(iStillTime) = std(redData) / sqrt(length(redData));
    end
end

green_means = greenMeans{pre};
red_means = redMeans{pre};
green_sem = greenSEM{pre};
red_sem = redSEM{pre};

all_means = [green_means, red_means];
all_stds = [green_sem, red_sem];

x_positions = [1:length(green_means), (length(green_means)+2):(2*length(green_means)+1)];

figure;
h = bar(x_positions, all_means, 'grouped');
h.FaceColor = 'flat';
h.CData(1:length(green_means), :) = repmat([215, 217, 177]/255, length(green_means), 1);
h.CData(length(green_means)+1:end, :) = repmat([163, 0, 0]/255, length(red_means), 1);

hold on;
errorbar(x_positions, all_means, all_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1);

xlabel('Seconds of stillness required');
ylabel('Run onset response (first 2 seconds), pre-DART');

xticks(x_positions);
xticklabels([string(stillTimeList), string(stillTimeList)]);

text(mean(1:length(green_means)), max(all_means + all_stds) + 0.1*range(all_means), '-HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(mean((length(green_means)+2):(2*length(green_means)+1)), max(all_means + all_stds) + 0.1*range(all_means), '+HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'TickDir', 'out');
grid off;
box off;
hold off;

%% find cells that are facilitated or suppressed by running onset
cellData_nonRed = mean_resp_data{pre, 1}(~red_ind_concat);
stdResp_nonRed = std(cellData_nonRed, 'omitnan');
upperThresh_nonRed = .5 * stdResp_nonRed;
lowerThresh_nonRed = -.5 * stdResp_nonRed;

highResponders_nonRed = cellData_nonRed > upperThresh_nonRed;
lowResponders_nonRed = cellData_nonRed < lowerThresh_nonRed;

cellData_red = mean_resp_data{pre, 1}(red_ind_concat);
stdResp_red = std(cellData_red, 'omitnan');
upperThresh_red = .5 * stdResp_red;
lowerThresh_red = -.5 * stdResp_red;

highResponders_red = cellData_red > upperThresh_red;
lowResponders_red = cellData_red < lowerThresh_red;

fprintf('=== NON-RED CELLS ===\n');
fprintf('Std response: %.4f\n', stdResp_nonRed);
fprintf('Upper threshold (+0.5 std): %.4f\n', upperThresh_nonRed);
fprintf('Lower threshold (-0.5 std): %.4f\n', lowerThresh_nonRed);
fprintf('Total non-red cells: %d\n', sum(~red_ind_concat));
fprintf('High responders (>+0.5 std above 0): %d\n', sum(highResponders_nonRed));
fprintf('Low responders (<-0.5 std below 0): %d\n', sum(lowResponders_nonRed));

fprintf('\n=== RED CELLS ===\n');
fprintf('Std response: %.4f\n', stdResp_red);
fprintf('Upper threshold (+0.5 std): %.4f\n', upperThresh_red);
fprintf('Lower threshold (-0.5 std): %.4f\n', lowerThresh_red);
fprintf('Total red cells: %d\n', sum(red_ind_concat));
fprintf('High responders (>+0.5 std above 0): %d\n', sum(highResponders_red));
fprintf('Low responders (<-0.5 std below 0): %d\n', sum(lowResponders_red));

nonRedIndices = find(~red_ind_concat);
redIndices = find(red_ind_concat);
highResponderIndices_nonRed = nonRedIndices(highResponders_nonRed);
lowResponderIndices_nonRed = nonRedIndices(lowResponders_nonRed);
highResponderIndices_red = redIndices(highResponders_red);
lowResponderIndices_red = redIndices(lowResponders_red);

%% barchart of means and SEM for high responders
greenMeans = cell(1, 2);
redMeans = cell(1, 2);
greenSEM = cell(1, 2);
redSEM = cell(1, 2);

for id = 1:2
    for iStillTime = 1:nStillTime
        greenData = mean_resp_data{id, iStillTime}(highResponderIndices_nonRed);
        greenMeans{id}(iStillTime) = mean(greenData);
        greenSEM{id}(iStillTime) = std(greenData) / sqrt(length(greenData));
        
        redData = mean_resp_data{id, iStillTime}(highResponderIndices_red);
        redMeans{id}(iStillTime) = mean(redData);
        redSEM{id}(iStillTime) = std(redData) / sqrt(length(redData));
    end
end

green_means = greenMeans{pre};
red_means = redMeans{pre};
green_sem = greenSEM{pre};
red_sem = redSEM{pre};

all_means = [green_means, red_means];
all_stds = [green_sem, red_sem];

x_positions = [1:length(green_means), (length(green_means)+2):(2*length(green_means)+1)];

figure;
h = bar(x_positions, all_means, 'grouped');
h.FaceColor = 'flat';
h.CData(1:length(green_means), :) = repmat([215, 217, 177]/255, length(green_means), 1);
h.CData(length(green_means)+1:end, :) = repmat([163, 0, 0]/255, length(red_means), 1);

hold on;
errorbar(x_positions, all_means, all_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1);

xlabel('Seconds of stillness required');
ylabel('Run onset response (first 2 seconds), pre-DART');

xticks(x_positions);
xticklabels([string(stillTimeList), string(stillTimeList)]);

text(mean(1:length(green_means)), max(all_means + all_stds) + 0.1*range(all_means), '-HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(mean((length(green_means)+2):(2*length(green_means)+1)), max(all_means + all_stds) + 0.1*range(all_means), '+HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'TickDir', 'out');
grid off;
box off;
hold off;

%% plot timecourses
redMeansTC = cell(1, 2);
greenMeansTC = cell(1, 2);
greenSEM_TC = cell(1, 2);
redSEM_TC = cell(1, 2);

for id = 1:2
    for iStillTime = 1:nStillTime
        greenData = tc_data{id, iStillTime}(:, ~red_ind_concat);
        greenMeansTC{id}(iStillTime, :) = mean(greenData, 2);
        greenSEM_TC{id}(iStillTime, :) = std(greenData, [], 2) / sqrt(sum(~red_ind_concat));
        
        redData = tc_data{id, iStillTime}(:, red_ind_concat);
        redMeansTC{id}(iStillTime, :) = mean(redData, 2);
        redSEM_TC{id}(iStillTime, :) = std(redData, [], 2) / sqrt(sum(red_ind_concat));
        clear greenData redData
    end
end

timeVec = 1:length(greenMeansTC{id}(iStillTime, :));
timeVec = (timeVec ./ frame_rate_double) - 5;  % Changed from -2 to -5

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
%legend({'-HTP', '+HTP'});
xline(0, '--', 'Color', [0.5 0.5 0.5]);  % Optional: mark the onset time

subplot(2,1,2)
shadedErrorBar(timeVec, redMeansTC{pre}(1, :), redSEM_TC{pre}(1, :), 'k');
hold on;
shadedErrorBar(timeVec, redMeansTC{post}(1, :), redSEM_TC{post}(1, :), 'b');
set(gca, 'TickDir', 'out');
grid off;
box off;
xlabel('Time from run onset (s)');
ylabel('dF/F');
title("-HTP")
%legend({'-HTP', '+HTP'});
xline(0, '--', 'Color', [0.5 0.5 0.5]);  % Optional: mark the onset time