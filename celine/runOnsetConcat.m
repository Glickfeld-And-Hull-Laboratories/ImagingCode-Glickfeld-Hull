%read in the run_onset_analsis .mat files
ds = 'DART_V1_atropine_Celine';
pre = 2; post = 1;
rc = behavConstsDART;
eval(ds);
nd = 2;

expList=[66,57,53,49,44];
nExp =length(expList);
stillTimeList=[5];
nStillTime=length(stillTimeList)

tc_data=cell(nd,nStillTime);
mean_resp_data=cell(nd,nStillTime);
red_ind_concat=[];
nOnsets=zeros(2,nExp);

for iExp = 1:nExp
    day_id=expList(iExp)
    pre_day = expt(day_id).multiday_matchdays;

    mouse = expt(day_id).mouse;
    experimentFolder = expt(day_id).exptType;
    if expt(day_id).multiday_timesincedrug_hours > 0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(rc.achAnalysis, experimentFolder, mouse, ['multiday_' dart_str]);
    
    for iStillTime = 1:nStillTime
        
        thisStillTime=stillTimeList(iStillTime)
        fn_still = fullfile(fn_multi,[num2str(thisStillTime),'sec_noITI'])
        load(fullfile(fn_still,"run_onset_analysis.mat"))
        onsetSummaryTable=readtable(fullfile(fn_still,'onset_filtering_summary.csv'))
        %this is hardcoded to maintain the same pre/post order that the
        %rest of the data has
        nOnsets(:,iExp)=onsetSummaryTable.Count(strcmp(onsetSummaryTable.Reason, 'kept'));
       
    
        for id=1:nd
            if iExp==1
                tc_data{id,iStillTime}=data_dfof_runOnset_match{id};
                mean_resp_data{id,iStillTime}=mean_resp_runOnset_match{id};
            else
                tc_data{id,iStillTime}=[tc_data{id,iStillTime},data_dfof_runOnset_match{id}];
                mean_resp_data{id,iStillTime}=[mean_resp_data{id,iStillTime},mean_resp_runOnset_match{id}];
            end
        end


    end
red_ind_concat=[red_ind_concat,red_ind_match];
end
red_ind_concat=logical(red_ind_concat);
clear data_dfof_runOnset_match mean_resp_runOnset_match red_ind_match
%% visualize pre-DART running respone by stillness criterion, for each cell type
figure;subplot(2,1,1)
histogram(mean_resp_data{id,iStillTime}(~red_ind_concat))
xlim([-.2 .2])
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(~red_ind_concat)),' -HTP cells'])
vline(0.036)
subplot(2,1,2)
histogram(mean_resp_data{id,iStillTime}(red_ind_concat))
xlim([-.2 .2])
set(gca, 'TickDir', 'out');
grid off;
box off;
title([num2str(sum(red_ind_concat)),' +HTP cells'])
sgtitle(sprintf('%d mice with %s onsets', ...
    nExp, strjoin(arrayfun(@num2str, nOnsets(pre,:), 'UniformOutput', false), ', ')));
vline(0.036)
%% barchart of means and SEM forall cells
greenMeans = cell(1,2);
redMeans = cell(1,2);
greenSEM = cell(1,2);
redSEM = cell(1,2);

for id = 1:2
    for iStillTime = 1:nStillTime
        % Green (non-red) cells
        greenData = mean_resp_data{id,iStillTime}(~red_ind_concat);
        greenMeans{id}(iStillTime) = mean(greenData);
        greenSEM{id}(iStillTime) = std(greenData) / sqrt(length(greenData));
        
        % Red cells
        redData = mean_resp_data{id,iStillTime}(red_ind_concat);
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

% Create positions with a gap between groups
x_positions = [1:length(green_means), (length(green_means)+2):(2*length(green_means)+1)];

% Create the bar chart
figure;
h = bar(x_positions, all_means, 'grouped');
h.FaceColor = 'flat';
h.CData(1:length(green_means),:) = repmat([215, 217, 177]/255, length(green_means), 1); % Green color
h.CData(length(green_means)+1:end,:) = repmat([163, 0, 0]/255, length(red_means), 1);   % Red color

% Add error bars
hold on;
errorbar(x_positions, all_means, all_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Formatting
xlabel('Seconds of stillness required');
ylabel('Run onset response (first 2 seconds), pre-DART');

% Custom x-axis labels
xticks(x_positions);
xticklabels([string(stillTimeList), string(stillTimeList)]);

% Add group labels - positioned above the bars
text(mean(1:length(green_means)), max(all_means + all_stds) + 0.1*range(all_means), '-HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(mean((length(green_means)+2):(2*length(green_means)+1)), max(all_means + all_stds) + 0.1*range(all_means), '+HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'TickDir', 'out');
grid off;
box off;
hold off;
%% find cells that are facilitated or suppressed by running onset
% Analysis for NON-RED cells
cellData_nonRed = mean_resp_data{id,iStillTime}(~red_ind_concat);

% Calculate standard deviation (comparing against 0, not mean)
stdResp_nonRed = std(cellData_nonRed, 'omitnan');

% Define thresholds relative to zero
upperThresh_nonRed = .5*(stdResp_nonRed);     % +1 std above zero
lowerThresh_nonRed = -.5*(stdResp_nonRed);    % -1 std below zero

% Identify high and low responding cells
highResponders_nonRed = cellData_nonRed > upperThresh_nonRed;
lowResponders_nonRed = cellData_nonRed < lowerThresh_nonRed;

% Analysis for RED cells
cellData_red = mean_resp_data{id,iStillTime}(red_ind_concat);

% Calculate standard deviation (comparing against 0, not mean)
stdResp_red = std(cellData_red, 'omitnan');

% Define thresholds relative to zero
upperThresh_red = .5*(stdResp_red);     % +1 std above zero
lowerThresh_red = -.5*(stdResp_red);    % -1 std below zero

% Identify high and low responding cells
highResponders_red = cellData_red > upperThresh_red;
lowResponders_red = cellData_red < lowerThresh_red;

% Display results
fprintf('=== NON-RED CELLS ===\n');
fprintf('Std response: %.4f\n', stdResp_nonRed);
fprintf('Upper threshold (+1 std): %.4f\n', upperThresh_nonRed);
fprintf('Lower threshold (-1 std): %.4f\n', lowerThresh_nonRed);
fprintf('Total non-red cells: %d\n', sum(~red_ind_concat));
fprintf('High responders (>+1 std above 0): %d\n', sum(highResponders_nonRed));
fprintf('Low responders (<-1 std below 0): %d\n', sum(lowResponders_nonRed));

fprintf('\n=== RED CELLS ===\n');
fprintf('Std response: %.4f\n', stdResp_red);
fprintf('Upper threshold (+1 std): %.4f\n', upperThresh_red);
fprintf('Lower threshold (-1 std): %.4f\n', lowerThresh_red);
fprintf('Total red cells: %d\n', sum(red_ind_concat));
fprintf('High responders (>+1 std above 0): %d\n', sum(highResponders_red));
fprintf('Low responders (<-1 std below 0): %d\n', sum(lowResponders_red));

% Optional: Get indices if needed
nonRedIndices = find(~red_ind_concat);
redIndices = find(red_ind_concat);
highResponderIndices_nonRed = nonRedIndices(highResponders_nonRed);
lowResponderIndices_nonRed = nonRedIndices(lowResponders_nonRed);
highResponderIndices_red = redIndices(highResponders_red);
lowResponderIndices_red = redIndices(lowResponders_red);


%% barchart of means and SEM for high responders
greenMeans = cell(1,2);
redMeans = cell(1,2);
greenSEM = cell(1,2);
redSEM = cell(1,2);

for id = 1:2
    for iStillTime = 1:nStillTime
        % Green (non-red) cells
        greenData = mean_resp_data{id,iStillTime}(highResponderIndices_nonRed);
        greenMeans{id}(iStillTime) = mean(greenData);
        greenSEM{id}(iStillTime) = std(greenData) / sqrt(length(greenData));
        
        % Red cells
        redData = mean_resp_data{id,iStillTime}(highResponderIndices_red);
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

% Create positions with a gap between groups
x_positions = [1:length(green_means), (length(green_means)+2):(2*length(green_means)+1)];

% Create the bar chart
figure;
h = bar(x_positions, all_means, 'grouped');
h.FaceColor = 'flat';
h.CData(1:length(green_means),:) = repmat([215, 217, 177]/255, length(green_means), 1); % Green color
h.CData(length(green_means)+1:end,:) = repmat([163, 0, 0]/255, length(red_means), 1);   % Red color

% Add error bars
hold on;
errorbar(x_positions, all_means, all_stds, 'k', 'LineStyle', 'none', 'LineWidth', 1);

% Formatting
xlabel('Seconds of stillness required');
ylabel('Run onset response (first 2 seconds), pre-DART');

% Custom x-axis labels
xticks(x_positions);
xticklabels([string(stillTimeList), string(stillTimeList)]);

% Add group labels - positioned above the bars
text(mean(1:length(green_means)), max(all_means + all_stds) + 0.1*range(all_means), '-HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(mean((length(green_means)+2):(2*length(green_means)+1)), max(all_means + all_stds) + 0.1*range(all_means), '+HTP', ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold');

set(gca, 'TickDir', 'out');
grid off;
box off;
hold off;
%% plot timecourses

redMeansTC = cell(1,2);
greenMeansTC = cell(1,2);
greenSEM_TC = cell(1,2);
redSEM_TC = cell(1,2);

for id = 1:2
    for iStillTime = 1:nStillTime
        % Green (non-red) cells
        greenData = tc_data{id,iStillTime}(:,highResponderIndices_nonRed);
        greenMeansTC{id}(iStillTime,:) = mean(greenData,2);
        greenSEM_TC{id}(iStillTime,:) = std(greenData') / sqrt(length(highResponderIndices_nonRed));
        
        % Red cells
        redData = tc_data{id,iStillTime}(:,highResponderIndices_red);
        redMeansTC{id}(iStillTime,:) = mean(redData,2);
        redSEM_TC{id}(iStillTime,:) = std(redData') / sqrt(length(highResponderIndices_red));
        clear greenData redData
    end
end
timeVec=1:length(greenMeansTC{id}(iStillTime,:));
timeVec=(timeVec./frame_rate_double)-2;
figure;
shadedErrorBar(timeVec,greenMeansTC{pre}(1,:),greenSEM_TC{pre}(1,:),'--k')
hold on
shadedErrorBar(timeVec,redMeansTC{pre}(1,:),redSEM_TC{pre}(1,:),'k')

%% 

