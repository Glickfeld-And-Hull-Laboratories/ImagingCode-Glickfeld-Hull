
clear all; clear global; close all; clc

% Configuration
ds = 'DART_V1_YM90K_Celine';
dataStructLabels = {'contrastxori'};
experimentFolder = 'SST_YM90K';
sess_list = [8 10 14 20 22];
targetCon = [.25 .5 1];
frame_rate = 15;
nd = 2;

% Setup directories and experiment data
rc = behavConstsDART;
eval(ds);

% Get reference session info
[pre, post] = getReferenceDays();

% Initialize output directory
fnout = createOutputDirectory(rc, experimentFolder, sess_list, expt);
cd(fnout);

% Main analysis
[concatData, sessionInfo, nSize, nCon] = concatenateSessionData(rc, experimentFolder, sess_list, expt, targetCon, nd, pre, post);
[cellIndices, cellCounts] = analyzeCellSelection(concatData, nd, sess_list, nSize, nCon);
saveResults(fnout, cellCounts, sessionInfo.mouseNames);



function [pre, post] = getReferenceDays()
prompt = 'Which session was used as reference for matching: 0- baseline, 1- post-DART';
x = input(prompt);
switch x
    case 0
        pre = 1;
        post = 2;
        disp("baseline used as reference");
    case 1
        pre = 2;
        post = 1;
        disp("post-DART used as reference");
end
end

function fnout = createOutputDirectory(rc, experimentFolder, sess_list, expt)
nSess = length(sess_list);
sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title, '_', string(sess_list(iSess)));
end

d = string(datetime('today'));
if nSess == 1
    if expt(sess_list(1)).multiday_timesincedrug_hours > 0
        dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fnout = fullfile(rc.achAnalysis, experimentFolder, expt(sess_list(1)).mouse, ['multiday_' dart_str], d);
else
    fnout = fullfile(rc.achAnalysis, experimentFolder, strcat('concat', sess_title), d);
end
mkdir(fnout);
end

function [concatData, sessionInfo, nSize, nCon] = concatenateSessionData(rc, experimentFolder, sess_list, expt, targetCon, nd, pre, post)
nSess = length(sess_list);
nCon = length(targetCon);

% Load first session to determine nSize
firstSess = sess_list(1);
tempSessionData = loadSessionData(rc, experimentFolder, expt, firstSess);
[~, ~] = processSessionData(tempSessionData, targetCon, nd);
tSize_match = celleqel2mat_padded(tempSessionData.input(1).tGratingDiameterDeg(1:tempSessionData.nTrials(1)));
sizes = unique(tSize_match);
nSize = length(sizes);

% Initialize concatenated data structures
concatData = initializeConcatStructures(nd, nCon, nSize);
sessionInfo = struct('mice', [], 'drug', {cell(1, nSess)}, 'mouseNames', []);

for iSess = 1:nSess
    thisSess = sess_list(iSess);
    mouse = expt(thisSess).mouse;
    sessionInfo.mice = [sessionInfo.mice; mouse];
    sessionInfo.drug{iSess} = expt(thisSess).drug;
    
    % Load session data
    sessionData = loadSessionData(rc, experimentFolder, expt, thisSess);
    
    % Process and concatenate
    [processedData, sharedCon] = processSessionData(sessionData, targetCon, nd);
    concatData = concatenateData(concatData, processedData, nd, sharedCon, nSize);
    
    fprintf('Session %d completed\n', iSess);
end

% Create mouse indices
sessionInfo.mouseInds = createMouseIndices(concatData.nKeep_concat);
sessionInfo.mouseNames = string(sessionInfo.mice);
end

function concatData = initializeConcatStructures(nd, nCon, nSize)
concatData = struct();

% Initialize cell arrays for each day
dayFields = {'tc_trial_avrg_stat_concat', 'tc_trial_avrg_loc_concat', 'conBySize_resp_stat_concat', ...
    'conBySize_resp_loc_concat', 'h_concat', 'data_resp_concat', 'resp_keep_concat', ...
    'resp_max_keep_concat', 'pref_responses_loc_concat', 'pref_responses_stat_concat', ...
    'pref_peak_stat_concat', 'pref_peak_loc_concat', 'RIx_concat', 'wheel_corr_concat', ...
    'meanF_concat', 'norm_dir_resp_stat_concat', 'norm_dir_resp_loc_concat', ...
    'pref_nonPref_stat_concat', 'pref_nonPref_loc_concat', 'pref_dir_concat', ...
    'noiseCorr_concat', 'sigCorr_concat', 'nonPref_trial_avrg_stat_concat', ...
    'nonPref_trial_avrg_loc_concat', 'data_dfof_runOnset_concat'};

for field = dayFields
    concatData.(field{1}) = cell(1, nd);
end

% Initialize other arrays
concatData.red_concat = [];
concatData.green_concat = [];
concatData.nKeep_concat = [];
concatData.dirs_concat = [];
concatData.cons_concat = [];
concatData.dfof_max_diff_concat = [];
concatData.red_fluor_concat = [];
concatData.green_fluor_concat = [];
concatData.pref_allTrials_stat_concat = cell(nCon, nSize, nd);
concatData.pref_allTrials_loc_concat = cell(nCon, nSize, nd);
concatData.pref_allTrials_largePupil_concat = cell(nCon, nSize, nd);
concatData.pref_allTrials_smallPupil_concat = cell(nCon, nSize, nd);
end

function sessionData = loadSessionData(rc, experimentFolder, expt, thisSess)
mouse = expt(thisSess).mouse;
if expt(thisSess).multiday_timesincedrug_hours > 0
    dart_str = [expt(thisSess).drug '_' num2str(expt(thisSess).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end

fn_multi = fullfile(rc.achAnalysis, experimentFolder, mouse, ['multiday_' dart_str]);

% Load required files
files = {'tc_keep.mat', 'resp_keep.mat', 'input.mat', 'locomotion.mat', 'HT_pyr_relationship.mat'};
sessionData = struct();

for i = 1:length(files)
    data = load(fullfile(fn_multi, files{i}));
    fields = fieldnames(data);
    for j = 1:length(fields)
        sessionData.(fields{j}) = data.(fields{j});
    end
end
end

function [processedData, sharedCon] = processSessionData(sessionData, targetCon, nd)
% Extract trial information
tCon_match = cell(1, nd);
tDir_match = cell(1, nd);
tSize_match = cell(1, nd);

for id = 1:nd
    tCon_match{id} = celleqel2mat_padded(sessionData.input(id).tGratingContrast(1:sessionData.nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(sessionData.input(id).tGratingDirectionDeg(1:sessionData.nTrials(id)));
    tSize_match{id} = celleqel2mat_padded(sessionData.input(id).tGratingDiameterDeg(1:sessionData.nTrials(id)));
end

dirs = unique(tDir_match{2}); % using post day
cons = unique(tCon_match{2});
sharedCon = find(ismember(cons, targetCon));
sizes = unique(tSize_match{1});

processedData = sessionData;
processedData.dirs = dirs;
processedData.cons = cons;
processedData.sizes = sizes;
processedData.sharedCon = sharedCon;
end

function concatData = concatenateData(concatData, processedData, nd, sharedCon, nSize)
% Concatenate basic arrays
concatData.dirs_concat = [concatData.dirs_concat, processedData.dirs];
concatData.cons_concat = [concatData.cons_concat, processedData.cons(sharedCon)];
concatData.red_concat = [concatData.red_concat, processedData.red_keep_logical];
concatData.green_concat = [concatData.green_concat, processedData.green_keep_logical];
concatData.nKeep_concat = [concatData.nKeep_concat, size(processedData.tc_trial_avrg_stat{2}, 2)];

% Concatenate day-specific data
for id = 1:nd
    concatData.tc_trial_avrg_stat_concat{id} = cat(2, concatData.tc_trial_avrg_stat_concat{id}, ...
        processedData.tc_trial_avrg_stat{id}(:,:,:,:));
    concatData.tc_trial_avrg_loc_concat{id} = cat(2, concatData.tc_trial_avrg_loc_concat{id}, ...
        processedData.tc_trial_avrg_loc{id}(:,:,sharedCon,:));
    concatData.nonPref_trial_avrg_stat_concat{id} = cat(2, concatData.nonPref_trial_avrg_stat_concat{id}, ...
        processedData.nonPref_trial_avrg_stat{id}(:,:,sharedCon,:));
    concatData.nonPref_trial_avrg_loc_concat{id} = cat(2, concatData.nonPref_trial_avrg_loc_concat{id}, ...
        processedData.nonPref_trial_avrg_loc{id}(:,:,sharedCon,:));
    
    % Concatenate response data
    responseFields = {'resp_keep_concat', 'pref_responses_loc_concat', 'pref_responses_stat_concat', ...
        'pref_peak_stat_concat', 'pref_peak_loc_concat', 'norm_dir_resp_stat_concat', ...
        'norm_dir_resp_loc_concat', 'h_concat', 'conBySize_resp_stat_concat', ...
        'conBySize_resp_loc_concat', 'data_resp_concat'};
    
    for field = responseFields
        concatData = concatenateResponseField(concatData, processedData, field{1}, id, sharedCon);
    end
    
    % Handle special concatenations
    concatData.resp_max_keep_concat{id} = cat(1, concatData.resp_max_keep_concat{id}, ...
        processedData.resp_max_keep{id}(:,sharedCon,:));
    concatData.RIx_concat{id} = cat(1, concatData.RIx_concat{id}, sum(processedData.RIx{id}));
    concatData.wheel_corr_concat{id} = cat(2, concatData.wheel_corr_concat{id}, processedData.wheel_corr{id});
    
    meanF = mean(processedData.fullTC_keep{id}, 1);
    concatData.meanF_concat{id} = cat(2, concatData.meanF_concat{id}, meanF);
    
    concatData.pref_dir_concat{id} = cat(2, concatData.pref_dir_concat{id}, processedData.pref_dir_keep{id});
    concatData.noiseCorr_concat{id} = cat(2, concatData.noiseCorr_concat{id}, processedData.noiseCorr{id});
    concatData.sigCorr_concat{id} = cat(2, concatData.sigCorr_concat{id}, processedData.sigCorr{id});
    
    % Concatenate trial data by condition
    for i = 1:length(sharedCon)
        iCon = sharedCon(i);
        for iSize = 1:length(processedData.sizes)
            concatData.pref_allTrials_stat_concat{i,iSize,id} = [concatData.pref_allTrials_stat_concat{i,iSize,id}, ...
                processedData.pref_allTrials_stat{iCon,iSize,id}];
            concatData.pref_allTrials_loc_concat{i,iSize,id} = [concatData.pref_allTrials_loc_concat{i,iSize,id}, ...
                processedData.pref_allTrials_loc{iCon,iSize,id}];
        end
    end
end

concatData.dfof_max_diff_concat = cat(1, concatData.dfof_max_diff_concat, ...
    processedData.dfof_max_diff(:,sharedCon,:));
end

function concatData = concatenateResponseField(concatData, processedData, fieldName, id, sharedCon)
baseFieldName = strrep(fieldName, '_concat', '');
if contains(fieldName, {'pref_responses_loc', 'pref_responses_stat', 'pref_peak_stat', 'pref_peak_loc'})
    concatData.(fieldName){id} = cat(1, concatData.(fieldName){id}, ...
        processedData.(baseFieldName){id}(:,sharedCon,:));
elseif contains(fieldName, {'conBySize_resp_stat', 'conBySize_resp_loc'})
    concatData.(fieldName){id} = cat(1, concatData.(fieldName){id}, ...
        processedData.([baseFieldName '_keep']){id});
else
    concatData.(fieldName){id} = cat(1, concatData.(fieldName){id}, ...
        processedData.(baseFieldName){id});
end
end

function mouseInds = createMouseIndices(nKeep_concat)
nSess = length(nKeep_concat);
mouseInds = cell(1, nSess);
start = 1;
for iMouse = 1:nSess
    mouseInds{iMouse} = start:(start-1)+nKeep_concat(iMouse);
    start = start + nKeep_concat(iMouse);
end
end

function [cellIndices, cellCounts] = analyzeCellSelection(concatData, nd, sess_list, nSize, nCon)
nSess = length(sess_list);
pre = 1; post = 2; % Assuming these values

% Find cells with complete data
[haveRunning_both, haveStat_both] = findCellsWithCompleteData(concatData, nd, nCon, nSize);
runningCells = intersect(haveStat_both, haveRunning_both);

% Find responsive cells
[respToSmall, respToLarge] = findResponsiveCells(concatData, nd, nSize);
includeCells = findIncludedCells(concatData, nd, nSize, haveRunning_both);

% Get cell type indices
red_ind_concat = find(concatData.red_concat);
green_ind_concat = find(concatData.green_concat);

runningGreen = intersect(runningCells, green_ind_concat);
runningRed = intersect(runningCells, red_ind_concat);

cellIndices = struct('runningCells', runningCells, 'runningGreen', runningGreen, ...
    'runningRed', runningRed, 'red_ind_concat', red_ind_concat, ...
    'green_ind_concat', green_ind_concat, 'respToSmall', respToSmall, ...
    'respToLarge', respToLarge, 'includeCells', includeCells);

% Count cells by mouse
cellCounts = countCellsByMouse(cellIndices, concatData, nSess);
end

function [haveRunning_both, haveStat_both] = findCellsWithCompleteData(concatData, nd, nCon, nSize)
pre = 1; post = 2;

haveRunning_pre = sum(squeeze(sum(~isnan(concatData.conBySize_resp_loc_concat{pre}), 2)), 2) == nSize*nCon;
haveRunning_post = sum(squeeze(sum(~isnan(concatData.conBySize_resp_loc_concat{post}), 2)), 2) == nSize*nCon;
haveRunning_both = find(haveRunning_pre .* haveRunning_post);

haveStat_pre = sum(squeeze(sum(~isnan(concatData.conBySize_resp_stat_concat{pre}), 2)), 2) == nSize*nCon;
haveStat_post = sum(squeeze(sum(~isnan(concatData.conBySize_resp_stat_concat{post}), 2)), 2) == nSize*nCon;
haveStat_both = find(haveStat_pre .* haveStat_post);
end

function [respToSmall, respToLarge] = findResponsiveCells(concatData, nd, nSize)
pre = 1; post = 2;

% Small size (size 1)
respToSizeBothDays = cell(1, nd);
for id = 1:nd
    respToSizeBothDays{id} = sum(squeeze(sum(concatData.h_concat{id}(:,:,:,1), 2)), 2);
end
respToSmall = logical(respToSizeBothDays{pre} + respToSizeBothDays{post});

% Large size (nSize)
respToSizeBothDays = cell(1, nd);
for id = 1:nd
    respToSizeBothDays{id} = sum(squeeze(sum(concatData.h_concat{id}(:,:,:,nSize), 2)), 2);
end
respToLarge = logical(respToSizeBothDays{pre} + respToSizeBothDays{post});
end

function includeCells = findIncludedCells(concatData, nd, nSize, haveRunning_both)
pre = 1;
responCriteria = cell(1, nd);
for id = 1:nd
    responseCheck = sum(squeeze(sum(concatData.h_concat{id}(:,:,:,nSize), 2)), 2);
    responCriteria{id} = find(logical(responseCheck));
end
includeCells = intersect(responCriteria{pre}, find(haveRunning_both));
end

function cellCounts = countCellsByMouse(cellIndices, concatData, nSess)
mouseInds = createMouseIndices(concatData.nKeep_concat);

cellCounts = struct();
cellCounts.red = nan(nSess, 2);
cellCounts.green = nan(nSess, 2);

for iMouse = 1:nSess
    cellCounts.red(iMouse, 1) = length(intersect(cellIndices.runningRed', mouseInds{iMouse}));
    cellCounts.red(iMouse, 2) = length(intersect(cellIndices.red_ind_concat', mouseInds{iMouse}));
    cellCounts.green(iMouse, 1) = length(intersect(cellIndices.runningGreen', mouseInds{iMouse}));
    cellCounts.green(iMouse, 2) = length(intersect(cellIndices.green_ind_concat', mouseInds{iMouse}));
end
end

function saveResults(fnout, cellCounts, mouseNames)
cellCountTable = table(cellCounts.red, 'RowNames', cellstr(mouseNames));
cellCountTableGreen = table(cellCounts.green, 'RowNames', cellstr(mouseNames));
writetable(cellCountTable, fullfile(fnout, 'cellCounts.csv'), 'WriteRowNames', true);
writetable(cellCountTableGreen, fullfile(fnout, 'cellCountsGreen.csv'), 'WriteRowNames', true);
end
%% plot timecourses

%first plotting all