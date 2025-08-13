%% set paths and read in data
clear all; clear global; clc

prompt = 'Enter ds (e.g., DART_V1_YM90K_Celine): ';
ds = input(prompt, 's');
clear x prompt

dataStructLabels = {'contrastxori'};
rc =  behavConstsDART; %directories
eval(ds);
day_id = input('Enter reference day id ');
pre_day = expt(day_id).multiday_matchdays;
allDays = [day_id,pre_day];

nd=2; %hardcoding the number of days for now
mouse = expt(day_id).mouse;
experimentFolder = expt(day_id).exptType;
fnout = fullfile(rc.achAnalysis,mouse);

if expt(day_id).multiday_timesincedrug_hours>0
    dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
else
    dart_str = 'control';
end
fn_multi = fullfile(rc.achAnalysis,experimentFolder,mouse,['multiday_' dart_str]);

%check which day is pre and which is post, set variables accordingly
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baeline session, used as reference, is in the 1st position
                    post=2;
                    "baseline used as reference"
                case 1
                  pre=2;
                  post=1; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

% load the data
cd(fn_multi)
load(fullfile(fn_multi,'timecourses.mat'))
load(fullfile(fn_multi,'multiday_alignment.mat'))
load(fullfile(fn_multi,'input.mat'))
% Rename input to inputStructure to avoid conflict with MATLAB's input function
inputStructure = input;
clear input

% load pupil data - optional because some experiments don't have this
prompt = 'Include pupil data? Eye analysis must be done for this experiment. Y = yes, N = no'
    x=input(prompt)
    switch x
        case 'Y'
            
            pupil=cell(1,nd);
            for id=1:nd
                % Get the appropriate folder path for each day
                mouse = expt(allDays(id)).mouse;
                date = expt(allDays(id)).date;
                imgFolder = expt(allDays(id)).contrastxori_runs{1};
                dayPath = fullfile(rc.achAnalysis,experimentFolder, mouse, date, imgFolder);
                % Load pupil data from the correct day's folder
                pupil{id}=load(fullfile(dayPath,'pupil.mat'));
            end
        case 'N' 
            fprint('Not including pupil data')

    end


%% BEGIN DATA EXTRACTION
%% Convert data to trials
prompt = 'How would you like to process trials: 0- using photoFrameFinder, 1- using counterValcorrect_noPhotodiode? ';
x = input(prompt);
switch x
    case 0
        % photodiode 

        for id = 1:nd %cycle through days
            mouse = expt(allDays(id)).mouse;
            date = expt(allDays(id)).date;
            imgFolder = expt(allDays(id)).contrastxori_runs{1};
            imgMatFile = [imgFolder '_000_000.mat'];
            dataPath = fullfile(rc.achData, mouse, date, imgFolder);
            load(fullfile(dataPath,imgMatFile));
            [cStimOn stimOffs] = photoFrameFinder_Sanworks(info.frame);
        end
        
    case 1
        
    correctedInputStructure=NaN(1,nd);
        for id = 1:nd
               correctedInputStructure(:,id)=counterValCorrect_noPhotodiode(inputStructure(:,id));
        end
        
end
clear x prompt data_f_match cellstd 

        data_dfof_trial_mfatch = cell(1,nd); %make an empty array that is 1 by however many days there are (1X2 usually)

        fractTimeActive_match = cell(1,nd);
        cellstd_match = cell(1,nd);
            nTrials(id) = length(cStimOn);
            [nFrames nCells] = size(cellTCs_match{id});
            
            data_trial_match = nan(nOn+nOff,nTrials(id),nCells);
            
            for itrial = 1:nTrials(id)
              if ~isnan(cStimOn(itrial)) & (cStimOn(itrial)+nOn+nOff/2)<nFrames
                data_trial_match(:,itrial,:) = cellTCs_match{id}(cStimOn(itrial)-nOff/2:cStimOn(itrial)-1+nOn+nOff/2,:);
              end
            end


            fractTimeActive_match{id} = zeros(1,nCells);
           

            data_f_match = mean(data_trial_match(1: (nOff/2),:,:),1);
            data_dfof_trial_match{id} = bsxfun(@rdivide,bsxfun(@minus,data_trial_match,data_f_match),data_f_match);
            meansub_match = cellTCs_match{id}-nanmean(cellTCs_match{id},1);
            cellstd = nanstd(meansub_match,[],1);
            cellstd_match{id}=cellstd;
            for iCell = 1:nCells
                fractTimeActive_match{id}(:,iCell) = length(find(meansub_match(:,iCell)>3.*cellstd(1,iCell)))./nFrames;
            end
        end
        
        stimStart = nOff/2;
        stimEnd = stimStart+nOn;
        
        % Store original time courses in case needed later
        cellTCs_match_OG = cellTCs_match;
%% Get stimulus parameters
nOn = inputStructure(1).nScansOn;
nOff = inputStructure(1).nScansOff;

%tells the contrast, direction and orientation for each trial each day
tCon_match = cell(1,nd);
tDir_match = cell(1,nd);
tOri_match = cell(1,nd);
tSize_match = cell(1,nd);

%find the contrasts, directions and orientations for each day
%in case of instances where the number of trails actually collected was not
%consistent with the number mWorks thinks occured, only take 1:nTrials as
%dictated above based on the number of frames recorded
for id = 1:nd
    nTrials(id) = size(cellTCs_match{id},1)/(nOn+nOff); %to account for times 
%when there is a disruption before the full set of trials is collcted, I'm 
%determining the number of trials each day by how many frames of data I 
%have divided by the number of frames per trial
    tCon_match{id} = celleqel2mat_padded(inputStructure(id).tGratingContrast(1:nTrials(id)));
    tDir_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDirectionDeg(1:nTrials(id)));
    tOri_match{id} = tDir_match{id};
    tOri_match{id}(find(tDir_match{id}>=180)) = tDir_match{id}(find(tDir_match{id}>=180))-180;
    tSize_match{id} = celleqel2mat_padded(inputStructure(id).tGratingDiameterDeg(1:nTrials(id)));
end
oris = unique(tOri_match{1});
dirs = unique(tDir_match{1});
cons = unique(tCon_match{1});
sizes = unique(tSize_match{1});
nOri = length(oris);
nCon = length(cons);
nDir = length(dirs);
nSize = length(sizes);
%% Find cells that were responsive on either day ("keep") and identify red cells within that subset

[h_match, p_match, resp_match] = findResponsiveCells(data_dfof_trial_match, ...
    tCon_match, tDir_match, tSize_match, stimStart, stimEnd, nTrials, ...
    nCells, nDir, nCon, nSize, dirs, cons, sizes);

% User parameters
remove_outliers = input('Remove outliers? (1=yes, 0=no): ');
if remove_outliers
    std_threshold = input('Standard deviation threshold for outlier removal (e.g., 3): ');
end

% Find all cells responsive on at least one day (regardless of red/green)
resp_either_day = resp_match{1} | resp_match{2};
keep_cells_temp = find(resp_either_day);

% Count initial responsive cells by type
initial_red_resp = sum(red_ind_match(keep_cells_temp));
initial_green_resp = sum(~red_ind_match(keep_cells_temp));
initial_total_resp = length(keep_cells_temp);

% Total cells by type (for reference)
total_red = sum(red_ind_match);
total_green = sum(~red_ind_match);
total_cells = length(red_ind_match);

% Optional outlier removal
if remove_outliers
    outliers_all = [];
    for id = 1:nd
        data_temp = data_resp_match{id}(keep_cells_temp,:,:,:,1);
        resp_max_temp = max(squeeze(max(data_temp(:,:,:,1),[],3)),[],2);
        thresh = nanmean(resp_max_temp) + std_threshold*std(resp_max_temp);
        outliers_all = union(outliers_all, keep_cells_temp(resp_max_temp > thresh));
    end
    
    fprintf('Removing %d outlier cells (>%.1f SD from mean)\n', length(outliers_all), std_threshold);
    keep_cells = setdiff(keep_cells_temp, outliers_all);
else
    fprintf('Skipping outlier removal\n');
    keep_cells = keep_cells_temp;
    outliers_all = [];
end

% Count final keep cells by type
final_red_keep = sum(red_ind_match(keep_cells_final));
final_green_keep = sum(~red_ind_match(keep_cells));
final_total_keep = length(keep_cells);

% Within the keep cells, identify which ones are red
red_cells_keep = red_ind_match(keep_cells_final);

% Create summary table
cell_summary = table(...
    {'Total'; 'Red (HTP+)'; 'Green (HTP-)'}, ...
    [total_cells; total_red; total_green], ...
    [initial_total_resp; initial_red_resp; initial_green_resp], ...
    [final_total_keep; final_red_keep; final_green_keep], ...
    [length(outliers_all); sum(red_ind_match(outliers_all)); sum(~red_ind_match(outliers_all))], ...
    'VariableNames', {'CellType', 'Total_Cells', 'Responsive_Cells', 'Final_Kept', 'Outliers_Removed'});

% Add percentage columns
cell_summary.Responsive_Percent = round(100 * cell_summary.Responsive_Cells ./ cell_summary.Total_Cells, 1);
cell_summary.Kept_Percent = round(100 * cell_summary.Final_Kept ./ cell_summary.Responsive_Cells, 1);

% Display summary
fprintf('\n--- Cell Filtering Summary ---\n');
disp(cell_summary);

% Save summary table
save('cell_filtering_summary.mat', 'cell_summary');
writetable(cell_summary, 'cell_filtering_summary.csv');

fprintf('\nSummary saved to: cell_filtering_summary.mat and cell_filtering_summary.csv\n');

% Clean up temporary variables
clear initial_red_resp initial_green_resp initial_total_resp
clear final_red_keep final_green_keep final_total_keep
clear total_red total_green total_cells outliers_all

%% For keep cells only, find preferred direction
%for the keep cells subset, find the preferred direction for each cell,
%then average over all trials at that direction for all sizes and contrasts
%to find each cell's trial-averaged timecourse ad each cell's
%trial-averaged mean response (ie, the mean over the respons window, so
%dropping the time dimension) for each size and contrast.
%% Get wheel speed, running vs. stationary trials
%% Get running onsets


%% Find large vs small pupil trials - this should also be optional

%%

%% Fluoresence level








