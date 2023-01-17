%% clear all; establish parameters and experiment
clc; clear all; close all
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'RFMapping', 'heatmaps');

ds = 'RFMapping_15Hz_ExptList_SG';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

iexp = 3; 

%% verify mouse & date
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).rfFolder;
time = expt(iexp).rfTime;
nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002';

fprintf(['2P imaging eye analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load data

load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

%% make new indices throwing out "bad" pupil trials

max_dist = 2;
seed = rng;
nCells = size(resp_cell{end,end,end},1);
nTrials = size(data_dfof_tc,3);
trial_n = zeros(nPhas,nEl,nAz);
trialInd = cell(nPhas,nEl,nAz);
    
for ip = 1:nPhas
    ind_Phas = find(stimPhas == Phas(ip));
    for iE = 1:nEl
        ind_El = find(stimEl == El(iE));
        for iA = 1:nAz
            ind_Az = find(stimAz == Az(iA));
            ind = intersect(ind_Phas,intersect(ind_Az,ind_El));
            trialInd{ip,iE,iA} = ind;
            trial_n(ip,iE,iA) = length(ind);
        end
    end
end
    
eye_n_all = nan(1,nPhas,nEl,nAz);
resp_all = [];
stim_all = [];
resp_avg = cell(nPhas, nEl,nAz);
resp_avg_m = nan(nPhas, nEl, nAz);

for ip = 1:nPhas
    for iE = 1:nEl
            [memb ind_El] = ismember(trialInd{1,iE,1},find(centroid_dist<max_dist));
        for iA = 1:nAz
            [memb ind_Az] = ismember(trialInd{1,1,iA},find(centroid_dist<max_dist));
            [memb ind] = ismember(trialInd{ip,iE,iA},find(centroid_dist<max_dist));
            resp_all = [resp_all resp_cell{ip,iE,iA}(:,find(ind))];
            trialInd{ip,iE,iA} = ind;
            resp_avg{ip,iE,iA} = [resp_avg{ip,iE,iA} resp_cell{ip,iE,iA}(:,find(ind))];
        end 
    end
end


%%  make heatmaps

if exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_heatmapData.mat']),'file')
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_heatmapData.mat']));
else
    for ic = 1:nCells
        resp_avg_m =  cellfun(@(x) mean(x, 2), resp_avg, 'UniformOutput', false);
    end
        
    for ip = 1:nPhas
        if ip == 1
            cellOn = squeeze(resp_avg(ip,:,:)); %6x6 cell array with nCells by nTrials in each cell 
            cellOnr = reshape(cellOn, nEl*nAz, []); %reshape into 36x1
            [nrow, ncols] = cellfun(@size, cellOnr); %how many rows and columns are in each cell
            maxLength = max(ncols);

            cellOnMat = zeros(maxLength, nEl*nAz);
            p_resp = zeros(nCells, 1);  %anova outputs one pvalue per run (some group is significantly different from some other group)
            for i = 1:nCells       %pad the column of each cell with NaNs such that all cells are now the same size
                cellOnPad = cellfun(@(x)[(x(i,1:end)) nan(1,maxLength-size(x,2))],cellOnr,'UniformOutput', false);
                for ii = 1:nEl*nAz      %pull each cell and turn it into columns in a matrix (number of trials by number of stimulus locations)
                cellOnMat(:,ii) = [cellOnPad{ii}'];
                end    
                [p_resp(i,:) tbl] = anova1(cellOnMat);    %anova1 assumes columns are groups unless told otherwise
                fprintf([num2str(i) ' of ' num2str(nCells) '\n'])
                close all
            end
            indOn05 = find(p_resp<0.05);
            indOn01 = find(p_resp<0.01);
            indOn001 = find(p_resp<0.001);
        end
        if ip == 2
            cellOff = squeeze(resp_avg(ip,:,:)); %6x6 cell array with nCells by nTrials in each cell 
            cellOffr = reshape(cellOff, nEl*nAz, []); %reshape into 36x1
            [nrow, ncols] = cellfun(@size, cellOffr); %how many rows and columns are in each cell
            maxLength = max(ncols);

            cellOffMat = zeros(maxLength, nEl*nAz);
            p_resp = zeros(nCells, 1);  %anova outputs one pvalue per run (some group is significantly different from some other group)
            for i = 1:nCells       %pad the column of each cell with NaNs such that all cells are now the same size
                cellOffPad = cellfun(@(x)[(x(i,1:end)) nan(1,maxLength-size(x,2))],cellOffr,'UniformOutput', false);
                for ii = 1:nEl*nAz      %pull each cell and turn it into columns in a matrix (number of trials by number of stimulus locations)
                cellOffMat(:,ii) = [cellOffPad{ii}'];
                end    
                [p_resp(i,:) tbl] = anova1(cellOffMat);    %anova1 assumes columns are groups unless told otherwise
                fprintf([num2str(i) ' of ' num2str(nCells) '\n'])
                close all
            end
            indOff05 = find(p_resp<0.05);
            indOff01 = find(p_resp<0.01);
            indOff001 = find(p_resp<0.001);        
        end
    end
end

indIntA = intersect(indOn05, indOff05);

cellOn = squeeze(resp_avg_m(1,:,:));
cellOn = cat(1,cellOn);
cellOn = reshape(cell2mat(cellOn), nCells, []);
dfofInd = [];
for i = 1:length(cellOn)
    x = reshape(squeeze(cellOn(i,:,:)), nEl*nAz, []);
    if max(x) >= 0.05
        dfofInd = [dfofInd; i];
    else
    end
end
index = intersect(dfofInd, indIntA);     
cellOn = cellOn(index,:,:);
cellOnMax = max(cellOn, [], 2:3);

cellOff = squeeze(resp_avg_m(2,:,:));
cellOff = cat(1,cellOff);
cellOff = reshape(cell2mat(cellOff), nCells, []);
cellOff = cellOff(index,:,:);
cellOffMax = max(cellOff, [], 2:3);

cellMax = max(cellOnMax,cellOffMax);
% cellNorm = (cellOn./abs(cellOnMax)) - (cellOff./abs(cellOffMax));
cellNorm = (cellOn - cellOff)./abs(cellMax);
% cellNorm = cellOn - cellOff;


figure;
s = 1;
for i = 1:40
    subplot(5,9,s)
    imagesc(squeeze(cellNorm(i,:,:)), [-2 2]);
    colormap(redblue)
    title = ({['Cell #', num2str(i)]});
    movegui('center')
    s = s + 1;
end


save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_heatmapData.mat']), 'resp_avg', 'indOn05', 'indOff05', 'indOn01', 'indOff01', 'indOn001', 'indOff001');

%% 
close all; clc;
nonz = cellfun(@nnz, trialInd);
nonzOn = squeeze(nonz(1,:,:));
nonzOff = squeeze(nonz(2,:,:));

figure;
subplot(1,3,1)
    histogram(nonzOff)
    hold on
    histogram(nonzOn)
    xlabel('number of trial repeats')
    ylabel('number of locations (81 total)')
    h.Title = 'Frequency of trial repeats';
    legend('nOff', 'nOn')
subplot(1,3,2)
    heatmap(nonzOff)
%     title('nOff')
subplot(1,3,3)
    heatmap(nonzOn)
%     title('nOn')

