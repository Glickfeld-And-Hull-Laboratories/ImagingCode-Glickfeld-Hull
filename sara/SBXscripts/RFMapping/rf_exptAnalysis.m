clc; clear all; close all;
ds = 'RFMapping_15Hz_ExptList_SG';

rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);

frame_rate = 15;

iexp  = 3;
%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).rfFolder;
time = expt(iexp).rfTime;
nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002';

base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];

fprintf([mouse ' ' date '\n'])

%% Test stim analysis
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
%data_tc, np_tc, npSub_tc, nCells, sz
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
%stim_list, listSqs, nOn, nOff, nTrials, stimEl, stimAz, stimPhas, El, Az, Phas, nEl, nAz, nPhas, nStim, frame_rate, nTrials
load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%% Find significant responses to stimuli
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'; 
summaryDir = fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);

fprintf(['Expt analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\n'])

data_resp = nan(nOff+nOn,nCells,nTrials);
data_f = nan(1,nCells,nTrials);

cStimOff = celleqel2mat_padded(input.counter);
nStimCon = 1:nEl*nAz;


for itrial = 1:nTrials
    if cStimOff(itrial) < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOff(itrial)-(nOff + nOn - 1):cStimOff(itrial),:);
    end
end


data_f = mean(data_resp(1:nOff,:,:),1);
data_dfof_tc = (data_resp-data_f)./data_f;


base_cell = cell(nPhas, nEl, nAz); 
resp_cell= cell(nPhas, nEl, nAz);
trialInd = cell(nPhas, nEl, nAz);
trialsperstim = zeros(nPhas, nEl, nAz);
base_win = 1:nOff;
resp_win = nOff:(nOn+nOff);

data_dfof_con_tc_avg = nan(nOn+nOff, nCells, nPhas, nEl, nAz);

resp = zeros(nCells, nEl*nAz);
h_resp = zeros(nCells, nPhas, nEl, nAz);
p_resp = zeros(nCells, nEl, nAz);
indOn = [];
indOff = [];
indResp = [];
indInt = [];
hOn = [];
hOff = [];

start = 1;
for ip = 1:nPhas
    ind_p = find(stimPhas == Phas(ip));
    for iE = 1:nEl
        ind_El = find(stimEl == El(iE));
        for iA = 1:nAz
            ind_Az = find(stimAz == Az(iA));
            ind = intersect(ind_p,intersect(ind_Az,ind_El));
            trialsperstim(ip,iE,iA) = length(ind);
            stim_list(start,:) = [Phas(ip) El(iE) Az(iA)];
            resp_cell{ip,iE,iA} = squeeze(mean(data_dfof_tc(resp_win,:,ind),1));
            base_cell{ip,iE,iA} = squeeze(mean(data_dfof_tc(base_win,:,ind),1));
            resp_cell_m{ip,iE,iA} = nanmean(squeeze(mean(data_dfof_tc(resp_win,:,ind),1)),2);
            data_dfof_con_tc_avg(:,:,ip,iE,iA,1) = squeeze(nanmean(data_dfof_tc(:,:,ind),3));
            data_dfof_con_tc_avg(:,:,ip,iE,iA,2) = squeeze(nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind)));
            [h_resp(:,ip,iE,iA) p_resp] =  ttest2(resp_cell{ip,iE,iA},base_cell{ip,iE,iA},'dim',2,'tail','right');
            trialInd{ip,iE,iA} = ind;
            start = start+1;
        end
    end
    if ip == 1
        hOn = h_resp(:,ip,:,:);
    end
    if ip == 2
        hOff = h_resp(:,ip,:,:);
    end
    indOn = find(sum(hOn,[3 4]));
    indOff = find(sum(hOff,[3 4]));
    indResp = unique([indOn; indOff]);
    indInt = intersect(indOn,indOff);
end

%% looking at RFs of cells that were significantly responsive (ttest) to at least one location on either subfield n=270

cellOn = squeeze(resp_cell_m(1,:,:));
cellOn = cat(3,cellOn{:});
cellOn = reshape(cellOn, nCells, nEl, nAz);
cellOn = cellOn(indResp,:,:);

cellOff = squeeze(resp_cell_m(2,:,:));
cellOff = cat(3,cellOff{:});
cellOff = reshape(cellOff, nCells, nEl, nAz);
cellOff = cellOff(indResp,:,:);

a = squeeze(cellOn(i,:,:));

figure;
s = 1;
for i = 1:3
    subplot(1,6,s)
    h = heatmap(squeeze(cellOn(i,:,:)), 'ColorMap', flipud(autumn), 'CellLabelColor', 'none');
    h.Title = 'On';
    h.GridVisible ='off';
    

    subplot(1,6,s+1)
    h = heatmap(squeeze(cellOff(i,:,:)), 'ColorMap', flipud(summer), 'CellLabelColor', 'none');
    h.Title = 'Off';
    h.GridVisible = 'off';

    movegui('center')
    s=s+2;
end 

    imagesc(squeeze(cellOn(15,:,:)))


%% looking cells that were significantly responsive (ttest) to at least one location on both subfields n=151



cellOn = squeeze(resp_cell_m(1,:,:));
cellOn = cat(3,cellOn{:});
cellOn = reshape(cellOn, nCells, nEl, nAz);
cellOn = cellOn(indInt,:,:);

cellOff = squeeze(resp_cell_m(2,:,:));
cellOff = cat(3,cellOff{:});
cellOff = reshape(cellOff, nCells, nEl, nAz);
cellOff = cellOff(indInt,:,:);


figure;
s = 1;
for i = 150
    subplot(1,2,s)
    h = heatmap(squeeze(cellOn(i,:,:)), 'ColorMap', parula);
    h.Title = 'On';
    
    subplot(1,2,s+1)
    h = heatmap(squeeze(cellOff(i,:,:)), 'ColorMap', parula);
    h.Title = 'Off';
    movegui('center')
    s=s+2;
end
    print(fullfile(summaryDir, ['_heatmaps.pdf']),'-dpdf', '-fillpage')


figure;
s = 1;
for i = 56
    subplot(1,2,s)
    a = pcolor(squeeze(cellOn(i,:,:)));
    
    subplot(1,2,s+1)
    a = imagesc(squeeze(cellOn(i,:,:)));
%     a.FaceColor = 'interp';
    s=s+2;
end 


%% testing significance with an ANOVA

if exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']),'file')
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']));
else
    for ip = 1:nPhas
        if ip == 1
            cellOn = squeeze(resp_cell(ip,:,:)); %6x6 cell array with nCells by nTrials in each cell 
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
            cellOff = squeeze(resp_cell(ip,:,:)); %6x6 cell array with nCells by nTrials in each cell 
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

indIntA = intersect(indOn01, indOff01);

cellOn = squeeze(resp_cell_m(1,:,:));
cellOn = cat(3,cellOn{:});
cellOn = reshape(cellOn, nCells, nEl, nAz);
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

cellOff = squeeze(resp_cell_m(2,:,:));
cellOff = cat(3,cellOff{:});
cellOff = reshape(cellOff, nCells, nEl, nAz);
cellOff = cellOff(index,:,:);
cellOffMax = max(cellOff, [], 2:3);

cellMax = max(cellOnMax,cellOffMax);
% cellNorm = (cellOn./abs(cellOnMax)) - (cellOff./abs(cellOffMax));
cellNorm = (cellOn - cellOff)./abs(cellOffMax);
% cellNorm = cellOn - cellOff;


figure;
s = 1;
for i = 1:45
    subplot(5,9,s)
%     myfilter = fspecial('gaussian',[20 20], 0.5);
%     cellNormGauss = imfilter(squeeze(cellNorm(i,:,:)), myfilter);
%     imagesc(cellNormGauss);
    imagesc(squeeze(cellNorm(i,:,:)), [-1.5 1.5]);
    colormap(redblue)
    title = ({['Cell #', num2str(i)]});
    movegui('center')
    s = s + 1;
end
    print(fullfile(summaryDir, ['_ANOVA_heatmaps.pdf']),'-dpdf', '-fillpage')




%% save variables

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']), 'data_dfof_con_tc_avg', 'resp_cell','base_cell', 'data_dfof_tc', 'indResp', 'indInt',  'indOn05', 'indOff05', 'indOn01', 'indOff01', 'indOn001', 'indOff001', 'frame_rate');
save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'resp_win', 'base_win', 'cStimOff', 'nStimCon', 'trialsperstim', 'stim_list', 'El', 'Az', 'Phas', 'nEl', 'nAz', 'nPhas', 'stimEl', 'stimAz', 'stimPhas');


%% looking into mworks randomization limit

[nrow, ncols] = cellfun(@size, trialInd);
ncolsOn = squeeze(ncols(1,:,:));
ncolsOff = squeeze(ncols(2,:,:));

figure;
subplot(1,3,1)
    histogram(ncolsOff)
    hold on
    histogram(ncolsOn)
    xlabel('number of trial repeats')
    ylabel('number of locations (81 total)')
    title('Frequency of trial repeats')
    legend('nOff', 'nOn')
subplot(1,3,2)
    heatmap(ncolsOff)
    title('nOff')
subplot(1,3,3)
    heatmap(ncolsOn)
    title('nOn')


