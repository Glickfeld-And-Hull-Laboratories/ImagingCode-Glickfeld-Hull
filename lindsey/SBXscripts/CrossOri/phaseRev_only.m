close all;clear all;clc;

ds = 'CrossOriRandDir_ExptList';
eval(ds)
nexp = length(expt);
rc = behavConstsAV;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';


%find new expt
for iexp = 1:nexp
    if ~isempty(expt(iexp).prFolder)
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).prFolder;
        time = expt(iexp).prTime;
        nrun = length(ImgFolder);
        frameRateHz = params.frameRate;
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        outDir = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str '_prOnly']);
        if ~exist(outDir)
            fprintf(['Found new experiment:\nMouse: ' mouse '\nDate: ' date '\n'])
            break
        end
    end
    if iexp == nexp
        fprintf('No new experiments')
    end
end

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    %CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);

    temp(irun) = input;
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = ntrials*(nOn+nOff);


    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,nframes);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    fprintf('Complete')
end
input = concatenateDataBlocks(temp);
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp
clear temp

toc
        
%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); colormap gray; clim([0 3000]); end
movegui('center')
data_avg = mean(data(:,:,20001:20500),3);
%% Register data
[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile(outDir))
save(fullfile(outDir, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(outDir, [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data out

%%
nOn = double(input.nScansOn);
nOff = double(input.nScansOff);
phaseCyc = double(input.nScansPhaseCyc);
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
sz = size(data_reg);

data_resp = reshape(data_reg, [sz(1) sz(2) nOn+nOff nTrials]);
data_f = mean(data_resp(:,:,nOff-15:nOff,:),3);
data_resp_dfof = (data_resp-data_f)./data_f;

clear data_resp data_f

dirs = unique(dir_mat);
nDir = length(dirs);

tStimNum = celleqel2mat_padded(input.tStimulusNumber);
nPhase = input.gratingStartingPhaseStepN;
phase_mat = zeros(1,ntrials);
for itrial = 1:ntrials
    if tStimNum(itrial) < (nDir.*nPhase)
        temp = tStimNum(itrial);
    else
        temp = mod(tStimNum(itrial),nDir.*nPhase);
    end
    if temp < nPhase
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*temp);
    else
        phase_mat(itrial) = input.gratingStartingPhaseDeg + (input.gratingStartingPhaseStepDeg.*(mod(temp,input.gratingStartingPhaseStepN)));
    end
end
phases = unique(phase_mat);

nStim = nPhase.*nDir*2;
[n n2] = subplotn(nStim);
data_dfof_stim = nan(sz(1),sz(2),nPhase,nDir,2);
start = 1;
for iPhase = 1:nPhase
    ind_p = find(phase_mat == phases(iPhase));
    for iDir = 1:nDir
        ind_d = find(dir_mat == dirs(iDir));
        ind = intersect(ind_p,ind_d);
        data_dfof_stim(:,:,iPhase,iDir,1) = mean(mean(data_resp_dfof(:,:,nOff+1:nOff+frameRateHz,ind),4),3);
        data_dfof_stim(:,:,iPhase,iDir,2) = mean(mean(data_resp_dfof(:,:,nOff+1+phaseCyc:nOff+phaseCyc+frameRateHz,ind),4),3);
        subplot(n,n2,start)
        imagesc(data_dfof_stim(:,:,iPhase,iDir,1))
        subplot(n,n2,start+1)
        imagesc(data_dfof_stim(:,:,iPhase,iDir,2))
        start = start+2;
    end
end
    
data_dfof_stim_all = reshape(data_dfof_stim, [sz(1) sz(2) nPhase.*nDir.*2]);
data_dfof_max = max(data_dfof_stim_all,[],3);

