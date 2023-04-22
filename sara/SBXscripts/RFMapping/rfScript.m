clc; clear all; close all;
ds = 'RFMapping_15Hz_ExptList_SG';
iexp = 3; 
eval(ds)

frame_rate = params.frameRate;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};

ImgFolder = expt(iexp).rfFolder;
time = expt(iexp).rfTime;

nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002'

base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];


fprintf(['2P imaging cross ori analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cStimOneOn') 
        if irun>1
            ntrials = size(input.cStimOneOn,2);
            for itrial = 1:ntrials
                temp(irun).cStimOneOn{itrial} = temp(irun).cStimOneOn{itrial}+offset;
                temp(irun).cStimOneOff{itrial} = temp(irun).cStimOneOff{itrial}+offset;
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateStructuresLG(temp);
clear data_temp
clear temp
toc

% Choose register interval
regIntv = 5000;
nep = floor(size(data,3)./regIntv);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); colormap gray; clim([0 3000]); end
movegui('center')
%% Register data
data_avg = mean(data(:,:,20001:20500),3);

if exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg] = stackRegister_MA(data,[],[],out);
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end

% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3)); title([num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]); end
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
movegui('center')
figure; imagesq(mean(data_reg(:,:,1:regIntv),3)); truesize;
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

i = 1;
sz = size(data_reg);
rg = zeros(sz(1),sz(2),3);
first = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3);
rg(:,:,1) = first./max(first(:));
i = nep; 
last = mean(data_reg(:,:,1+((i-1)*regIntv):500+((i-1)*regIntv)),3);
rg(:,:,2) = last./max(last(:));
figure; image(rg)
movegui('center')

%% find activated cells
clear data out
nOn = input.nScansOn;
nOff = input.nScansOff;
nTrials = length(input.tTrialsDoneSinceStart);
sz = size(data_reg);

data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, nTrials]);

data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_resp = mean(data_tr(:,:,nOff+1:nOn+nOff,:),3); %are these the windows we want?

data_resp_dfof = (data_resp-data_f)./data_f;


stimEl = celleqel2mat_padded(input.tGratingElevationDeg);
El = unique(stimEl);
nEl = length(El);
stimAz = celleqel2mat_padded(input.tGratingAzimuthDeg);
Az = unique(stimAz);
nAz = length(Az);
stimPhas = celleqel2mat_padded(input.tGratingStartingPhaseDeg);
Phas = unique(stimPhas);
nPhas = length(Phas);




nStim = nEl*nAz*nPhas;
data_dfof = nan(sz(1),sz(2),nStim);
stim_list = nan(nStim,3);
start = 1;

for ip = 1:nPhas
    ind_p = find(stimPhas == Phas(ip));
    figure;
    t = 1;
    for iE = 1:nEl
        ind_El = find(stimEl == El(iE));
        for iA = 1:nAz
            ind_Az = find(stimAz == Az(iA));
            ind_all = intersect(ind_p,intersect(ind_Az,ind_El));
            data_dfof(:,:,start) = nanmean(data_resp_dfof(:,:,ind_all),3);
            subplot(nAz,nEl,t)
            imagesc(data_dfof(:,:,start))
            stim_list(start,:) = [Phas(ip) El(iE) Az(iA)];
            start = start+1;
            t = t+1;
        end
    end
    sgtitle(num2str(Phas(ip)))
    movegui('center')
end

szSq = 3; %size of square side length
nSq = (nAz*nEl)/(szSq*szSq);
if rem(nAz, szSq) ~= 0
    fprintf('*Error: mismatched square size*\n')
    stop
elseif rem(nEl, szSq) ~= 0
    fprintf('*Error: mismatched square size*\n')
    stop
end



stimSq = flip(reshape(1:(nEl*nAz), nEl, nAz).',1);
stimlist = [];
ii = 0;
for i = 1:szSq:nEl
    stimlist = [stimlist, stimSq(i:(szSq + ii),:)];
    ii = ii + szSq;
end
stimSqs = reshape(stimlist, szSq, szSq, []);
listSqs = reshape(stimSqs, szSq*szSq, []);
listSqs = [listSqs listSqs + nEl*nAz];

figure; 
for i = 1:nPhas*nSq
    data_dfof(:,:,i) = mean(data_dfof(:,:,listSqs(:,i)),3); 
    subplot(nPhas,nSq,i)
    imagesc(mean(data_dfof(:,:,listSqs(:,i)),3))
    movegui('center')
end
       
data_dfof(:,:,isnan(mean(mean(data_dfof,1),2))) = [];
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
movegui('center')
imagesc(data_dfof_max)
print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_maxdfof.pdf']),'-dpdf')

data_dfof = cat(3, data_dfof, data_dfof_max);

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'nOn', 'nOff', 'stimEl', 'stimAz', 'stimPhas', 'El', 'Az', 'Phas', 'nEl', 'nAz', 'nPhas', 'nStim', 'stim_list', 'listSqs', 'frame_rate', 'nTrials')

%% cell segmentation 


mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

    for iStim = 1:size(data_dfof,3)   
        mask_data_temp = mask_data(:,:,end+1-iStim);
        mask_data_temp(find(mask_exp >= 1)) = 0;
        bwout = imCellEditInteractiveLG(mask_data_temp);
        mask_all = mask_all+bwout;
        mask_exp = imCellBuffer(mask_all,3)+mask_all;
        close all
    end
    
    mask_cell= bwlabel(mask_all);
    figure; movegui('center')
    imagesc(mask_cell)
    

clear data_adapt data_adapt_dfof data_test data_test_dfof data_test_avg data_resp data_resp_dfof bwout
%% neuropil subtraction

    mask_np = imCellNeuropil(mask_cell, 3, 5);
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 


down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell


 