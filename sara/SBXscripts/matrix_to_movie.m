clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
iexp = 52; 
doPhaseAfterDir = 0;
doDirAfterPass = 0;
eval(ds)

frame_rate = params.frameRate;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
if doPhaseAfterDir
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    ref_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).copFolder;
    time = expt(iexp).copTime;
elseif doDirAfterPass
    ImgFolder = expt(iexp).passFolder;
    nrun = length(ImgFolder);
    ref_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
else
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
end
nrun = length(ImgFolder);
run_str = catRunName(ImgFolder, nrun);

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
data_avg = mean(data(:,:,30001:30500),3);
if doPhaseAfterDir || doDirAfterPass
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
elseif exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
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
movegui('center')
figure; imagesq(mean(data_reg(:,:,1:regIntv),3)); truesize;

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

%%
data_reg_mov = uint8(data_reg/256);

nMovieFrames = 50;
cmap=colormap("winter");

for in = 1:nMovieFrames
    movie_matrix = data_reg_mov(:,:,in);
    imwrite(movie_matrix,fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'movie_test_folder',  '2P_movie_test.tif'), 'Compression','none','WriteMode','append');
end


%% failed video writing code

T = uint8(data_reg/256);

VidObj = VideoWriter('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Analysis\2P\230127_i1373\movie_test2.avi', 'Uncompressed AVI'); %set your file name and video compression
VidObj.FrameRate = 15; %set your frame rate
open(VidObj);

for f = 45:110   %1:size(T, 3)
    writeVideo(VidObj, T(:, :, f));
end
close(VidObj);




