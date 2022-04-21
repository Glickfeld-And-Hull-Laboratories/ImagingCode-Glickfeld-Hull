clear all;clc;
mouse = 'WK19';
date = '220418';
time = char('1630');
ImgFolder = char('003');
RetImgFolder = char('001');



doFromRef = 1;
ref = char('001');

nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

fprintf(['2p_analysis imaging size tuning analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    % load 2p_analysis imaging data
  CD = ['Z:\home\Celine\Data\2p_data\' mouse '\' date '\' ImgFolder];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    % load behavior/experimental data
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    % read in frames with sbxread
    nframes = min([input.counterValues{end}(end) info.config.frames]);

    fprintf(['\nReading run ' num2str(irun) ' - ' num2str(nframes) ' frames \n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    toc
    
    fprintf('Reshaping data...\n')
    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
for i=1:length(input.tGratingContrast) % replace int64(con=1) with double
    if ~(class(input.tGratingContrast{i})=="double")
        input.tGratingContrast{i} = double(input.tGratingContrast{i});
    end
end
clear data_temp
clear temp

%% Choose register interval
regIntv = 1000;
nep = floor(size(data,3)./regIntv);
fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

% plot 500 frame means at each register interval
[n, n2] = subplotn(nep);
figure(1);clf;
colormap(gray)
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

%% Register data

chooseInt = 2; 

fprintf('\nBegin registering...\n')
if exist(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\2p_analysis_analysis\2p_analysis', mouse, date, ImgFolder), 'dir')
    % checks if analysis already present
    % load reg_shifts.mat (out, data_avg) and save the current input file
    fprintf('Found previous analysis! Loading...\n')
    
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_reg_shifts.mat']))
    
    % register
    fprintf('stackRegister_MA, using shifts from previous registration\n')
    % uses previous registration shifts (out) to re-register data quickly
    [outs, data_reg]=stackRegister_MA(data,[],[],double(out));
    fprintf('Previous registration loaded...\n')
    
    % save new input
    save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_input.mat']), 'input')
    
elseif doFromRef
    % if doFromRef specified, along with ref (ref needs to exist, no error catch)
    % load from ref folder:
    % reg_shifts.mat (out, data_avg)
    % mask_cell.mat ()
    % trialData.mat ()
    % then create new directory and save analysis
    fprintf(['Reference specified: ' ref '\n'])
    
    ref_str = ['runs-' ref];
    
    % someone put this to use multiple refs?
    %     ref_str = ['runs-' ref(1,:)];
    %     if size(ref,1)>1
    %         ref_str = [ref_str '-' ref(end,:)];
    %     end
    
    % load from folder specified by ref_str
    fprintf(['Loading from folder: ' ref_str '\n'])
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ref, [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ref, [date '_' mouse '_' ref_str '_mask_cell.mat']))
    
    % register
    fprintf('stackRegister with reference image\n')
    [out, data_reg] = stackRegister(data,data_avg);
    
    % save
    fprintf('Registration complete, now saving...\n')
    mkdir(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder))
    save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    % else means no previous analysis present
    % use data_avg selected above (could move here?)
    % then create new directory and save analysis
    fprintf('\nCreating new analysis!')
    
    meanrng = regIntv*(chooseInt)+(1:500);
    data_avg = mean(data(:,:,meanrng),3);
    fprintf(['\nRegister frame averaged from ' num2str(meanrng(1)) ' - ' num2str(meanrng(end)) '\n'])
    
    % register
    fprintf('stackRegister\n')
    [out, data_reg] = stackRegister(data,data_avg);
    
    % save
    fprintf('Registration complete, now saving...\n')
    mkdir(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder))
    save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data % depending on memory
%% check stability

regIntv = 2000;
nep = floor(size(data,3)./regIntv);
fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

% plot 500 frame means at each register interval
[n, n2] = subplotn(nep);

% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

% figure 3 shows imagesq (scaled?) of the registered data 1-10000, then prints to FOV_avg.mat
fprintf('Examine FOV\n')
figure(3);
if nframes>10000
    imagesq(mean(data_reg(:,:,1:10000),3));
else
    imagesq(mean(data_reg,3));
end
truesize;
% print to file
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_FOV_avg.pdf']))



%% import masks from retinotopy run and apply to the registered data

nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);
fprintf(['Loading masks from retinotopy runs: ' ret_str '\n'])

% loads 'mask_cell', 'mask_np'
load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse,date,RetImgFolder, [date '_' mouse '_' ret_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')



figure;
imagesc(mean(data_reg(:,:,:),3)); hold on;
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); 
colormap(gray)
% % extract dfof max image
% figure;
% imagesc(data_dfof_max); hold on;
% bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
% plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); 
% colormap(gray)

%% extract timecourses using the imported masks

 sz = size(data_reg);
fprintf('\nBegin time course extraction...\n')
down = 5; % downsample rate
data_reg_down  = stackGroupProject(data_reg,down);
fprintf(['Downsampling at M=' num2str(down) '\n'])

fprintf('Extracting cell signal...\n')
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
fprintf([num2str(nCells) ' total cells extracted\n'])

fprintf('Extracting neuropil signal...\n')
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
    np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
    fprintf(['Cell #' num2str(i) ' / ' num2str(nCells) '\n'])
end

fprintf('Subtract neuropil signal, maximizing skewness\n')
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
fprintf(['Maximum skews: ' num2str(max_skew) '\n'])
fprintf(['at inds: ' num2str(ind) '\n'])

np_w = 0.01*ind;
fprintf(['np_w = ' num2str(np_w) '\n'])

npSub_tc = data_tc-(tcRemoveDC(np_tc).*np_w);
clear data_reg_down

fprintf('Neuropil subtraction complete, saving data...\n')

fprintf('Calculating dF/F...\n')
%get dF/F
nCells = size(npSub_tc,2);
tc_mat = reshape(npSub_tc,[nOn+nOff nCells ntrials]);
% reshape into trials, could do this in one line with reshape?

tc_f = mean(tc_mat(nOff/2:nOff,:,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
clear tc_mat tc_f

fprintf('Time course extraction complete.\n')


save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc','tc_dfof')
save(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\celine\Analysis\2p_analysis', mouse, date, ImgFolder, [date '_' mouse '_' run_str '_input.mat']), 'input')




%% for each cell, find preferred orientation
% get average response to each orientation for each cell
% find max for each cell 


nOn = input.nScansOn;
nOff = input.nScansOff;

tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;

oris = unique(tOri);
nOri = length(oris);





