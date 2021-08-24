%get path names
date = '200825';
ImgFolder = strvcat('003');
run = strvcat('003');
time = strvcat('1248');
mouse = 'i1328';
doFromRef = 0;
ref = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

%open TC data
load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\test\' date '_' mouse '\' date '_' mouse '_' run_str '\' date '_' mouse '_' run_str '_TCs.mat']);
load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\test\' date '_' mouse '\' date '_' mouse '_' run_str '\' date '_' mouse '_' run_str '_input.mat']);

%create necessary variables
tGratingDir = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(tGratingDir);
nDir = length(dirs);
nOn = input.nScansOn;
nOff = input.nScansOff;
nFrames = nOn+nOff;
nCells = size(npSub_tc,2);
nTrials = size(tGratingDir,2);

%transform into dF/F for each trial
trial_tc = reshape(npSub_tc,[nFrames nTrials nCells]);
trial_f = mean(trial_tc(nOff/2:nOff,:,:),1);
trial_dfof = (trial_tc-trial_f)./trial_f;

%divide into trial types and plot timecourses
iCell = 1;
[n n2] = subplotn(nDir);
dfof_dir = zeros(nFrames,nCells,nDir);

figure;
for idir = 1:nDir
    subplot(2.*n, n2,idir)
    ind = find(tGratingDir==dirs(idir));
    dfof_dir(:,:,idir) = squeeze(mean(trial_dfof(:,ind,:),2));
    plot(dfof_dir(:,iCell,idir))
    title(num2str(dirs(idir)))
    ylabel('dF/F')
    xlabel('Frames')
end  

%plot average response for each trial type
subplot(2,1,2)
plot(dirs,squeeze(mean(dfof_dir(nOff+1:end,iCell,:),1)),'-o')
xlabel('Direction (deg)')
ylabel('dF/F')
suptitle(['Cell #' num2str(iCell)])

print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' run_str '\directionTuning\directionTuning_Cell' num2str(iCell)],'-dpdf','-fillpage') 
