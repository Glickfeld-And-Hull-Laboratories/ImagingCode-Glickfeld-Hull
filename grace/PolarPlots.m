% Used to somehow make polar plots of direction tuning. I haven't really
% fine tuned this yet, but it could be useful.
% load and define variables
date1 = '111419';
date2 = '111819';
run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

day1DT = load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date1 '_' mouse '\' date1 '_' mouse '_' run_str '\' date1 '_' mouse '_' run_str '_DirTuning.mat'],'DT');
dat2DT = load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date2 '_' mouse '\' date2 '_' mouse '_' run_str '\' date2 '_' mouse '_' run_str '_DirTuning.mat'],'DT');

y_fit1 = day1DT.DT.y_fit;
dfof_dir1 = dat1DT.DT.dfof_dir;
dirs1 = day1DT.DT.dirs;
idx1 = day1DT.DT.idx;
dfof_polar1 = day1DT.DT.dfof_polar;

y_fit2 = day2DT.DT.y_fit;
dfof_dir2 = dat2DT.DT.dfof_dir;
dirs2 = day2DT.DT.dirs;
idx2 = day2DT.DT.idx;
dfof_polar2 = day2DT.DT.dfof_polar;

% plot day 2 over day 1
for i_cell = 1:length(idx1)
    iCell=idx1(i_cell);
 figure;
    if min(dfof_dir1(:,iCell,1)) < 0
        dfof_polar1(:,iCell,1) = squeeze(dfof_dir1(:,iCell,1))-squeeze(min(dfof_dir1(:,iCell,1)));
        polarplot([theta 2*pi],[(dfof_polar1(:,iCell)); (dfof_polar1(1,iCell))])
        suptitle(['Cell #' num2str(iCell)])
    else
        polarplot([theta 2*pi],[squeeze(dfof_dir1(:,iCell,1)); squeeze(dfof_dir1(1,iCell,1))])
        suptitle(['Cell #' num2str(iCell)])
    end
hold on
    if min(dfof_dir2(:,iCell,1)) < 0
        dfof_polar2(:,iCell,1) = squeeze(dfof_dir2(:,iCell,1))-squeeze(min(dfof_dir2(:,iCell,1)));
        polarplot([theta 2*pi],[(dfof_polar2(:,iCell)); (dfof_polar2(1,iCell))])
        suptitle(['Cell #' num2str(iCell)])
    else
        polarplot([theta 2*pi],[squeeze(dfof_dir2(:,iCell,1)); squeeze(dfof_dir2(1,iCell,1))])
        suptitle(['Cell #' num2str(iCell)])
    end
end

