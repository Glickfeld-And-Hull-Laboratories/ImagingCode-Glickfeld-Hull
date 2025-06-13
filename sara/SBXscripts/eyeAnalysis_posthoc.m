clc; clear all; close all;
doRedChannel = 0;
% ds = 'CrossOriRandDirFourPhase_ExptList_SG';
ds = 'CrossOriRandPhase_15Hz_ExptList_SG';
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
eval(ds)
nexp = length(expt);

max_dist = 5;
frame_rate = 15;
seed = rng;

start=1;

for iexp = [92]
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';

    fprintf([mouse ' ' date '\n'])

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    
    figure;
    subplot(2,1,1)
        plot(centroid_dist)
        refline(0,5)
        subtitle('Centroid dist')
        xlabel('frames')
        ylabel('centroid dist')
    subplot(2,1,2)
        plot(Area)
        subtitle('Pupil area')
        xlabel('frames')
        ylabel('pupil area')

    
    
    sgtitle([ mouse ' ' date ' (exp 00' num2str(iexp) ')'])
    print(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],  [date '_' mouse '_' run_str '_EyeAnalysis.pdf']), '-dpdf','-fillpage')
end


%