%% Plotting timecourses, gettting preferred phase...

clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;
max_dist = 5;
    

%% plot timecourses with shaded error bar; one cell per page

for iexp = [4]; 
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = 'runs-002';

    base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara'];

    fprintf([mouse ' ' date '\n'])

% Load timecourse data, stim data
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))

    clear data_tc data_dfof_con_ph_tc_avg base_cell np_tc npSub_tc Val

    fprintf(['\nMouse: ' mouse '\nDate: ' date '\n'])
    fprintf(['nCells = ' num2str(nCells) '\n'])
    fprintf(['nCells resp = ' num2str(length(resp_ind)) '\n'])

    phases = [0:(360/nMaskPhas):360];

    %Plot individual and average time courses for all phases for cells
    
    for ic = [1]; %cell number %[rand_ind']
        figure;
        s=1;
        for ip = 1:nMaskPhas
            ptrials = trialInd{2,2,ip}; %set which specific plaid (1-8) you want to look at
            ptrials = setdiff(ptrials, find(centroid_dist>max_dist));
                X = [1:30]; %taken from how data_dfof_tc is calculated in crossori_expt script (should equal on + off frames)
                dfof_trialavg = mean(data_dfof_tc(:,ic,ptrials),3);
                numtrials = size(squeeze(data_dfof_tc(:,ic,ptrials)),2);
                dfof_trialstd = std(squeeze(data_dfof_tc(:,ic,ptrials)),0,2);
                dfof_trialsem = dfof_trialstd ./ sqrt(numtrials);
           subplot(5,12,ip)
               plot(X,dfof_trialavg, 'b','LineWidth',1);
               shadedErrorBar(X,dfof_trialavg,dfof_trialsem);
               xline(15,'--b')
                ylim([-0.2 1])
               xlabel('frames')
               ylabel('df/f')
               set(gca,'TickDir','out'); box off; 
               title(['cell #' num2str(ic), ' ' num2str(ip)])
               movegui('center')
               s=s+1;
        end
        print(fullfile(base, 'Analysis\2P\CrossOri\timecourses', [date '_' mouse '_' run_str '_TCs_maxdist' num2str(max_dist) '_Cell' num2str(ic) '.pdf']), '-dpdf','-fillpage')
        close all
    end
end

