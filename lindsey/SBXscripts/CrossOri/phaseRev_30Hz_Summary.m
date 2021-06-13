clear all;
close all;
clc
doRedChannel = 0;
ds = 'CrossOriRandPhase_lowSF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'PhaseRev');
phaseRev = struct;
    
totCells = 0;
resp_ind_all = [];
f1_all = [];
f2_all = [];
f2overf1_all = [];
mouse_list = [];

for iexp = 1:nexp
    if ~isempty(expt(iexp).prFolder)
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        ImgFolder = expt(iexp).coFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        fprintf([mouse ' ' date '\n'])

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        nCells = size(data_dfof_tc,2);
        resp_ind_all = [resp_ind_all resp_ind'+totCells];

        ImgFolder = expt(iexp).prFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']))
        f1_all = [f1_all f1];
        f2_all = [f2_all f2];
        f2overf1_all = [f2overf1_all f2overf1];
        totCells = totCells+nCells;
    end
end
phaseRev.nExpt = size(mouse_list,1);
phaseRev.nMice = size(unique(mouse_list,'rows'),1);
phaseRev.totCells = totCells;
phaseRev.respCells = length(resp_ind_all);
phaseRev.threshCells = length(find(f1_all(resp_ind_all)>0.02));

save(fullfile(summaryDir,['phaseRev30Hz_Summary.mat']),'resp_ind_all', 'f1_all','f2_all','f2overf1_all','mouse_list')

phaseRevTable = struct2table(phaseRev);
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
stim.frameRate = frame_rate;
stim.nOnSec = double(input.nScansOn)./stim.frameRate;
stim.nOffSec = double(input.nScansOff)./stim.frameRate;
stim.phaseRevHz = stim.frameRate./double(input.nScansPhaseCyc.*2);
stim.gratingDiameter = input.gratingDiameterDeg;
stim.gratingSF = input.gratingSpatialFreqCPD;
stim.gratingOriStep = input.gratingDirectionStepDeg;
stim.gratingOriStepN = input.gratingDirectionStepN;
stim.gratingPhaseStep = input.gratingStartingPhaseStepDeg;
stim.gratingPhaseStepN = input.gratingStartingPhaseStepN;

save(fullfile(summaryDir,['phaseRevStats_V1only_30Hz.mat']),'stim','phaseRev');
