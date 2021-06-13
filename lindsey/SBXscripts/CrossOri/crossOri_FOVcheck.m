close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%combine all RandomPhase at 0.05 cpd
ds = ['CrossOriRandPhase_lowSF_ExptList'];
eval(ds);
expt1 = expt;
clear expt
ds = ['CrossOriRandPhase_ExptList'];
eval(ds);
expt2 = expt;

expt_all = combineStructures(expt1,expt2);
mouse_cell = {expt_all.mouse};
mice  = unique(mouse_cell);
nmice = length(mice);
depth_mat = cell2mat({expt_all.z});

%check to make sure FOVs are unique
for imouse = 1:length(mice)
    ind = find(strcmp(mouse_cell, mice(imouse)));
    nexp = length(ind);
    if nexp>1
        mouse = cell2mat(mice(imouse));
        [n n2] = subplotn(nexp);
        figure;
        for iexp = 1:nexp
            subplot(n,n2,iexp);
            i = ind(iexp);
            date = expt_all(i).date;
            ImgFolder = expt_all(i).coFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']);
            load(fn)
            imagesc(data_avg)
            title([date ' z = ' num2str(expt_all(i).z)])
        end
        suptitle(mouse)
    end
end