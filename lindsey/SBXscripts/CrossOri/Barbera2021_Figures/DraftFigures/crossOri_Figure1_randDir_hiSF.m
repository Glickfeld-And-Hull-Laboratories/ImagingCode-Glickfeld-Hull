close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F1= fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure1');

ds = ['CrossOriRandDir_ExptList'];
eval(ds);
area = 'V1';

nexp = size(expt,2);
totCells = 0;
totExp = 0;
exptInd = [];
resp_any_ind = [];
resp_ind_all = [];
resp_ind_45_all = [];
plaid_SI_all = [];
plaid_SI_45_all = [];
resp_stim_prefDir_all = [];
resp_mask_prefDir_all = [];
resp_plaid_prefDir_all = [];
resp_stim_45Dir_all = [];
resp_mask_45Dir_all = [];
resp_plaid_45Dir_all = [];
avg_resp_dir_all = [];
Zc_all = [];
Zp_all = [];
OSI_all = [];

for iexp = 1:nexp
    
     if strcmp(expt(iexp).img_loc,area)
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).coFolder;
        time = expt(iexp).coTime;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);

        fprintf([mouse ' ' date '\n'])

        %% load data

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))

      
        plaid_SI_all = [plaid_SI_all plaid_SI];
        plaid_SI_45_all = [plaid_SI_45_all plaid_SI_45];
        OSI_all = [OSI_all stim_OSI];
        
        resp_stim_prefDir_all = [resp_stim_prefDir_all resp_stim_prefDir];
        resp_mask_prefDir_all = [resp_mask_prefDir_all resp_mask_prefDir];
        resp_plaid_prefDir_all = [resp_plaid_prefDir_all resp_plaid_prefDir];
        resp_stim_45Dir_all = [resp_stim_45Dir_all resp_stim_45Dir];
        resp_mask_45Dir_all = [resp_mask_45Dir_all resp_mask_45Dir];
        resp_plaid_45Dir_all = [resp_plaid_45Dir_all resp_plaid_45Dir];
        
        avg_resp_dir_all = [avg_resp_dir_all; avg_resp_dir];
        [max_val max_dir] = max(avg_resp_dir(:,:,1),[],2);
        resp_ind = [];
        resp_ind_45 = [];
        for iCell = 1:nCells
            dir_mask = max_dir(iCell)+4;
            if dir_mask>nStimDir
                dir_mask = dir_mask-nStimDir;
            end
            if sum(h_resp(iCell,[max_dir(iCell) dir_mask],1),2)
                resp_ind = [resp_ind iCell];
            end
            dir_45 = max_dir(iCell)+2;
            if dir_45>nStimDir
                dir_45 = dir_45-nStimDir;
            end
            dir_mask_45 = max_dir(iCell)+6;
            if dir_mask_45>nStimDir
                dir_mask_45 = dir_mask_45-nStimDir;
            end
            if sum(h_resp(iCell,[dir_45 dir_mask_45],1),2)
                resp_ind_45 = [resp_ind_45 iCell];
            end
        end
        
        resp_any = find(sum(h_resp(:,:,1),2))';
        resp_any_ind = [resp_any_ind resp_any+totCells];
        resp_ind_all = [resp_ind_all resp_ind+totCells];
        resp_ind_45_all = [resp_ind_45_all resp_ind_45+totCells];
        Zc_all = [Zc_all Zc];
        Zp_all = [Zp_all Zp];
        
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
     end
end

figure;
subplot(3,2,1)
scatter(resp_stim_prefDir_all(resp_ind_all)+resp_mask_prefDir_all(resp_ind_all),resp_plaid_prefDir_all(resp_ind_all),'ok')
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Test + Mask')
ylabel('Plaid')
title(['Pref dir: n = ' num2str(length(resp_ind_all))])
subplot(3,2,2)
scatter(resp_stim_45Dir_all(resp_ind_45_all)+resp_mask_45Dir_all(resp_ind_45_all),resp_plaid_45Dir_all(resp_ind_45_all),'ok')
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Test + Mask')
ylabel('Plaid')
title(['45 deg from pref: n = ' num2str(length(resp_stim_45Dir_all(resp_ind_45_all)))])
subplot(3,2,3)
histogram(plaid_SI_all(resp_ind_all),-1:0.2:1)
xlim([-1 1])
xlabel('Selectivity index')
ylabel('Number of cells')
subplot(3,2,4)
histogram(plaid_SI_45_all(resp_ind_45_all),-1:0.2:1)
xlim([-1 1])
xlabel('Selectivity index')
ylabel('Number of cells')
subplot(3,2,5)
cdfplot(plaid_SI_all(resp_ind_all))
xlim([-1 1])
xlabel('Selectivity index')
ylabel('Fraction of cells')
title('')
subplot(3,2,6)
cdfplot(plaid_SI_45_all(resp_ind_45_all))
xlim([-1 1])
xlabel('Selectivity index')
ylabel('Fraction of cells')
title('')
expts = unique(exptInd);
nexp_area = length(expts);
mouse_list = [{expt.mouse}];
mice = unique(mouse_list(expts));
nmice = length(mice);
suptitle([num2str(nexp_area) ' expts; ' num2str(nmice) ' mice; ' num2str(length(resp_ind_all)) ' cells'])
print(fullfile(summaryDir_F1, 'Figure1_SI_prefV45_hiSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_SI_prefV45_hiSF.fig'))

