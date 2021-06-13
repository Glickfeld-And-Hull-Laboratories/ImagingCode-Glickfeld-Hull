close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F1= fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure1');

ds = ['CrossOriRandDirTwoPhase_ExptList'];
eval(ds);
area = 'V1';
ind = find([expt.SF] == 0.05);

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
max_dir_all = [];

firstexp = 1;
for i = 1:length(ind)
    iexp = ind(i);
    
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
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
      
        plaid_SI_all = [plaid_SI_all plaid_SI(1,:)];
        plaid_SI_45_all = [plaid_SI_45_all plaid_SI_45(1,:)];
        
        resp_stim_prefDir_all = [resp_stim_prefDir_all resp_stim_prefDir];
        resp_mask_prefDir_all = [resp_mask_prefDir_all resp_mask_prefDir];
        resp_plaid_prefDir_all = [resp_plaid_prefDir_all resp_plaid_prefDir(1,:)];
        resp_stim_45Dir_all = [resp_stim_45Dir_all resp_stim_45Dir];
        resp_mask_45Dir_all = [resp_mask_45Dir_all resp_mask_45Dir];
        resp_plaid_45Dir_all = [resp_plaid_45Dir_all resp_plaid_45Dir(1,:)];
        
        avg_resp_dir_all = [avg_resp_dir_all; avg_resp_dir];
        [max_val max_dir] = max(avg_resp_dir(:,:,1,1),[],2);
        max_dir_all = [max_dir_all max_dir'];
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
        
        resp_any = find(sum(h_resp(:,:,1,1),2))';
        resp_any_ind = [resp_any_ind resp_any+totCells];
        resp_ind_all = [resp_ind_all resp_ind+totCells];
        resp_ind_45_all = [resp_ind_45_all resp_ind_45+totCells];
        
        if firstexp == 1
            prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
            postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn))+unique(celleqel2mat_padded(input.tItiWaitFrames));
            tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);
            data_resp = nan(prewin_frames+postwin_frames,nCells,nTrials-1);
            data_f = nan(1,nCells,nTrials-1);

            for itrial = 1:nTrials
                if cStimOn(itrial) + postwin_frames < sz(3)
                    data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
                end
            end
            data_f = mean(data_resp(1:prewin_frames,:,:),1);
            data_dfof_tc = (data_resp-data_f)./data_f;
            ind_stimAlone = intersect(find(stimCon_all),find(maskCon_all==0));
            ind_maskAlone = intersect(find(stimCon_all==0),find(maskCon_all));
            ind_plaid = intersect(find(stimCon_all),find(maskCon_all));
            ind_blank = intersect(find(stimCon_all==0),find(maskCon_all==0));
            resp_tc_blank = data_dfof_tc(:,:,ind_blank(1:end-1));
            ind_p = find(maskPhas_all == maskPhas(1));
            resp_tc_cell = cell(nStimDir,2);
            for iDir = 1:nStimDir
                ind_stimdir = find(stimDir_all == stimDirs(iDir));
                ind_maskdir = find(maskDir_all == stimDirs(iDir));
                ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
                ind_dirplaid = [intersect(ind_stimdir, ind_plaid)];
                resp_tc_cell{iDir,1} = data_dfof_tc(:,:,ind_diralone);
                ind_dpplaid = intersect(ind_dirplaid,ind_p);
                resp_tc_cell{iDir,2} = data_dfof_tc(:,:,ind_dpplaid);
            end
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']));
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']));
        end
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
        firstexp = 0;
     end
end
%%
plaid_SI_avg = mean(plaid_SI_all(resp_ind_all));
plaid_SI_sem = std(plaid_SI_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all));
plaid_SI_45_avg = mean(plaid_SI_45_all(resp_ind_45_all));
plaid_SI_45_sem = std(plaid_SI_45_all(resp_ind_45_all),[],2)./sqrt(length(resp_ind_45_all));

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
print(fullfile(summaryDir_F1, 'Figure1_SI_prefV45_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_SI_prefV45_lowSF.fig'))

ind_SI{1} = intersect(resp_ind_all,find(plaid_SI_all>0.3));
ind_SI{2} = intersect(resp_ind_all,find(plaid_SI_all<-0.5));
ind_SI{3} = intersect(resp_ind_all,intersect(find(plaid_SI_all>-0.1),find(plaid_SI_all<0.1)));
str = {'high','low','zero'};
for i =1:3
figure;
    for iC = 1:5
        iCell = ind_SI{i}(iC);
        subplot(5,4,1+((iC-1).*4))
        shadedErrorBar(tt, mean(resp_tc_blank(:,iCell,:),3), std(resp_tc_blank(:,iCell,:),[],3)./sqrt(size(resp_tc_blank,3)));
        ylim([-0.2 2.5])
        title(num2str(iCell))
        subplot(5,4,2+((iC-1).*4))
        shadedErrorBar(tt, mean(resp_tc_cell{max_dir_all(iCell),1}(:,iCell,:),3), std(resp_tc_cell{max_dir_all(iCell),1}(:,iCell,:),[],3)./sqrt(size(resp_tc_cell{max_dir_all(iCell),1},3)));
        ylim([-0.2 2.5])
        title('test')
        mask = max_dir_all(iCell)+4;
        if mask>nStimDir
            mask = mask-nStimDir;
        end
        subplot(5,4,3+((iC-1).*4))
        shadedErrorBar(tt, mean(resp_tc_cell{mask,1}(:,iCell,:),3), std(resp_tc_cell{mask,1}(:,iCell,:),[],3)./sqrt(size(resp_tc_cell{mask,1},3)));
        ylim([-0.2 2.5])
        title('mask')
        subplot(5,4,4+((iC-1).*4))
        shadedErrorBar(tt, mean(resp_tc_cell{max_dir_all(iCell),2}(:,iCell,:),3), std(resp_tc_cell{max_dir_all(iCell),2}(:,iCell,:),[],3)./sqrt(size(resp_tc_cell{max_dir_all(iCell),1},3)));
        ylim([-0.2 2.5])
        title(num2str(chop(plaid_SI_all(iCell),2)))
    end
    suptitle([str{i} ' SI'])
    print(fullfile(summaryDir_F1, ['Figure1_exampleCells_prefDir_lowSF_' str{i} ' SI.pdf']),'-dpdf','-bestfit')
    savefig(fullfile(summaryDir_F1, ['Figure1_exampleCells_prefDir_lowSF_' str{i} ' SI.fig']))
end   
    
date = expt(1).date;
mouse = expt(1).mouse;
[x y z] = scalebarCalib(str2num(date),'16x',[],1.7);
nCells = max(mask_cell(:));
stats = regionprops(mask_cell);
centroids = reshape([stats.Centroid],[2 nCells])';
figure; movegui('center')
imagesc(data_avg)
truesize
caxis([0 2500])
colormap(gray)
c1 = 13;
c2 = 21;
hold on
plot(centroids(c1,1), centroids(c1,2),'or','MarkerSize',10)
plot(centroids(c2,1), centroids(c2,2),'ob','MarkerSize',10)
title([date ' ' mouse '- x = ' num2str(chop(x,3)) 'um; y = ' num2str(chop(y,3)) 'um'])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCells_FOV.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCells_FOV.fig'))

