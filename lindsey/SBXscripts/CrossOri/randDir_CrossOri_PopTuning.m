close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F5= fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure5');

ds = ['CrossOriRandDirTwoPhase_ExptList'];
eval(ds);
area = 'V1';
ind = find([expt.SF] == 0.05);

nexp = size(expt,2);
totCells = 0;
totExp = 0;
exptInd = [];
resp_ind_all = [];
avg_resp_dir_all = [];
OSI_all = [];
max_dir_all = [];

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
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriFit.mat']))

        OSI_all = [OSI_all stim_OSI];
        
        avg_resp_dir_all = [avg_resp_dir_all; avg_resp_dir];
        [max_val max_dir] = max(avg_resp_dir(:,:,1,1),[],2);
        max_dir_all = [max_dir_all; max_dir];
        resp_ind = [];
        for iCell = 1:nCells
            dir_mask = max_dir(iCell)+4;
            if dir_mask>nStimDir
                dir_mask = dir_mask-nStimDir;
            end
            if sum(h_resp(iCell,[max_dir(iCell) dir_mask],1),2)
                resp_ind = [resp_ind iCell];
            end
        end
        
        resp_ind_all = [resp_ind_all resp_ind+totCells];
        
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
     end
end

%% by OSI
pop_tuning_all_align = zeros(nStimDir/2,3,2,4);
pop_tuning_highOSI_align = zeros(nStimDir/2,3,2,4);
pop_tuning_medOSI_align = zeros(nStimDir/2,3,2,4);
pop_tuning_lowOSI_align = zeros(nStimDir/2,3,2,4);
for ii = 1:4
pop_tuning_all = zeros(nStimDir/2,3,2);
pop_tuning_highOSI = zeros(nStimDir/2,3,2);
pop_tuning_medOSI = zeros(nStimDir/2,3,2);
pop_tuning_lowOSI = zeros(nStimDir/2,3,2);
for iDir = 1:nStimDir/2
    ind = intersect(resp_ind_all,find(max_dir_all == iDir));
    ind_n(iDir) = length(ind);
    pop_tuning_all(iDir,1,1) = mean(avg_resp_dir_all(ind,ii,1,1,1),1);
    pop_tuning_all(iDir,1,2) = std(avg_resp_dir_all(ind,ii,1,1,1),[],1)./sqrt(length(ind));
    pop_tuning_all(iDir,2,1) = mean(avg_resp_dir_all(ind,ii+4,1,1,1),1);
    pop_tuning_all(iDir,2,2) = std(avg_resp_dir_all(ind,ii+4,1,1,1),[],1)./sqrt(length(ind));
    pop_tuning_all(iDir,3,1) = mean(avg_resp_dir_all(ind,ii,1,2,1),1);
    pop_tuning_all(iDir,3,2) = std(avg_resp_dir_all(ind,ii,1,2,1),[],1)./sqrt(length(ind));
    indhigh = intersect(ind,find(OSI_all>0.7));
    indhigh_n(iDir) = length(indhigh);
    pop_tuning_highOSI(iDir,1,1) = mean(avg_resp_dir_all(indhigh,ii,1,1,1),1);
    pop_tuning_highOSI(iDir,1,2) = std(avg_resp_dir_all(indhigh,ii,1,1,1),[],1)./sqrt(length(indhigh));
    pop_tuning_highOSI(iDir,2,1) = mean(avg_resp_dir_all(indhigh,ii+4,1,1,1),1);
    pop_tuning_highOSI(iDir,2,2) = std(avg_resp_dir_all(indhigh,ii+4,1,1,1),[],1)./sqrt(length(indhigh));
    pop_tuning_highOSI(iDir,3,1) = mean(avg_resp_dir_all(indhigh,ii,1,2,1),1);
    pop_tuning_highOSI(iDir,3,2) = std(avg_resp_dir_all(indhigh,ii,1,2,1),[],1)./sqrt(length(indhigh));
    indmed = intersect(ind,intersect(find(OSI_all>=0.3),find(OSI_all<=0.7)));
    indmed_n(iDir) = length(indmed);
    pop_tuning_medOSI(iDir,1,1) = mean(avg_resp_dir_all(indmed,ii,1,1,1),1);
    pop_tuning_medOSI(iDir,1,2) = std(avg_resp_dir_all(indmed,ii,1,1,1),[],1)./sqrt(length(indmed));
    pop_tuning_medOSI(iDir,2,1) = mean(avg_resp_dir_all(indmed,ii+4,1,1,1),1);
    pop_tuning_medOSI(iDir,2,2) = std(avg_resp_dir_all(indmed,ii+4,1,1,1),[],1)./sqrt(length(indmed));
    pop_tuning_medOSI(iDir,3,1) = mean(avg_resp_dir_all(indmed,ii,1,2,1),1);
    pop_tuning_medOSI(iDir,3,2) = std(avg_resp_dir_all(indmed,ii,1,2,1),[],1)./sqrt(length(indmed));
    indlow = intersect(ind,find(OSI_all<0.3));
    indlow_n(iDir) = length(indlow);
    pop_tuning_lowOSI(iDir,1,1) = mean(avg_resp_dir_all(indlow,ii,1,1,1),1);
    pop_tuning_lowOSI(iDir,1,2) = std(avg_resp_dir_all(indlow,ii,1,1,1),[],1)./sqrt(length(indlow));
    pop_tuning_lowOSI(iDir,2,1) = mean(avg_resp_dir_all(indlow,ii+4,1,1,1),1);
    pop_tuning_lowOSI(iDir,2,2) = std(avg_resp_dir_all(indlow,ii+4,1,1,1),[],1)./sqrt(length(indlow));
    pop_tuning_lowOSI(iDir,3,1) = mean(avg_resp_dir_all(indlow,ii,1,2,1),1);
    pop_tuning_lowOSI(iDir,3,2) = std(avg_resp_dir_all(indlow,ii,1,2,1),[],1)./sqrt(length(indlow));
end

figure;
for i = 1:3
    subplot(3,4,1+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_all(:,i,1),pop_tuning_all(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    
    if i == 1
        title('All')
        ylabel([num2str(0+((ii-1).*22.5)) ' deg'])
    end
    if i == 2
        ylabel([num2str(90+((ii-1).*22.5)) ' deg'])
    end
    if i == 3
        ylabel('Plaid')
    end
    subplot(3,4,2+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_highOSI(:,i,1),pop_tuning_highOSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('OSI>0.7')
    end
    subplot(3,4,3+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_medOSI(:,i,1),pop_tuning_medOSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('0.3<OSI<0.7')
    end
    subplot(3,4,4+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_lowOSI(:,i,1),pop_tuning_lowOSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('OSI<0.3')
    end
    if i == 3
        subplot(3,4,1+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_all(:,1,1)+pop_tuning_all(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_all(:,1,1)+pop_tuning_all(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,2+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_highOSI(:,1,1)+pop_tuning_highOSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_highOSI(:,1,1)+pop_tuning_highOSI(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,3+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_medOSI(:,1,1)+pop_tuning_medOSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_medOSI(:,1,1)+pop_tuning_medOSI(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,4+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_lowOSI(:,1,1)+pop_tuning_lowOSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_lowOSI(:,1,1)+pop_tuning_lowOSI(:,2,1))/2,'-r')
        xlabel('Ori')
    end
end
pop_tuning_all_align(:,:,:,ii) = circshift(pop_tuning_all,1-ii,1);
pop_tuning_highOSI_align(:,:,:,ii) = circshift(pop_tuning_highOSI,1-ii,1);
pop_tuning_medOSI_align(:,:,:,ii) = circshift(pop_tuning_medOSI,1-ii,1);
pop_tuning_lowOSI_align(:,:,:,ii) = circshift(pop_tuning_lowOSI,1-ii,1);

end
pop_tuning_all_avg(:,:,1) = mean(pop_tuning_all_align(:,:,1,:),4);
pop_tuning_highOSI_avg(:,:,1) = mean(pop_tuning_highOSI_align(:,:,1,:),4);
pop_tuning_medOSI_avg(:,:,1) = mean(pop_tuning_medOSI_align(:,:,1,:),4);
pop_tuning_lowOSI_avg(:,:,1) = mean(pop_tuning_lowOSI_align(:,:,1,:),4);
pop_tuning_all_avg(:,:,2) = mean(pop_tuning_all_align(:,:,2,:),4);
pop_tuning_highOSI_avg(:,:,2) = mean(pop_tuning_highOSI_align(:,:,2,:),4);
pop_tuning_medOSI_avg(:,:,2) = mean(pop_tuning_medOSI_align(:,:,2,:),4);
pop_tuning_lowOSI_avg(:,:,2) = mean(pop_tuning_lowOSI_align(:,:,2,:),4);

figure;
for i = 1:3
    subplot(3,4,1+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_all_avg(:,i,1),pop_tuning_all_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    
    if i == 1
        title('All')
        ylabel([num2str(0+((ii-1).*22.5)) ' deg'])
    end
    if i == 2
        ylabel([num2str(90+((ii-1).*22.5)) ' deg'])
    end
    if i == 3
        ylabel('Plaid')
    end
    subplot(3,4,2+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_highOSI_avg(:,i,1),pop_tuning_highOSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('OSI>0.7')
    end
    subplot(3,4,3+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_medOSI_avg(:,i,1),pop_tuning_medOSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('0.3<OSI<0.7')
    end
    subplot(3,4,4+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_lowOSI_avg(:,i,1),pop_tuning_lowOSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('OSI<0.3')
    end
    if i == 3
        subplot(3,4,1+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_all_avg(:,1,1)+pop_tuning_all_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_all_avg(:,1,1)+pop_tuning_all_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,2+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_highOSI_avg(:,1,1)+pop_tuning_highOSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_highOSI_avg(:,1,1)+pop_tuning_highOSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,3+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_medOSI_avg(:,1,1)+pop_tuning_medOSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_medOSI_avg(:,1,1)+pop_tuning_medOSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,4+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_lowOSI_avg(:,1,1)+pop_tuning_lowOSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_lowOSI_avg(:,1,1)+pop_tuning_lowOSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
    end
end

%% by SI
avg_resp_dir_all_rect = avg_resp_dir_all;
avg_resp_dir_all_rect(find(avg_resp_dir_all<0)) = 0;

pop_tuning_all_align = zeros(nStimDir/2,3,2,4);
pop_tuning_highSI_align = zeros(nStimDir/2,3,2,4);
pop_tuning_medSI_align = zeros(nStimDir/2,3,2,4);
pop_tuning_lowSI_align = zeros(nStimDir/2,3,2,4);
for ii = 1:4
pop_tuning_all = zeros(nStimDir/2,3,2);
pop_tuning_highSI = zeros(nStimDir/2,3,2);
pop_tuning_medSI = zeros(nStimDir/2,3,2);
pop_tuning_lowSI = zeros(nStimDir/2,3,2);
stimSI_temp = abs((avg_resp_dir_all_rect(:,ii,1,1,1)-avg_resp_dir_all_rect(:,ii+4,1,1,1))./(avg_resp_dir_all_rect(:,ii,1,1,1)+avg_resp_dir_all_rect(:,ii+4,1,1,1)));
for iDir = 1:nStimDir/2
    ind = intersect(resp_ind_all,find(max_dir_all == iDir));
    ind_n(iDir) = length(ind);
    pop_tuning_all(iDir,1,1) = mean(avg_resp_dir_all(ind,ii,1,1,1),1);
    pop_tuning_all(iDir,1,2) = std(avg_resp_dir_all(ind,ii,1,1,1),[],1)./sqrt(length(ind));
    pop_tuning_all(iDir,2,1) = mean(avg_resp_dir_all(ind,ii+4,1,1,1),1);
    pop_tuning_all(iDir,2,2) = std(avg_resp_dir_all(ind,ii+4,1,1,1),[],1)./sqrt(length(ind));
    pop_tuning_all(iDir,3,1) = mean(avg_resp_dir_all(ind,ii,1,2,1),1);
    pop_tuning_all(iDir,3,2) = std(avg_resp_dir_all(ind,ii,1,2,1),[],1)./sqrt(length(ind));
    indhigh = intersect(ind,find(stimSI_temp>0.7));
    indhigh_n(iDir) = length(indhigh);
    pop_tuning_highSI(iDir,1,1) = mean(avg_resp_dir_all(indhigh,ii,1,1,1),1);
    pop_tuning_highSI(iDir,1,2) = std(avg_resp_dir_all(indhigh,ii,1,1,1),[],1)./sqrt(length(indhigh));
    pop_tuning_highSI(iDir,2,1) = mean(avg_resp_dir_all(indhigh,ii+4,1,1,1),1);
    pop_tuning_highSI(iDir,2,2) = std(avg_resp_dir_all(indhigh,ii+4,1,1,1),[],1)./sqrt(length(indhigh));
    pop_tuning_highSI(iDir,3,1) = mean(avg_resp_dir_all(indhigh,ii,1,2,1),1);
    pop_tuning_highSI(iDir,3,2) = std(avg_resp_dir_all(indhigh,ii,1,2,1),[],1)./sqrt(length(indhigh));
    indmed = intersect(ind,intersect(find(stimSI_temp>=0.3),find(stimSI_temp<=0.7)));
    indmed_n(iDir) = length(indmed);
    pop_tuning_medSI(iDir,1,1) = mean(avg_resp_dir_all(indmed,ii,1,1,1),1);
    pop_tuning_medSI(iDir,1,2) = std(avg_resp_dir_all(indmed,ii,1,1,1),[],1)./sqrt(length(indmed));
    pop_tuning_medSI(iDir,2,1) = mean(avg_resp_dir_all(indmed,ii+4,1,1,1),1);
    pop_tuning_medSI(iDir,2,2) = std(avg_resp_dir_all(indmed,ii+4,1,1,1),[],1)./sqrt(length(indmed));
    pop_tuning_medSI(iDir,3,1) = mean(avg_resp_dir_all(indmed,ii,1,2,1),1);
    pop_tuning_medSI(iDir,3,2) = std(avg_resp_dir_all(indmed,ii,1,2,1),[],1)./sqrt(length(indmed));
    indlow = intersect(ind,find(stimSI_temp<0.3));
    indlow_n(iDir) = length(indlow);
    pop_tuning_lowSI(iDir,1,1) = mean(avg_resp_dir_all(indlow,ii,1,1,1),1);
    pop_tuning_lowSI(iDir,1,2) = std(avg_resp_dir_all(indlow,ii,1,1,1),[],1)./sqrt(length(indlow));
    pop_tuning_lowSI(iDir,2,1) = mean(avg_resp_dir_all(indlow,ii+4,1,1,1),1);
    pop_tuning_lowSI(iDir,2,2) = std(avg_resp_dir_all(indlow,ii+4,1,1,1),[],1)./sqrt(length(indlow));
    pop_tuning_lowSI(iDir,3,1) = mean(avg_resp_dir_all(indlow,ii,1,2,1),1);
    pop_tuning_lowSI(iDir,3,2) = std(avg_resp_dir_all(indlow,ii,1,2,1),[],1)./sqrt(length(indlow));
end

figure;
for i = 1:3
    subplot(3,4,1+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_all(:,i,1),pop_tuning_all(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    
    if i == 1
        title('All')
        ylabel([num2str(0+((ii-1).*22.5)) ' deg'])
    end
    if i == 2
        ylabel([num2str(90+((ii-1).*22.5)) ' deg'])
    end
    if i == 3
        ylabel('Plaid')
    end
    subplot(3,4,2+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_highSI(:,i,1),pop_tuning_highSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('SI>0.7')
    end
    subplot(3,4,3+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_medSI(:,i,1),pop_tuning_medSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('0.3<SI<0.7')
    end
    subplot(3,4,4+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_lowSI(:,i,1),pop_tuning_lowSI(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('SI<0.3')
    end
    if i == 3
        subplot(3,4,1+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_all(:,1,1)+pop_tuning_all(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_all(:,1,1)+pop_tuning_all(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,2+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_highSI(:,1,1)+pop_tuning_highSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_highSI(:,1,1)+pop_tuning_highSI(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,3+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_medSI(:,1,1)+pop_tuning_medSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_medSI(:,1,1)+pop_tuning_medSI(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,4+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_lowSI(:,1,1)+pop_tuning_lowSI(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_lowSI(:,1,1)+pop_tuning_lowSI(:,2,1))/2,'-r')
        xlabel('Ori')
    end
end
pop_tuning_all_align(:,:,:,ii) = circshift(pop_tuning_all,1-ii,1);
pop_tuning_highSI_align(:,:,:,ii) = circshift(pop_tuning_highSI,1-ii,1);
pop_tuning_medSI_align(:,:,:,ii) = circshift(pop_tuning_medSI,1-ii,1);
pop_tuning_lowSI_align(:,:,:,ii) = circshift(pop_tuning_lowSI,1-ii,1);

end
pop_tuning_all_avg(:,:,1) = mean(pop_tuning_all_align(:,:,1,:),4);
pop_tuning_highSI_avg(:,:,1) = mean(pop_tuning_highSI_align(:,:,1,:),4);
pop_tuning_medSI_avg(:,:,1) = mean(pop_tuning_medSI_align(:,:,1,:),4);
pop_tuning_lowSI_avg(:,:,1) = mean(pop_tuning_lowSI_align(:,:,1,:),4);
pop_tuning_all_avg(:,:,2) = mean(pop_tuning_all_align(:,:,2,:),4);
pop_tuning_highSI_avg(:,:,2) = mean(pop_tuning_highSI_align(:,:,2,:),4);
pop_tuning_medSI_avg(:,:,2) = mean(pop_tuning_medSI_align(:,:,2,:),4);
pop_tuning_lowSI_avg(:,:,2) = mean(pop_tuning_lowSI_align(:,:,2,:),4);

figure;
for i = 1:3
    subplot(3,4,1+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_all_avg(:,i,1),pop_tuning_all_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    
    if i == 1
        title('All')
        ylabel([num2str(0+((ii-1).*22.5)) ' deg'])
    end
    if i == 2
        ylabel([num2str(90+((ii-1).*22.5)) ' deg'])
    end
    if i == 3
        ylabel('Plaid')
    end
    subplot(3,4,2+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_highSI_avg(:,i,1),pop_tuning_highSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('SI>0.7')
    end
    subplot(3,4,3+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_medSI_avg(:,i,1),pop_tuning_medSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('0.3<SI<0.7')
    end
    subplot(3,4,4+((i-1)*4))
    errorbar(stimDirs(1:nStimDir/2),pop_tuning_lowSI_avg(:,i,1),pop_tuning_lowSI_avg(:,i,2),'o')
    ylim([0 0.6])
    xlim([-22.5 180])
    if i == 1
        title('SI<0.3')
    end
    if i == 3
        subplot(3,4,1+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_all_avg(:,1,1)+pop_tuning_all_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_all_avg(:,1,1)+pop_tuning_all_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,2+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_highSI_avg(:,1,1)+pop_tuning_highSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_highSI_avg(:,1,1)+pop_tuning_highSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,3+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_medSI_avg(:,1,1)+pop_tuning_medSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_medSI_avg(:,1,1)+pop_tuning_medSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
        subplot(3,4,4+((i-1)*4))
        hold on
        plot(stimDirs(1:nStimDir/2),pop_tuning_lowSI_avg(:,1,1)+pop_tuning_lowSI_avg(:,2,1),'-k')
        plot(stimDirs(1:nStimDir/2),(pop_tuning_lowSI_avg(:,1,1)+pop_tuning_lowSI_avg(:,2,1))/2,'-r')
        xlabel('Ori')
    end
end

