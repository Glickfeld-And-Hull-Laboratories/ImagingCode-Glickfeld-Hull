close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F8 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure8');

ds = ['CrossOriRandDirTwoPhase_ExptList'];
eval(ds);
area = 'V1';
ind = find([expt.SF] == 0.05);

nexp = size(expt,2);
totExp = 0;
exptInd = [];
Zc_all = [];
Zp_all = [];
Zc_phase2_all = [];
Zp_phase2_all = [];
plaid_SI_all = [];
h_plaid_SI_all = [];
totCells = 0;
resp_ind_all = [];
resp_ind_dir_all = [];
resp_ind_plaid_all = [];
avg_resp_dir_all = [];
pattern_all = [];
component_all = [];
OSI_all = [];
k1_all = [];

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
        
        fprintf(['n = ' num2str(nCells) '\n'])
        
        Zc_all = [Zc_all Zc(1,:)];
        Zp_all = [Zp_all Zp(1,:)];
        Zc_phase2_all = [Zc_phase2_all Zc(2,:)];
        Zp_phase2_all = [Zp_phase2_all Zp(2,:)];
        OSI_all = [OSI_all stim_OSI(1,:)];
        plaid_SI_all = [plaid_SI_all plaid_SI(1,:)];
        h_plaid_SI_all = [h_plaid_SI_all h_plaid_SI(1,:)];
        
        resp_ind = find(sum(sum(sum(h_resp,2),3),4));
        resp_ind_dir = find(sum(h_resp(:,:,1,1),2));
        resp_ind_plaid = find(sum(sum(h_resp(:,:,:,2),2),3));

        resp_ind_all = [resp_ind_all resp_ind'+totCells];
        resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
        resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];
        
        avg_resp_dir_all = cat(1,avg_resp_dir_all, avg_resp_dir); 
        % nCells nDir nPhase maskCon mean/std
        component_all = [component_all; component];
        pattern_all = [pattern_all; pattern];
                
        avg_resp_ori = mean(reshape(avg_resp_dir(:,:,1,1,1),[nCells, nStimDir/2, 2]),3);
        oris = stimDirs(1:nStimDir/2);
        b_ori = nan(1,nCells);
        k1_ori = nan(1,nCells);
        R1_ori = nan(1,nCells);
        u1_ori = nan(1,nCells);
        sse_ori = nan(1,nCells);
        R_square_ori = nan(1,nCells);
        range = 1:1:180;
        y_ori_fit = nan(nCells,length(range));
        
        for iCell = 1:nCells
            data = [avg_resp_ori(iCell,:) avg_resp_ori(iCell,1)];
            theta = [deg2rad(oris) pi];
            [b_ori(:,iCell),k1_ori(:,iCell),R1_ori(:,iCell),u1_ori(:,iCell),sse_ori(:,iCell),R_square_ori(:,iCell)] ...
                = miaovonmisesfit_ori(theta,data);
            y_ori_fit(iCell,:) = b_ori(:,iCell)+R1_ori(:,iCell).*exp(k1_ori(:,iCell).*(cos(2.*(deg2rad(range)-u1_ori(:,iCell)))-1));
        end
        save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriFit.mat']), 'b_ori','k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori', 'y_ori_fit','range')
        k1_all = [k1_all k1_ori];
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
     end
end
%%
avg_resp_dir_all_circ = squeeze(cat(2,avg_resp_dir_all(:,:,1,:,:), avg_resp_dir_all(:,1,1,:,:)));
avg_resp_dir_all_circ(find(avg_resp_dir_all_circ<0)) = 0;
avg_resp_dir_all_shift = circshift(avg_resp_dir_all,2,2);
avg_resp_dir_all_circ_shift = squeeze(cat(2,avg_resp_dir_all_shift(:,:,1,:,:), avg_resp_dir_all_shift(:,1,1,:,:)));
avg_resp_dir_all_circ_shift(find(avg_resp_dir_all_circ_shift<0)) = 0;
component_all_shift = circshift(component_all,2,2);
component_all_circ = cat(2,component_all_shift, component_all_shift(:,1));
component_all_circ(find(component_all_circ<0)) = 0;
stimDirs_circ = [stimDirs stimDirs(1)];

[max_val max_ind] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
max_val(find(max_val<0)) = 0;
nDir = size(avg_resp_dir_all,2);
null_ind = max_ind+(nDir/2);
null_ind(find(null_ind>nDir)) = null_ind(find(null_ind>nDir))-nDir;
null_val = avg_resp_dir_all(:,null_ind,1,1,1);
null_val(find(null_val<0)) = 0;
DSI_all = (max_val-null_val)./(max_val+null_val);

figure; 
movegui('center')
subplot(2,2,1)
scatter(OSI_all(resp_ind_all), Zc_all(resp_ind_all));
xlabel('OSI')
ylabel('Zc')
ylim([-5 10])
subplot(2,2,2)
scatter(k1_all(resp_ind_all), Zc_all(resp_ind_all));
xlabel('k')
ylabel('Zc')
ylim([-5 10])
subplot(2,2,3)
scatter(OSI_all(resp_ind_all), Zp_all(resp_ind_all));
xlabel('OSI')
ylabel('Zp')
ylim([-5 10])
subplot(2,2,4)
scatter(k1_all(resp_ind_all), Zp_all(resp_ind_all));
xlabel('k')
ylabel('Zp')
ylim([-5 10])
print(fullfile(summaryDir_F8, 'Figure8_ZpZc_vs_OSIk_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ZpZc_vs_OSIk_lowSF.fig'))


figure; 
movegui('center')
start = 1;
for i = 1:16
    subplot(4,4,i)
    iC = resp_ind_all(i);
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1) avg_resp_dir_all_circ(iC,:,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,1),'k')
    rlim([0 r_max])
    title(['Zc=' num2str(Zc_all(iC)) '; Zp=' num2str(Zp_all(iC))])
end
suptitle('Example cells- Black: plaid; Blue: pattern; Red: component')
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_Responsive_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_Responsive_lowSF.fig'))

Zc_use = intersect(resp_ind_all,intersect(find(Zc_all>1.28),find(Zc_all-Zp_all>1.28)));
[x Zc_ind] = sort(Zc_all(Zc_use),'descend');
figure; 
movegui('center')
if length(Zc_use) > 16
    n_use = 16;
else
    n_use = length(Zc_use);
end
[n n2] = subplotn(n_use);
for i = 1:n_use
    subplot(n,n2,i)
    iC = Zc_use(Zc_ind(i));
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1) avg_resp_dir_all_circ(iC,:,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,1),'k')
    rlim([0 r_max])
    title(['Zc=' num2str(Zc_all(iC)) '; Zp=' num2str(Zp_all(iC))])
end
suptitle('Most component-like- Black: plaid; Blue: pattern; Red: component')
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_ComponentLike_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_ComponentLike_lowSF.fig'))

Zp_use = intersect(resp_ind_all,intersect(find(Zp_all>1.28),find(Zp_all-Zc_all>1.28)));
[x Zp_ind] = sort(Zp_all(Zp_use),'descend');
if length(Zp_use) > 16
    n_use = 16;
else
    n_use = length(Zp_use);
end
[n n2] = subplotn(n_use);
figure; 
movegui('center')
for i = 1:n_use
    subplot(n,n2,i)
    iC = Zp_use(Zp_ind(i));
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1) avg_resp_dir_all_circ(iC,:,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,1),'k')
    rlim([0 r_max])
    title(['Zc=' num2str(Zc_all(iC)) '; Zp=' num2str(Zp_all(iC))])
end
suptitle('Most pattern-like- Black: plaid; Blue: pattern; Red: component')
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_PatternLike_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_PatternLike_lowSF.fig'))

ind_h = intersect(resp_ind_all,find(plaid_SI_all>0));
ind_l = intersect(resp_ind_all,find(plaid_SI_all<0));

figure;
movegui('center')
subplot(2,2,1)
scatter(Zc_all(Zc_use),Zp_all(Zc_use))
hold on
scatter(Zc_all(Zp_use),Zp_all(Zp_use))
scatter(Zc_all(setdiff(resp_ind_all,[Zp_use Zc_use])),Zp_all(setdiff(resp_ind_all,[Zp_use Zc_use])),'k')
ylabel('Zp')
xlabel('Zc')
hold on
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
legend({['Zc- ' num2str(length(Zc_use))],['Zp- ' num2str(length(Zp_use))],['N- ' num2str(length(setdiff(resp_ind_all,[Zp_use Zc_use])))]})


subplot(2,2,2)
scatter(Zc_all(ind_l),Zp_all(ind_l),'b')
hold on
scatter(Zc_all(ind_h),Zp_all(ind_h),'r')
plotZcZpBorders
ylabel('Zp')
xlabel('Zc')
xlim([-5 8])
ylim([-5 8])
axis square
legend({['SI<0- ' num2str(length(ind_l))],['SI>0- ' num2str(length(ind_h))]})


subplot(2,2,3)
cdfplot(Zc_all(ind_l))
hold on
cdfplot(Zc_all(ind_h))
xlabel('Zc')
ylabel('Fraction of cells')
xlim([-5 8])
legend({'SI<0','SI>0'},'location','southeast')
[h p] = kstest2(Zc_all(ind_l),Zc_all(ind_h));
title(['p = ' num2str(chop(p,2))])

subplot(2,2,4)
cdfplot(Zp_all(ind_l))
hold on
cdfplot(Zp_all(ind_h))
xlabel('Zp')
ylabel('Fraction of cells')
xlim([-5 8])
[h p] = kstest2(Zp_all(ind_l),Zp_all(ind_h));
title(['p = ' num2str(chop(p,2))])

expts = unique(exptInd);
nexp_area = length(expts);
mouse_list = [{expt.mouse}];
mice = unique(mouse_list(expts));
nmice = length(mice);
suptitle([num2str(nexp_area) ' expts; ' num2str(nmice) ' mice; ' num2str(length(resp_ind_all)) ' cells'])
print(fullfile(summaryDir_F8, 'Figure8_SI_ZcVZp_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_SI_ZcVZp_lowSF.fig'))

%% all cells Pop tuning
close all
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align = nan(nStimDir,nStimDir,2,2);
resp_comp_align = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(resp_ind_all,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,1,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(1)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(2)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(1)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(2)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(1)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align(:,:,ii,1),2), nanmean(resp_dir_align(:,:,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
    hold on
    figure(2)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,ii,1),2); nanmean(resp_dir_align(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(1)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align(:,:,1),2),nanmean(resp_comp_align(:,:,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        title(['All aligned'])
        figure(2)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align(:,:,1),2); nanmean(resp_comp_align(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(1)
suptitle({['Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_lowSF.fig'))

figure(2)
suptitle({['Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_lowSF.fig'))

%Zc only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zc = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zc = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
indZc = cell(nStimDir,2);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(Zc_use,find(max_dir == i));
        indZc{i,ii} = ind; 
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,1,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(3)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(4)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align_Zc(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(3)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(4)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align_Zc(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(3)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,ii,1),2), nanmean(resp_dir_align_Zc(:,:,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    figure(4)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,ii,1),2); nanmean(resp_dir_align_Zc(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(3)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc(:,:,1),2),nanmean(resp_comp_align_Zc(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title(['All aligned'])
        figure(4)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc(:,:,1),2); nanmean(resp_comp_align_Zc(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(3)
suptitle({['Zc cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc_lowSF.fig'))

figure(4)
suptitle({['Zc cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc_lowSF.fig'))

%Zp only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zp = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zp = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
indZp = cell(nStimDir,2);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(Zp_use,find(max_dir == i));
        indZp{i,ii} = ind;
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,1,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(5)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(6)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align_Zp(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(5)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(6)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align_Zp(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(5)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,ii,1),2), nanmean(resp_dir_align_Zp(:,:,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    figure(6)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,ii,1),2); nanmean(resp_dir_align_Zp(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(5)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp(:,:,1),2),nanmean(resp_comp_align_Zp(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title(['All aligned'])
        figure(6)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp(:,:,1),2); nanmean(resp_comp_align_Zp(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(5)
suptitle({['Zp cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp_lowSF.fig'))

figure(6)
suptitle({['Zp cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp_lowSF.fig'))

figure(7)
for ii = 1:2
    subplot(2,3,1)
    errorbar(stimDirs, nanmean(resp_dir_align(:,:,ii,1),2), nanmean(resp_dir_align(:,:,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
    hold on
    subplot(2,3,2)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,ii,1),2), nanmean(resp_dir_align_Zc(:,:,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    subplot(2,3,3)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,ii,1),2), nanmean(resp_dir_align_Zp(:,:,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    subplot(2,3,4)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,ii,1),2); nanmean(resp_dir_align(1,:,ii,1),2)])
    hold on
    subplot(2,3,5)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,ii,1),2); nanmean(resp_dir_align_Zc(1,:,ii,1),2)])
    hold on
    subplot(2,3,6)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,ii,1),2); nanmean(resp_dir_align_Zp(1,:,ii,1),2)])
    hold on
    if ii == 1
        subplot(2,3,1)
        errorbar(stimDirs, nanmean(resp_comp_align(:,:,1),2), nanmean(resp_comp_align(:,:,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        title(['All cells- n = ' num2str(length(resp_ind_all))])
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc(:,:,1),2), nanmean(resp_comp_align_Zc(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title(['Zc cells- n = ' num2str(length(Zc_use))])
        subplot(2,3,3)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp(:,:,1),2), nanmean(resp_comp_align_Zp(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title(['Zp cells- n = ' num2str(length(Zp_use))])
        subplot(2,3,4)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align(:,:,1),2); nanmean(resp_comp_align(1,:,1),2)])
        hold on
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc(:,:,1),2); nanmean(resp_comp_align_Zc(1,:,1),2)])
        hold on
        subplot(2,3,6)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp(:,:,1),2); nanmean(resp_comp_align_Zp(1,:,1),2)])
        hold on
    end
end
suptitle('Aligned Summary: Blue- pattern; Red- component; Yellow- plaid')
print(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar_lowSF.fig'))

%% phase 2 all cells Pop tuning
close all
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align = nan(nStimDir,nStimDir,2,2);
resp_comp_align = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(resp_ind_all,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,2,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(8)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(9)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(8)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(9)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(8)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align(:,:,ii,1),2), nanmean(resp_dir_align(:,:,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
    hold on
    figure(9)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,ii,1),2); nanmean(resp_dir_align(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(8)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align(:,:,1),2),nanmean(resp_comp_align(:,:,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        title(['All aligned'])
        figure(9)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align(:,:,1),2); nanmean(resp_comp_align(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(8)
suptitle({['Phase 2: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_lowSF_phase2.fig'))

figure(9)
suptitle({['Phase2: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_lowSF_phase2.fig'))

%Zc only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zc = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zc = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = indZc{i,ii};
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,2,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(10)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(4)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align_Zc(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(10)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(4)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align_Zc(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(10)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,ii,1),2), nanmean(resp_dir_align_Zc(:,:,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    figure(4)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,ii,1),2); nanmean(resp_dir_align_Zc(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(10)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc(:,:,1),2),nanmean(resp_comp_align_Zc(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title(['All aligned'])
        figure(4)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc(:,:,1),2); nanmean(resp_comp_align_Zc(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(10)
suptitle({['Phase 2 Zc cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc_lowSF_phase2.fig'))

figure(11)
suptitle({['Phase 2 Zc cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc_lowSF_phase2.fig'))

%Zp only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zp = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zp = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = indZp{i,ii};
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,2,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,1,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,1,ii,1) + circshift(avg_resp_dir_all(ind,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(12)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,ii,1), pop_resp_dir(:,i,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(13)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,ii,1); pop_resp_dir(1,i,ii,1)])
        hold on
        resp_dir_align_Zp(:,i,ii,:) = circshift(pop_resp_dir(:,i,ii,:),9-i,1);
        if ii == 1
            figure(12)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(13)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align_Zp(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(12)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,ii,1),2), nanmean(resp_dir_align_Zp(:,:,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    figure(13)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,ii,1),2); nanmean(resp_dir_align_Zp(1,:,ii,1),2)])
    hold on
    if ii == 1
        figure(12)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp(:,:,1),2),nanmean(resp_comp_align_Zp(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title(['All aligned'])
        figure(13)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp(:,:,1),2); nanmean(resp_comp_align_Zp(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(12)
suptitle({['Phase 2 Zp cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp_lowSF_phase2.fig'))

figure(13)
suptitle({['Phase2 Zp cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp_lowSF_phase2.fig'))

figure(14)
for ii = 1:2
    subplot(2,3,1)
    errorbar(stimDirs, nanmean(resp_dir_align(:,:,ii,1),2), nanmean(resp_dir_align(:,:,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
    hold on
    subplot(2,3,2)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,ii,1),2), nanmean(resp_dir_align_Zc(:,:,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    subplot(2,3,3)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,ii,1),2), nanmean(resp_dir_align_Zp(:,:,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    subplot(2,3,4)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,ii,1),2); nanmean(resp_dir_align(1,:,ii,1),2)])
    hold on
    subplot(2,3,5)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,ii,1),2); nanmean(resp_dir_align_Zc(1,:,ii,1),2)])
    hold on
    subplot(2,3,6)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,ii,1),2); nanmean(resp_dir_align_Zp(1,:,ii,1),2)])
    hold on
    if ii == 1
        subplot(2,3,1)
        errorbar(stimDirs, nanmean(resp_comp_align(:,:,1),2), nanmean(resp_comp_align(:,:,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        title(['All cells- n = ' num2str(length(resp_ind_all))])
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc(:,:,1),2), nanmean(resp_comp_align_Zc(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title(['Zc cells- n = ' num2str(length(Zc_use))])
        subplot(2,3,3)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp(:,:,1),2), nanmean(resp_comp_align_Zp(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title(['Zp cells- n = ' num2str(length(Zp_use))])
        subplot(2,3,4)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align(:,:,1),2); nanmean(resp_comp_align(1,:,1),2)])
        hold on
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc(:,:,1),2); nanmean(resp_comp_align_Zc(1,:,1),2)])
        hold on
        subplot(2,3,6)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp(:,:,1),2); nanmean(resp_comp_align_Zp(1,:,1),2)])
        hold on
    end
end
suptitle('Phase 2 Aligned Summary: Blue- pattern; Red- component; Yellow- plaid')
print(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar_lowSF_phase2.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar_lowSF_phase2.fig'))

%% Zc/Zp scatter by phase
Zc_phase2_use = intersect(resp_ind_all,intersect(find(Zc_phase2_all>1.28),find(Zc_phase2_all-Zp_phase2_all>1.28)));
Zp_phase2_use = intersect(resp_ind_all,intersect(find(Zp_phase2_all>1.28),find(Zp_phase2_all-Zc_phase2_all>1.28)));

figure;
subplot(2,2,1)
scatter(Zc_all(resp_ind_all),Zp_all(resp_ind_all),'ok')
hold on
scatter(Zc_all(Zc_phase2_use),Zp_all(Zc_phase2_use),'or')
ylabel('Zp')
xlabel('Zc')
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
title('Phase = 0')
legend({['all cells- n = ' num2str(length(resp_ind_all))],['Zc at 90 deg- n = ' num2str(length(Zc_phase2_use))]},'location','southeast')
axis square
subplot(2,2,2)
scatter(Zc_phase2_all(resp_ind_all),Zp_phase2_all(resp_ind_all),'ok')
hold on
scatter(Zc_phase2_all(Zc_use),Zp_phase2_all(Zc_use),'or')
ylabel('Zp')
xlabel('Zc')
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
title('Phase = 90')
legend({['all cells- n = ' num2str(length(resp_ind_all))],['Zc at 0 deg- n = ' num2str(length(Zc_use))]},'location','southeast')
subplot(2,2,3)
scatter(Zc_all(resp_ind_all),Zp_all(resp_ind_all),'ok')
hold on
scatter(Zc_all(Zp_phase2_use),Zp_all(Zp_phase2_use),'or')
ylabel('Zp')
xlabel('Zc')
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
title('Phase = 0')
legend({['all cells- n = ' num2str(length(resp_ind_all))],['Zp at 90 deg- n = ' num2str(length(Zp_phase2_use))]},'location','southeast')
axis square
subplot(2,2,4)
scatter(Zc_phase2_all(resp_ind_all),Zp_phase2_all(resp_ind_all),'ok')
hold on
scatter(Zc_phase2_all(Zp_use),Zp_phase2_all(Zp_use),'or')
ylabel('Zp')
xlabel('Zc')
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
title('Phase = 90')
legend({['all cells- n = ' num2str(length(resp_ind_all))],['Zp at 0 deg- n = ' num2str(length(Zp_use))]},'location','southeast')
suptitle('Cross-phase Zc/Zp comparison')
print(fullfile(summaryDir_F8, 'Figure8_crossPhaseCompScatter_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_crossPhaseCompScatter_lowSF.fig'))

figure;
errorbarxy(mean(Zc_all(Zc_phase2_use)),mean(Zp_all(Zc_phase2_use)),std(Zc_all(Zc_phase2_use))./sqrt(length(Zc_phase2_use)),std(Zp_all(Zc_phase2_use)./sqrt(length(Zc_phase2_use))),{'bo-', 'b', 'b'})
hold on
errorbarxy(mean(Zc_phase2_all(Zc_use)),mean(Zp_phase2_all(Zc_use)),std(Zc_phase2_all(Zc_use))./sqrt(length(Zc_use)),std(Zp_phase2_all(Zc_use))./sqrt(length(Zc_use)),{'co-', 'c', 'c'})
hold on
errorbarxy(mean(Zc_all(Zp_phase2_use)),mean(Zp_all(Zp_phase2_use)),std(Zc_all(Zp_phase2_use))./sqrt(length(Zp_phase2_use)),std(Zp_all(Zp_phase2_use))./sqrt(length(Zp_phase2_use)),{'ro-', 'r', 'r'})
hold on
errorbarxy(mean(Zc_phase2_all(Zp_use)),mean(Zp_phase2_all(Zp_use)),std(Zc_phase2_all(Zp_use))./sqrt(length(Zp_use)),std(Zp_phase2_all(Zp_use))./sqrt(length(Zp_use)),{'mo-', 'm', 'm'})
legend({'Zc- 90->0','Zc- 0->90','Zp- 90->0','Zp- 0->90'})
ylabel('Zp')
xlabel('Zc')
hold on
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
print(fullfile(summaryDir_F8, 'Figure8_crossPhaseCompSummary_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_crossPhaseCompSummary_lowSF.fig'))

figure;
subplot(2,2,1)
scatter(Zc_all(Zc_use),Zp_all(Zc_use))
hold on
scatter(Zc_all(Zp_use),Zp_all(Zp_use))
scatter(Zc_all(setdiff(resp_ind_all,[Zp_use Zc_use])),Zp_all(setdiff(resp_ind_all,[Zp_use Zc_use])),'k')
ylabel('Zp')
xlabel('Zc')
hold on
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
title('0 deg')
subplot(2,2,3)
scatter(Zc_phase2_all(Zc_use),Zp_phase2_all(Zc_use))
ylabel('Zp')
xlabel('Zc')
hold on
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
title('Zc cells- 90 deg')
subplot(2,2,4)
scatter(Zc_phase2_all(Zp_use),Zp_phase2_all(Zp_use))
ylabel('Zp')
xlabel('Zc')
hold on
plotZcZpBorders
xlim([-5 8])
ylim([-5 8])
axis square
title('Zp cells- 90 deg')


suptitle(['Zc- ' num2str(length(Zc_use)) '; Zp- ' num2str(length(Zp_use)) '; N- ' num2str(length(setdiff(resp_ind_all,[Zp_use Zc_use])))])
print(fullfile(summaryDir_F8, 'Figure8_crossPhaseZcZpScatter_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_crossPhaseZcZpScatter_lowSF.fig'))

