close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F8 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure8');

ds = ['CrossOriRandDir_ExptList'];
eval(ds);
area = 'V1';

nexp = size(expt,2);
totExp = 0;
exptInd = [];
stim_OSI_all = [];
plaid_OSI_all = [];
stim_DSI_all = [];
plaid_DSI_all = [];
Zc_all = [];
Zp_all = [];
plaid_SI_all = [];
h_plaid_SI_all = [];
totCells = 0;
resp_ind_all = [];
resp_ind_dir_all = [];
resp_ind_plaid_all = [];
avg_resp_dir_all = [];
f1_all = [];
f2_all = [];
f2overf1_all = [];
pattern_all = [];
component_all = [];

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

        stim_OSI_all = [stim_OSI_all stim_OSI];
        plaid_OSI_all = [plaid_OSI_all plaid_OSI];
        stim_DSI_all = [stim_DSI_all stim_DSI];
        plaid_DSI_all = [plaid_DSI_all plaid_DSI];
        Zc_all = [Zc_all Zc];
        Zp_all = [Zp_all Zp];
        plaid_SI_all = [plaid_SI_all plaid_SI];
        h_plaid_SI_all = [h_plaid_SI_all h_plaid_SI];
        
        resp_ind = find(sum(sum(h_resp,2),3));
        resp_ind_dir = find(sum(h_resp(:,:,1),2));
        resp_ind_plaid = find(sum(h_resp(:,:,2),2));

        resp_ind_all = [resp_ind_all resp_ind'+totCells];
        resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
        resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];
        
        avg_resp_dir_all = cat(1,avg_resp_dir_all, avg_resp_dir); 
        component_all = [component_all; component];
        pattern_all = [pattern_all; pattern];
        
        if ~isempty(expt(iexp).prFolder)
            ImgFolder = expt(iexp).prFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']))
            f1_all = [f1_all f1];
            f2_all = [f2_all f2];
            f2overf1_all = [f2overf1_all f2overf1];
        else
            f1_all = [f1_all nan(size(stim_OSI))];
            f2_all = [f2_all nan(size(stim_OSI))];
            f2overf1_all = [f2overf1_all nan(size(stim_OSI))];
        end
        
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
     end
end
%%

avg_resp_dir_all_circ = cat(2,avg_resp_dir_all, avg_resp_dir_all(:,1,:,:));
avg_resp_dir_all_circ(find(avg_resp_dir_all_circ<0)) = 0;
avg_resp_dir_all_shift = circshift(avg_resp_dir_all,2,2);
avg_resp_dir_all_circ_shift = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:));
avg_resp_dir_all_circ_shift(find(avg_resp_dir_all_circ_shift<0)) = 0;
component_all_shift = circshift(component_all,2,2);
component_all_circ = cat(2,component_all_shift, component_all_shift(:,1));
component_all_circ(find(component_all_circ<0)) = 0;
stimDirs_circ = [stimDirs stimDirs(1)];


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
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_Responsive.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_Responsive.fig'))

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
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_ComponentLike.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_ComponentLike.fig'))

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
print(fullfile(summaryDir_F8, 'Figure8_ExampleCells_PatternLike.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_ExampleCells_PatternLike.fig'))

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
print(fullfile(summaryDir_F8, 'Figure8_SI_ZcVZp.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_SI_ZcVZp.fig'))

resp_ind_all_f1f2 = intersect(resp_ind_all,find(f1_all>0.02));
ind_h = intersect(resp_ind_all_f1f2,find(f2overf1_all>0.75));
ind_l = intersect(resp_ind_all_f1f2,find(f2overf1_all<0.5));

figure;
movegui('center')
subplot(2,2,1)
histogram(f2overf1_all(resp_ind_all_f1f2))
xlabel('F2/F1')
xlim([0 2])

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
legend({['f2/f1<0.5- ' num2str(length(ind_l))],['f2/f1>0.75- ' num2str(length(ind_h))]})


subplot(2,2,3)
cdfplot(Zc_all(ind_l))
hold on
cdfplot(Zc_all(ind_h))
xlabel('Zc')
ylabel('Fraction of cells')
xlim([-5 8])
legend({'f2/f1<0.5','f2/f1>0.75'},'location','southeast')
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

pr_run = [{expt.prFolder}];
pr_mice = [{expt.mouse}];
expt_pr = intersect(expts,find(~cellfun(@isempty,pr_run)));
nexpt_pr = length(expt_pr);
nmice_pr = length(unique(pr_mice(expt_pr)));

suptitle([num2str(nexpt_pr) ' expts; ' num2str(nmice_pr) ' mice; ' num2str(length(resp_ind_all_f1f2)) ' cells'])
print(fullfile(summaryDir_F8, 'Figure8_f2f1_ZcVZp.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_f2f1_ZcVZp.fig'))

%% all cells Pop tuning
close all
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1),[],2);
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
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,ii,1) + circshift(avg_resp_dir_all(ind,:,ii,1),-nStimDir/4,2),nStimDir/8,2);
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
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar.fig'))

figure(2)
suptitle({['Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar.fig'))

%Zc only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zc = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zc = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(Zc_use,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,ii,1) + circshift(avg_resp_dir_all(ind,:,ii,1),-nStimDir/4,2),nStimDir/8,2);
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
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zc.fig'))

figure(4)
suptitle({['Zc cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zc.fig'))

%Zp only
[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1),[],2);
pop_resp_dir = nan(nStimDir,nStimDir,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align_Zp = nan(nStimDir,nStimDir,2,2);
resp_comp_align_Zp = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(Zp_use,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,ii,1);
            end

            pop_resp_dir(i,:,ii,1) = mean(temp_resp_dir(ind,:),1);
            pop_resp_dir(i,:,ii,2) = std(temp_resp_dir(ind,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,ii,1) + circshift(avg_resp_dir_all(ind,:,ii,1),-nStimDir/4,2),nStimDir/8,2);
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
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_errorbar_Zp.fig'))

figure(6)
suptitle({['Zp cells: Blue- pattern; Red- component; Yellow- plaid'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_polar_Zp.fig'))

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
print(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuningAligned_errorbar.fig'))

%% OSI, DSI, F2/F1 and Zp/Zc

figure; 
movegui('center')
subplot(3,2,1)
scatter(stim_OSI_all(resp_ind_all), Zc_all(resp_ind_all));
xlabel('OSI')
ylabel('Zc')
ylim([-10 10])
subplot(3,2,2)
scatter(stim_OSI_all(resp_ind_all), Zp_all(resp_ind_all));
xlabel('OSI')
ylabel('Zp')
ylim([-10 10])
subplot(3,2,3)
scatter(stim_DSI_all(resp_ind_all), Zc_all(resp_ind_all));
xlabel('DSI')
ylabel('Zc')
ylim([-10 10])
subplot(3,2,4)
scatter(stim_DSI_all(resp_ind_all), Zp_all(resp_ind_all));
xlabel('DSI')
ylabel('Zp')
ylim([-10 10])
subplot(3,2,5)
scatter(log10(f2overf1_all(resp_ind_all)), Zc_all(resp_ind_all));
xlabel('log(F2/F1)')
ylabel('Zc')
ylim([-10 10])
subplot(3,2,6)
scatter(log10(f2overf1_all(resp_ind_all)), Zp_all(resp_ind_all));
xlabel('log(F2/F1)')
ylabel('Zp')
ylim([-10 10])

