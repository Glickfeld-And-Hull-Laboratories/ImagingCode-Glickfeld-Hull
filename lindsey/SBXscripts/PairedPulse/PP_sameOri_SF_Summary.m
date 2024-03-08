close all
clear all
clc
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'Adaptation', 'SFSummary');

ds = 'AdaptSF_ExptList';
eval(ds);
nexp = size(expt,2);

area = 'LM';
area_ind = find(strcmp([expt.img_loc], area));
randDir = 1;
nDir = 6;
dir_ind = intersect(find([expt.randDir]==randDir),find([expt.nDir]==nDir));

expt_use = intersect(area_ind,dir_ind);

totCells = 0;
resp_dfof_stim_all = [];
h1_ori_all = [];
norm_sf = [];
norm_prefsf = [];
norm_dir = [];
norm_prefdir = [];

for iexp = expt_use
    if expt(iexp).randDir
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).adaptFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(ImgFolder, nrun);
        
        fprintf([mouse ' ' date '\n'])
        
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    
        resp_dfof_stim_all = cat(1, resp_dfof_stim_all, squeeze(resp_dfof_stim));
        h1_ori_all = cat(1, h1_ori_all, h1_ori);

        norm_sf_all(:,end) =  norm_sf_all(:,end) + totCells;
        norm_prefsf_all(:,end) =  norm_prefsf_all(:,end) + totCells;
        norm_dir_all(:,end) =  norm_dir_all(:,end) + totCells;
        norm_prefdir_all(:,end) =  norm_prefdir_all(:,end) + totCells;
    
        norm_sf = [norm_sf; norm_sf_all];
        norm_prefsf = [norm_prefsf; norm_prefsf_all];
        norm_dir = [norm_dir; norm_dir_all];
        norm_prefdir = [norm_prefdir; norm_prefdir_all];
    
        nCells = size(pref_sf,1);
        totCells = totCells+nCells;
    end
end
sfs = unique(norm_sf(:,2));
oris = unique(norm_sf(:,3));

resp_ind = find(sum(sum(h1_ori_all,2),3));

[max_val max_ori] = max(resp_dfof_stim_all(:,:,:,1),[],2);
max_all =  max(max(resp_dfof_stim_all(:,:,:,1),[],2),[],3);
resp_dfof_stim_all_shift = resp_dfof_stim_all;
resp_dfof_stim_all_norm = resp_dfof_stim_all;
resp_dfof_stim_all_norm_all = resp_dfof_stim_all;
for iC = 1:totCells
    for iSF = 1:length(sfs)
        resp_dfof_stim_all_shift(iC,:,iSF,:) = circshift(resp_dfof_stim_all(iC,:,iSF,:),5-max_ori(iC,1,iSF),2);
        resp_dfof_stim_all_norm(iC,:,iSF,:) = resp_dfof_stim_all_shift(iC,:,iSF,:)./max_val(iC,1,iSF);
    end
    resp_dfof_stim_all_norm_all(iC,:,:,:) = resp_dfof_stim_all_shift(iC,:,:,:)./max_all(iC,:);
end

figure
subplot(2,3,1)
for iSF = 1:length(sfs)
    ind_use = find(sum(h1_ori_all(:,:,iSF),2));
    errorbar(oris, squeeze(mean(resp_dfof_stim_all_norm(ind_use,:,iSF,1),1)), squeeze(std(resp_dfof_stim_all_norm(ind_use,:,iSF,1),[],1)./sqrt(length(ind_use))));
    hold on
    str{iSF} = [num2str(sfs(iSF)) 'cpd; n = ' num2str(length(ind_use))];
end
legend(str)
ylim([0 1.5])
xlabel('Orientation (deg)')
ylabel('Normalized dF/F')
subplot(2,3,2)
for iSF = 1:length(sfs)
    ind_use = find(sum(h1_ori_all(:,:,iSF),2));
    cdfplot(max_val(ind_use,iSF))
    hold on
end
xlabel('max dF/F')
OSI = (resp_dfof_stim_all_shift(:,5,:,1)-resp_dfof_stim_all_shift(:,1,:,1))./(resp_dfof_stim_all_shift(:,5,:,1)+resp_dfof_stim_all_shift(:,1,:,1));
subplot(2,3,3)
for iSF = 1:length(sfs)
    ind_use = find(sum(h1_ori_all(:,:,iSF),2));
    cdfplot(OSI(ind_use,iSF))
    hold on
end
xlabel('OSI')
h1_all = squeeze(sum(h1_ori_all,2)>0);
ind_match = find(sum(h1_all,2)==3);
subplot(2,3,4)
for iSF = 1:length(sfs)
    errorbar(oris, squeeze(mean(resp_dfof_stim_all_norm(ind_match,:,iSF,1),1)), squeeze(std(resp_dfof_stim_all_norm(ind_match,:,iSF,1),[],1)./sqrt(length(ind_match))));
    hold on
    str{iSF} = [num2str(sfs(iSF)) 'cpd; n = ' num2str(length(ind_match))];
end
legend(str)
xlabel('Orientation (deg)')
ylabel('Normalized dF/F')
ylim([0 1.5])
subplot(2,3,5)
for iSF = 1:length(sfs)
    cdfplot(max_val(ind_match,iSF))
    hold on
end
xlabel('max dF/F')
OSI = (resp_dfof_stim_all_shift(:,5,:,1)-resp_dfof_stim_all_shift(:,1,:,1))./(resp_dfof_stim_all_shift(:,5,:,1)+resp_dfof_stim_all_shift(:,1,:,1));
subplot(2,3,6)
for iSF = 1:length(sfs)
    cdfplot(OSI(ind_match,iSF))
    hold on
end
xlabel('OSI')
suptitle(['LM- ' [expt(expt_use).mouse]])
print(fullfile(summaryDir,['tuningBySF_' area '_' num2str(nDir) 'Dirs.pdf']),'-dpdf','-bestfit');


figure; 
subplot(2,3,1)
[n edges bin] = histcounts(norm_prefsf(:,2),[sfs; 1]);
histogram(bin)
set(gca,'XtickLabels', sfs(unique(bin)))
ylabel('Number of cells')
xlabel('Spatial frequency (cpd)')
title([num2str(size(norm_prefsf,1)) ' cells'])

%figure;
subplot(2,3,2)
swarmchart(norm_sf(:,2), norm_sf(:,1),'k')
set(gca,'XScale','log')
xticks(sfs)
ylabel('Normalized dF/F')
xlabel('Spatial frequency')
ylim([0 6])
xlim([0.03 1])
hline(1)
hold on
title('Responsive')
norm_sf_avg = zeros(length(sfs),2);
for is = 1:length(sfs)
    ind = find(norm_sf(:,2) == sfs(is));
    norm_sf_avg(is,:) = [mean(norm_sf(ind,1),1) std(norm_sf(ind,1),[],1)./sqrt(length(ind))];
end
errorbar(sfs,norm_sf_avg(:,1),norm_sf_avg(:,2),'or')

[h_all p_all stats_all] = anova1(norm_sf(:,1), norm_sf(:,2),'off');
table_all = multcompare(stats_all,'Display','off');
ind = find(table_all(:,end)<0.05);
groups_all = mat2cell(sfs(table_all(ind,1:2)),ones(1,size(table_all(ind,:),1)));
stats_all = table_all(ind,end);

sigstar(groups_all,stats_all);


subplot(2,3,3)
swarmchart(norm_prefsf(:,2), norm_prefsf(:,1),'k')
set(gca,'XScale','log')
xticks(sfs)
ylabel('Normalized dF/F')
xlabel('Spatial frequency')
hold on
ylim([0 6])
xlim([0.03 1])
hline(1)
title('Preferred')
norm_prefsf_avg = zeros(length(sfs),2);
for is = 1:length(sfs)
    ind = find(norm_prefsf(:,2) == sfs(is));
    norm_prefsf_avg(is,:) = [mean(norm_prefsf(ind,1),1) std(norm_prefsf(ind,1),[],1)./sqrt(length(ind))];
end
errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'or')
[h_pref p_pref stats_pref] = anova1(norm_prefsf(:,1), norm_prefsf(:,2),'off');
table_pref = multcompare(stats_pref,'Display','off');
ind = find(table_pref(:,end)<0.05);
groups_pref = mat2cell(sfs(table_pref(ind,1:2)),ones(1,size(table_all(ind,:),1)));
stats_pref = table_pref(ind,end);

sigstar(groups_pref,stats_pref);

subplot(2,2,3)
swarmchart(norm_dir(:,3), norm_dir(:,1),'k')
xticks(oris)
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
ylim([0 6])
xlim([-10 160])
hline(1)
hold on
title('Responsive')
norm_ori_avg = zeros(length(oris),2);
for is = 1:length(oris)
    ind = find(norm_dir(:,3) == oris(is));
    norm_ori_avg(is,:) = [mean(norm_dir(ind,1),1) std(norm_dir(ind,1),[],1)./sqrt(length(ind))];
end
errorbar(oris,norm_ori_avg(:,1),norm_ori_avg(:,2),'or')
[h_all p_all stats_all] = anova1(norm_dir(:,1), norm_dir(:,3),'off');
table_all = multcompare(stats_all,'Display','off');
ind = find(table_all(:,end)<0.05);
groups_all = mat2cell(oris(table_all(ind,1:2)),ones(1,size(table_all(ind,:),1)));
stats_all = table_all(ind,end);

sigstar(groups_all,stats_all);

subplot(2,2,4)
swarmchart(norm_prefdir(:,3), norm_prefdir(:,1),'k')
xticks(oris)
ylabel('Normalized dF/F')
xlabel('Orientation (deg)')
hold on
ylim([0 6])
xlim([-10 160])
hline(1)
title('Preferred')
norm_prefori_avg = zeros(length(oris),2);
for is = 1:length(oris)
    ind = find(norm_prefdir(:,3) == oris(is));
    norm_prefori_avg(is,:) = [mean(norm_prefdir(ind,1),1) std(norm_prefdir(ind,1),[],1)./sqrt(length(ind))];
end
errorbar(oris,norm_prefori_avg(:,1),norm_prefori_avg(:,2),'or')
[h_pref p_pref stats_pref] = anova1(norm_prefdir(:,1), norm_prefdir(:,3),'off');
table_pref = multcompare(stats_pref,'Display','off');
ind = find(table_pref(:,end)<0.05);
if length(ind)>1
    groups_pref = mat2cell(oris(table_pref(ind,1:2)),ones(1,size(table_all(ind,:),1)));
else
    groups_pref = mat2cell(oris(table_pref(ind,1:2))',ones(1,size(table_all(ind,:),1)));
end
stats_pref = table_pref(ind,end);

sigstar(groups_pref,stats_pref);
suptitle(area)
print(fullfile(summaryDir,['adaptationBySF_' area '_' num2str(nDir) 'Dirs.pdf']),'-dpdf','-bestfit');
save(fullfile(summaryDir,['adaptationDirBySF_' area '_' num2str(nDir) 'Dirs.mat']),'oris', 'sfs', 'resp_dfof_stim_all', 'resp_dfof_stim_all_norm', 'h1_ori_all', 'norm_sf', 'norm_prefsf', 'norm_dir', 'norm_prefdir', 'norm_sf_avg', 'norm_prefsf_avg', 'norm_ori_avg', 'norm_prefori_avg')
% uniqueCells = size(norm_prefsf,1);
% nboot = 1000;
% match_cell_diff = cell(1,nboot);
% nonmatch_cell_diff = cell(1,nboot);
% for iboot = 1:nboot
%     match_cell_diff{iboot} = [];
%     nonmatch_cell_diff{iboot} = [];
%     for iC = 1:uniqueCells
%         iCell = norm_prefsf(iC,4);
%         ind = intersect(find(norm_sf(:,2) ~= norm_prefsf(iC,2)),find(norm_sf(:,4) == iCell));
%         if exist('ind')
%             if iboot == 1;
%                 for i = 1:length(ind)
%                     match_cell_diff{iboot} = [match_cell_diff{iboot} abs(norm_sf(ind(i),1)-norm_prefsf(iC,1))];
%                 end
%             end
%             ind_rand = randi(uniqueCells,length(ind));
%             for i = 1:length(ind_rand)
%                 nonmatch_cell_diff{iboot} = [nonmatch_cell_diff{iboot} abs(norm_sf(ind_rand(i),1)-norm_prefsf(iC,1))];
%             end
%         end
%     end
% end

% figure(2);
% H1 = cdfplot(match_cell_diff{1});
% YData = unique(H1.YData);
% np = length(YData);
% XData = zeros(nboot,np);
% hold on
% for iboot = 1:nboot
%     H2 = cdfplot(nonmatch_cell_diff{iboot});
%     set(H2,'Color',[0.5 0.5 0.5])
%     [tempy ind] = unique(H2.YData);
%     tempx = H2.XData(ind);
%     XData(iboot,:) = interp1(tempy,tempx,YData); 
% end
% hold on
% plot(mean(XData,1),YData,'r')
% xlabel('Absolute diff in adaptation from pref SF')
% ylabel('Fraction of cells')
% title([area ': ' num2str(uniqueCells) ' cells; ' num2str(nexp) ' mice'])
% print(fullfile(summaryDir,['adaptationDiffFromPrefSF_' area '.pdf']),'-dpdf','-bestfit','-painters');