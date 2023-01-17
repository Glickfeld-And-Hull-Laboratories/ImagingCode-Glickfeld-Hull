close all
clear all
clc
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'Adaptation', 'SFSummary');

ds = 'AdaptSF_ExptList';
eval(ds);
nexp = size(expt,2);
totCells = 0;
norm_sf = [];
norm_prefsf = [];
norm_dir = [];
norm_prefdir = [];

for iexp = 1:nexp
    if expt(iexp).randDir
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).adaptFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(ImgFolder, nrun);
        
        fprintf([mouse ' ' date '\n'])
        
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    
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

figure(1)
subplot(2,2,1)
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


subplot(2,2,2)
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
norm_prefsf_avg = zeros(length(oris),2);
for is = 1:length(oris)
    ind = find(norm_prefdir(:,3) == oris(is));
    norm_prefori_avg(is,:) = [mean(norm_prefdir(ind,1),1) std(norm_prefdir(ind,1),[],1)./sqrt(length(ind))];
end
errorbar(oris,norm_prefori_avg(:,1),norm_prefori_avg(:,2),'or')
[h_pref p_pref stats_pref] = anova1(norm_prefdir(:,1), norm_prefdir(:,3),'off');
table_pref = multcompare(stats_pref,'Display','off');
ind = find(table_pref(:,end)<0.05);
groups_pref = mat2cell(oris(table_pref(ind,1:2)),ones(1,size(table_all(ind,:),1)));
stats_pref = table_pref(ind,end);

sigstar(groups_pref,stats_pref);

print(fullfile(summaryDir,'adaptationBySF.pdf'),'-dpdf','-bestfit');

uniqueCells = size(norm_prefsf,1);
nboot = 1000;
match_cell_diff = cell(1,nboot);
nonmatch_cell_diff = cell(1,nboot);
for iboot = 1:nboot
    match_cell_diff{iboot} = [];
    nonmatch_cell_diff{iboot} = [];
    for iC = 1:uniqueCells
        iCell = norm_prefsf(iC,4);
        ind = intersect(find(norm_sf(:,2) ~= norm_prefsf(iC,2)),find(norm_sf(:,4) == iCell));
        if exist('ind')
            if iboot == 1;
                for i = 1:length(ind)
                    match_cell_diff{iboot} = [match_cell_diff{iboot} abs(norm_sf(ind(i),1)-norm_prefsf(iC,1))];
                end
            end
            ind_rand = randi(uniqueCells,length(ind));
            for i = 1:length(ind_rand)
                nonmatch_cell_diff{iboot} = [nonmatch_cell_diff{iboot} abs(norm_sf(ind_rand(i),1)-norm_prefsf(iC,1))];
            end
        end
    end
end

figure(2);
H1 = cdfplot(match_cell_diff{1});
YData = unique(H1.YData);
np = length(YData);
XData = zeros(nboot,np);
hold on
for iboot = 1:nboot
    H2 = cdfplot(nonmatch_cell_diff{iboot});
    set(H2,'Color',[0.5 0.5 0.5])
    [tempy ind] = unique(H2.YData);
    tempx = H2.XData(ind);
    XData(iboot,:) = interp1(tempy,tempx,YData); 
end
hold on
plot(mean(XData,1),YData,'r')
xlabel('Absolute diff in adaptation from pref SF')
ylabel('Fraction of cells')
title([num2str(uniqueCells) ' cells; ' num2str(nexp) ' mice'])
print(fullfile(summaryDir,'adaptationDiffFromPrefSF.pdf'),'-dpdf','-bestfit','-painters');