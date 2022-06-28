close all
clear all
clc
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'Adaptation', 'SFSummary');

ds = 'AdaptPhaseRev_ExptList';
eval(ds);
nexp = size(expt,2);
totCells = 0;
norm_prefsf = [];
f2overf1_all = [];
pref_sf_all = [];
ind_use = [];
mice = [];
for iexp = 1:nexp
    if expt(iexp).prSF == 0.05
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).adaptFolder;
        nrun = length(ImgFolder);
        ad_run_str = catRunName(ImgFolder, nrun);
        ImgFolder = expt(iexp).prFolder;
        nrun = length(ImgFolder);
        pr_run_str = catRunName(ImgFolder, nrun);
        
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_input.mat']))
        
        if ~input.doRandDir
            fprintf([mouse ' ' date '\n'])
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_dfofData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ad_run_str], [date '_' mouse '_' ad_run_str '_stimData.mat']))     
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' pr_run_str], [date '_' mouse '_' pr_run_str '_f1f2.mat'])) 
        
        
    
            norm_prefsf_temp =  [norm_dfof_stim_pref sfs(pref_sf)'];
            resp_ind = find(sum(h1_ori,[2 3]));
            resp_dfof_pref = indOnly(resp_dfof_stim,pref_sf);
            ind = intersect(resp_ind,find(resp_dfof_pref>0.05));
            resp_pr = intersect(ind, find(f1_dir>0.03));
            ind_use = [ind_use; resp_pr+totCells];
            norm_prefsf = [norm_prefsf; norm_prefsf_temp];
            f2overf1_all = [f2overf1_all; f2overf1];
            pref_sf_all = [pref_sf_all pref_sf];
        
            nCells = size(pref_sf,1);
            totCells = totCells+nCells;
            mice = [mice;{mouse}];
        end
    end
end

%%
figure(1)
scatter(f2overf1_all(ind_use),norm_prefsf(ind_use),[],sfs(pref_sf_all(ind_use)))
xlabel('F2/F1')
ylabel('Norm resp')
colorbar
sgtitle([num2str(size(unique(mice),1)) ' mice; ' num2str(length(ind_use)) ' cells'])
print(fullfile(summaryDir,'adaptationByF2overF1.pdf'),'-dpdf','-bestfit');           

% subplot(1,2,1)
% swarmchart(norm_sf(:,2), norm_sf(:,1),'k')
% set(gca,'XScale','log')
% ylabel('Normalized dF/F')
% xlabel('Spatial frequency')
% ylim([0 6])
% xlim([0.02 1])
% hline(1)
% hold on
% title('Responsive')
% norm_sf_avg = zeros(length(sfs),2);
% for is = 1:length(sfs)
%     ind = find(norm_sf(:,2) == sfs(is));
%     norm_sf_avg(is,:) = [mean(norm_sf(ind,1),1,'omitnan') std(norm_sf(ind,1),[],1,'omitnan')./sqrt(length(ind))];
% end
% errorbar(sfs,norm_sf_avg(:,1),norm_sf_avg(:,2),'or')
% 
% [h_all p_all stats_all] = anova1(norm_sf(:,1), norm_sf(:,2),'off');
% table_all = multcompare(stats_all,'Display','off');
% ind = find(table_all(:,end)<0.05);
% groups_all = mat2cell(sfs(table_all(ind,1:2)),ones(1,size(table_all(ind,:),1)));
% stats_all = table_all(ind,end);
% 
% %sigstar(groups_all,stats_all);
% 
% 
% subplot(1,2,2)
% swarmchart(norm_prefsf(ind_use,2), norm_prefsf(ind_use,1),'k')
% set(gca,'XScale','log')
% ylabel('Normalized dF/F')
% xlabel('Spatial frequency')
% hold on
% ylim([0 6])
% xlim([0.02 1])
% hline(1)
% title('Preferred')
% norm_prefsf_avg = zeros(length(sfs),2);
% for is = 1:length(sfs)
%     ind = intersect(ind_use,find(norm_prefsf(:,2) == sfs(is)));
%     norm_prefsf_avg(is,:) = [mean(norm_prefsf(ind,1),1,'omitnan') std(norm_prefsf(ind,1),[],1,'omitnan')./sqrt(length(ind))];
% end
% errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'or')
% [h_pref p_pref stats_pref] = anova1(norm_prefsf(ind_use,1), norm_prefsf(ind_use,2),'off');
% table_pref = multcompare(stats_pref,'Display','off');
% ind = find(table_pref(:,end)<0.05);
% groups_pref = mat2cell(sfs(table_pref(ind,1:2)),ones(1,size(table_all(ind,:),1)));
% stats_pref = table_pref(ind,end);
% 
% %sigstar(groups_pref,stats_pref);
% sgtitle([num2str(size(unique(mice),1)) ' mice; ' num2str(length(ind_use)) ' cells'])
% print(fullfile(summaryDir,'adaptationBySFonly.pdf'),'-dpdf','-bestfit');

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
% 
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
% title([num2str(uniqueCells) ' cells; ' num2str(nexp) ' mice'])
% print(fullfile(summaryDir,'adaptationDiffFromPrefSF.pdf'),'-dpdf','-bestfit','-painters');