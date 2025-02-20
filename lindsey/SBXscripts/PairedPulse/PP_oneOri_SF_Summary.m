close all
clear all
clc
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'Adaptation', 'SFSummary');

ds = 'AdaptSF_ExptList';
eval(ds);
nexp = size(expt,2);

area = 'V1_GC6s_WT';

expt_use = 26;
totCells = 0;
norm_sf = [];
norm_prefsf = [];
resp_pref = [];
resp_var = [];
ind_use = [];
mice = [];
for iexp = expt_use
    if ~expt(iexp).randDir
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).adaptFolder;
        nrun = length(ImgFolder);
        run_str = catRunName(ImgFolder, nrun);
        
        fprintf([mouse ' ' date '\n'])
        
        if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
        else
            fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat'])
        end
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    
        norm_sf_temp =  [reshape(norm_dfof_stim_nan, [nCells.*nSF 1]) reshape(repmat(sfs,[nCells 1]), [nCells.*nSF 1])];
        norm_prefsf_temp =  [norm_dfof_stim_pref sfs(pref_sf)'];
        resp_ind = find(sum(h1_ori,[2 3]));
        resp_dfof_pref = indOnly(resp_dfof_stim(:,:,1),pref_sf);
        %resp_var_pref = indOnly(resp_dfof_var(:,:,1),pref_sf);
        resp_pref_temp = [resp_dfof_pref sfs(pref_sf)'];
        %resp_var_temp = [resp_var_pref sfs(pref_sf)'];
        ind = intersect(resp_ind,find(resp_dfof_pref>0.05));
        ind_use = [ind_use; ind+totCells];
        norm_sf = [norm_sf; norm_sf_temp];
        norm_prefsf = [norm_prefsf; norm_prefsf_temp];
        resp_pref = [resp_pref; resp_pref_temp];
        %resp_var = [resp_var; resp_var_temp];
    
        nCells = size(pref_sf,1);
        totCells = totCells+nCells;
        mice = [mice;{mouse}];
    end
end

%%
save(fullfile(summaryDir,['adaptationBySF_' area '.mat']),'norm_sf', 'norm_prefsf', 'resp_pref', 'resp_var', 'ind_use', 'sfs')


figure(1)
subplot(3,2,1)
swarmchart(norm_sf(:,2), norm_sf(:,1),'k')
set(gca,'XScale','log')
ylabel('Normalized dF/F')
xlabel('Spatial frequency')
ylim([0 4])
xlim([0.02 1])
hline(1)
hold on
title('Responsive')
norm_sf_avg = zeros(length(sfs),2);
for is = 1:length(sfs)
    ind = find(norm_sf(:,2) == sfs(is));
    norm_sf_avg(is,:) = [mean(norm_sf(ind,1),1,'omitnan') std(norm_sf(ind,1),[],1,'omitnan')./sqrt(length(ind))];
end
errorbar(sfs,norm_sf_avg(:,1),norm_sf_avg(:,2),'or')

[h_all p_all stats_all] = anova1(norm_sf(:,1), norm_sf(:,2),'off');
table_all = multcompare(stats_all,'Display','off');
ind = find(table_all(:,end)<0.05);
groups_all = mat2cell(sfs(table_all(ind,1:2)),ones(1,size(table_all(ind,:),1)));
stats_all = table_all(ind,end);

%sigstar(groups_all,stats_all);


subplot(3,2,2)
swarmchart(norm_prefsf(ind_use,2), norm_prefsf(ind_use,1),'k')
set(gca,'XScale','log')
ylabel('Normalized dF/F')
xlabel('Spatial frequency')
hold on
ylim([0 4])
xlim([0.02 1])
hline(1)
title('Preferred')
norm_prefsf_avg = zeros(length(sfs),2);
for is = 1:length(sfs)
    ind = intersect(ind_use,find(norm_prefsf(:,2) == sfs(is)));
    norm_prefsf_avg(is,:) = [mean(norm_prefsf(ind,1),1,'omitnan') std(norm_prefsf(ind,1),[],1,'omitnan')./sqrt(length(ind))];
end
errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'or')
[h_pref p_pref stats_pref] = anova1(norm_prefsf(ind_use,1), norm_prefsf(ind_use,2),'off');
table_pref = multcompare(stats_pref,'Display','off');
ind = find(table_pref(:,end)<0.05);
groups_pref = mat2cell(sfs(table_pref(ind,1:2)),ones(1,size(table_pref(ind,:),1)));
stats_pref = table_pref(ind,end);

subplot(3,2,3)
swarmchart(resp_pref(ind_use,2), resp_pref(ind_use,1),'k')
set(gca,'XScale','log')
ylabel('Resp 1 dF/F')
xlabel('Spatial frequency')
hold on
ylim([0 2])
xlim([0.02 1])
title('Preferred')
resp_pref_avg = zeros(length(sfs),2);
for is = 1:length(sfs)
    ind = intersect(ind_use,find(resp_pref(:,2) == sfs(is)));
    resp_pref_avg(is,:) = [mean(resp_pref(ind,1),1,'omitnan') std(resp_pref(ind,1),[],1,'omitnan')./sqrt(length(ind))];
end
errorbar(sfs,resp_pref_avg(:,1),resp_pref_avg(:,2),'or')
[h_resp p_resp stats_resp] = anova1(resp_pref(ind_use,1), resp_pref(ind_use,2),'off');
table_resp = multcompare(stats_resp,'Display','off');
ind = find(table_resp(:,end)<0.05);
groups_resp = mat2cell(sfs(table_resp(ind,1:2)),ones(1,size(table_resp(ind,:),1)));
stats_resp = table_resp(ind,end);

subplot(3,2,4)
errorbar(resp_pref_avg(:,1),norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),norm_prefsf_avg(:,2),resp_pref_avg(:,2),resp_pref_avg(:,2),'o')
xlabel('Resp 1 dF/F')
ylabel('Normalized dF/F')
xlim([0 1])
ylim([0 1.5])

% subplot(3,2,5)
% swarmchart(resp_var(ind_use,2), resp_var(ind_use,1)./(resp_pref(ind_use,1).^2),'k')
% set(gca,'XScale','log')
% ylabel('Resp 1 variance/mean^2')
% xlabel('Spatial frequency')
% hold on
% ylim([0 10])
% xlim([0.02 1])
% title('Preferred')
% resp_snr_avg = zeros(length(sfs),2);
% for is = 1:length(sfs)
%     ind = intersect(ind_use,find(resp_var(:,2) == sfs(is)));
%     resp_snr_avg(is,:) = [mean(resp_var(ind,1)./(resp_pref(ind,1).^2),1,'omitnan') std(resp_var(ind,1)./(resp_pref(ind,1).^2),[],1,'omitnan')./sqrt(length(ind))];
% end
% errorbar(sfs,resp_snr_avg(:,1),resp_snr_avg(:,2),'or')

% resp_var_avg = zeros(length(sfs),2);
% for is = 1:length(sfs)
%     ind = intersect(ind_use,find(resp_var(:,2) == sfs(is)));
%     resp_var_avg(is,:) = [mean(resp_var(ind,1),1,'omitnan') std(resp_var(ind,1),[],1,'omitnan')./sqrt(length(ind))];
% end
% errorbar(resp_var_avg(:,1),norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),norm_prefsf_avg(:,2),resp_var_avg(:,2),resp_var_avg(:,2),'o')
% xlabel('Resp 1 variance')
% ylabel('Normalized dF/F')
% xlim([0 0.2])
% ylim([0 1.5])

% subplot(3,2,6)
% errorbar(resp_snr_avg(:,1),norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),norm_prefsf_avg(:,2),resp_snr_avg(:,2),resp_snr_avg(:,2),'o')
% xlabel('Resp 1 variance/mean^2')
% ylabel('Normalized dF/F')
% xlim([0 5])
% ylim([0 1.5])
% 
% lm = fitlm(resp_snr_avg(:,1), norm_prefsf_avg(:,1));

%sigstar(groups_pref,stats_pref);
subplot(3,2,6)
histogram(norm_prefsf(ind_use,2),'Normalization','probability')
hold on
cdfplot(norm_prefsf(ind_use,2))
sgtitle([area ' ' num2str(size(unique(mice),1)) ' mice; ' num2str(length(ind_use)) ' cells'])
print(fullfile(summaryDir,['adaptationBySFonly_' area '.pdf']),'-dpdf','-bestfit');

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

figure; 
subplot(1,2,1)
errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'or')
set(gca,'XScale','log')
ylabel('Normalized dF/F')
xlabel('Pref. Spatial frequency')
ylim([0 1.4])
xlim([0.02 1])
hline(1)

subplot(1,2,2)
errorbar(sfs,resp_pref_avg(:,1),resp_pref_avg(:,2),'or')
set(gca,'XScale','log')
ylabel('R1 dF/F')
xlabel('Pref. Spatial frequency')
ylim([0 0.4])
xlim([0.02 1])

sgtitle(mice')

print(fullfile(summaryDir,['adaptationBySFAvgonly_' area '.pdf']),'-dpdf','-bestfit');