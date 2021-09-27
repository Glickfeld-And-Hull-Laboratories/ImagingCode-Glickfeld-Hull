close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F6 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure6');

ds = ['CrossOriSingleStimAdapt_ExptList'];
eval(ds);

nexp = size(expt,2);
totCells = 0;
resptest_ind_all = [];
respmask_ind_all = [];
respplaid_ind_all = [];
preftest_ind_all = [];
prefmask_ind_all = [];
respmask_ad_ind_all = [];
noadapt_resp_tc_all = [];
singadapt_resp_tc_all = [];
noadapt_resp_cell_all= cell(5,5);
singadapt_resp_cell_all= cell(5,5);
expt_ind = [];
h_MI_all = [];
p_MI_all = [];
h_MI_ad_all = [];
p_MI_ad_all = [];

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    run_str = ['runs-' cell2mat(expt(iexp).coFolder)];
    dir_run_str = ['runs-' cell2mat(expt(iexp).dirFolder)];
    
    fprintf([mouse ' ' date '\n'])
    fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
    dir_fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str]);
    load(fullfile(fn,[date '_' mouse '_' run_str '_allCellResp.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_dataStim.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_TCs.mat']));
    nMask = length(maskCons);
    nTest = length(testCons);
    
    expt_ind = [expt_ind iexp.*ones(1,nCells)];
    for im = 1:nMask
        for it = 1:nTest
            noadapt_resp_cell_all{im,it} = [noadapt_resp_cell_all{im,it}; nanmean(noadapt_resp_cell{im,it},2)];
            singadapt_resp_cell_all{im,it} = [singadapt_resp_cell_all{im,it}; nanmean(singadapt_resp_cell{im,it},2)];
        end
    end
    test_resp = mean(noadapt_resp_cell_all{1,end},2);
    test_resp(find(test_resp<0)) = 0;
    mask_resp = mean(noadapt_resp_cell_all{end,1},2);
    mask_resp(find(mask_resp<0)) = 0;
    
    ad_test_resp = mean(singadapt_resp_cell_all{1,end},2);
    ad_test_resp(find(ad_test_resp<0)) = 0;
    ad_mask_resp = mean(singadapt_resp_cell_all{end,1},2);
    ad_mask_resp(find(ad_mask_resp<0)) = 0;
    
    h_resptest = ttest2(noadapt_resp_cell{1,end},noadapt_resp_cell{1,1},'dim',2, 'tail','right');
    h_respmask = ttest2(noadapt_resp_cell{end,1},noadapt_resp_cell{1,1},'dim',2, 'tail','right');
    h_respplaid = ttest2(noadapt_resp_cell{end,end},noadapt_resp_cell{1,1},'dim',2, 'tail','right');
    h_respmask_ad = ttest2(singadapt_resp_cell{end,1},singadapt_resp_cell{1,1},'dim',2, 'tail','right');
    
    h_preftest = ttest2(noadapt_resp_cell{1,end},noadapt_resp_cell{end,1},'dim',2, 'tail','right');
    h_prefmask = ttest2(noadapt_resp_cell{end,1},noadapt_resp_cell{1,end},'dim',2, 'tail','right');
    
    p_MI = zeros(nCells,1);
    h_MI = zeros(nCells,1);
    p_MI_ad = zeros(nCells,1);
    h_MI_ad = zeros(nCells,1);
    for iCell = 1:nCells
        [h_MI(iCell) p_MI(iCell)] = ttest(noadapt_resp_cell{end,end}(iCell,:),test_resp(iCell)+mask_resp(iCell));
        [h_MI_ad(iCell) p_MI_ad(iCell)] = ttest(singadapt_resp_cell{end,end}(iCell,:),ad_test_resp(iCell)+ad_mask_resp(iCell));
    end
    
    resptest_ind_all = [resptest_ind_all; find(h_resptest)+totCells];
    respmask_ind_all = [respmask_ind_all; find(h_respmask)+totCells];
    respplaid_ind_all = [respplaid_ind_all; find(h_respplaid)+totCells];
    respmask_ad_ind_all = [respmask_ad_ind_all; find(h_respmask_ad)+totCells];
    
    preftest_ind_all = [preftest_ind_all; intersect(find(h_resptest),find(h_preftest))+totCells];
    prefmask_ind_all = [prefmask_ind_all; intersect(find(h_respmask),find(h_prefmask))+totCells];
    
    h_MI_all = [h_MI_all; h_MI];
    p_MI_all = [p_MI_all; p_MI];
    h_MI_ad_all = [h_MI_ad_all; h_MI_ad];
    p_MI_ad_all = [p_MI_ad_all; p_MI_ad];
    
    prewin_frames = 15;
    postwin_frames = 45;
    tt = (-prewin_frames:postwin_frames-1).*(1/frame_rate);
    data_test = nan(prewin_frames+postwin_frames,nCells,nTrials);

    for itrial = 1:nTrials
        if cTest(itrial) + postwin_frames < sz(3)
            data_test(:,:,itrial) = npSub_tc(cTest(itrial)-prewin_frames:cTest(itrial)+postwin_frames-1,:);
        end
    end

    data_f = mean(data_test(1:prewin_frames,:,:),1);
    data_dfof_tc = (data_test-data_f)./data_f;
    noadapt_resp_tc = nan(prewin_frames+postwin_frames,nCells,nTest,nTest,2);
    singadapt_resp_tc = nan(prewin_frames+postwin_frames,nCells,nTest,nTest,2);
    ind_singadapt = setdiff(1:nTrials,ind_noadapt);
    for im = 1:nTest
        ind_mask = find(maskCon == maskCons(im));
        for it = 1:nTest
            ind_test = find(testCon == testCons(it));
            ind = intersect(ind_noadapt,intersect(ind_test,ind_mask));
            noadapt_resp_tc(:,:,im,it,1) = nanmean(data_dfof_tc(:,:,ind),3);
            noadapt_resp_tc(:,:,im,it,2) = nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind));
            ind = intersect(ind_singadapt,intersect(ind_test,ind_mask));
            singadapt_resp_tc(:,:,im,it,1) = nanmean(data_dfof_tc(:,:,ind),3);
            singadapt_resp_tc(:,:,im,it,2) = nanstd(data_dfof_tc(:,:,ind),[],3)./sqrt(length(ind));
        end
    end
    
    noadapt_resp_tc_all = cat(2, noadapt_resp_tc_all, noadapt_resp_tc);
    singadapt_resp_tc_all = cat(2, singadapt_resp_tc_all, singadapt_resp_tc);
    
    totCells = totCells+size(noadapt_resp_cell{1,1},1); 
end
%%

nmice = length([{expt.mouse}]);

noadapt_mask_resp = noadapt_resp_cell_all{5,1}-noadapt_resp_cell_all{1,1};
noadapt_mask_resp(find(noadapt_mask_resp<0)) = 0;
noadapt_test_resp = noadapt_resp_cell_all{1,5}-noadapt_resp_cell_all{1,1};
noadapt_test_resp(find(noadapt_test_resp<0)) = 0;
noadapt_plaid_resp =noadapt_resp_cell_all{5,5}-noadapt_resp_cell_all{1,1};
noadapt_plaid_resp(find(noadapt_plaid_resp<0)) = 0;

noadapt_SI_all = abs(noadapt_mask_resp-noadapt_test_resp)./(noadapt_mask_resp+noadapt_test_resp);
noadapt_MI_all = (noadapt_plaid_resp-(noadapt_mask_resp+noadapt_test_resp))./(noadapt_plaid_resp+(noadapt_mask_resp+noadapt_test_resp));

singadapt_mask_resp = singadapt_resp_cell_all{5,1}-singadapt_resp_cell_all{1,1};
singadapt_mask_resp(find(singadapt_mask_resp<0)) = 0;
singadapt_test_resp = singadapt_resp_cell_all{1,5}-singadapt_resp_cell_all{1,1};
singadapt_test_resp(find(singadapt_test_resp<0)) = 0;
singadapt_plaid_resp =singadapt_resp_cell_all{5,5}-singadapt_resp_cell_all{1,1};
singadapt_plaid_resp(find(singadapt_plaid_resp<0)) = 0;

singadapt_SI_all = abs(singadapt_mask_resp-singadapt_test_resp)./(singadapt_mask_resp+singadapt_test_resp);
singadapt_MI_all = (singadapt_plaid_resp-(singadapt_mask_resp+singadapt_test_resp))./(singadapt_plaid_resp+(singadapt_mask_resp+singadapt_test_resp));

figure;
subplot(2,2,1)
histogram(noadapt_MI_all(preftest_ind_all),[-1:0.1:1])
xlim([-1 1])
xlabel('Masking Index')
ylabel('Number of cells')
title(['Adapter preferring- n = ' num2str(length(preftest_ind_all)) ' cells'])
subplot(2,2,2)
histogram(noadapt_MI_all(prefmask_ind_all),[-1:0.1:1])
xlim([-1 1])
xlabel('Masking Index')
ylabel('Number of cells')
title(['Ortho preferring- n = ' num2str(length(prefmask_ind_all)) ' cells'])
subplot(2,2,3)
histogram(singadapt_MI_all(preftest_ind_all),[-1:0.1:1])
xlim([-1 1])
xlabel('Masking Index')
ylabel('Number of cells')
title(['Adapted'])
subplot(2,2,4)
histogram(singadapt_MI_all(prefmask_ind_all),[-1:0.1:1])
xlim([-1 1])
xlabel('Masking Index')
ylabel('Number of cells')
title(['Adapted'])
suptitle(['Single-stim adapt- n = ' num2str(nmice) ' mice'])
print(fullfile(summaryDir_F6, 'Figure6_MI_histograms.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_MI_histograms.fig'))

ind_h = intersect(find(noadapt_MI_all>0),prefmask_ind_all);
ind_l = intersect(find(noadapt_MI_all<0),prefmask_ind_all);
ind_h_sig = intersect(intersect(find(noadapt_MI_all>0),prefmask_ind_all),respmask_ad_ind_all);
ind_l_sig = intersect(intersect(find(noadapt_MI_all<0),prefmask_ind_all),respmask_ad_ind_all);
[temp ind_h_sigsub] = intersect(ind_h,ind_h_sig);
[temp ind_l_sigsub] = intersect(ind_l,ind_l_sig);


for ic = 41:45;
    iC  = ind_h(ic);
    figure;
iexp = expt_ind(iC);

movegui('center')
ii = 1;
for im = 1:nTest
    for it = 1:nTest
        subplot(nTest,nTest,ii)
        shadedErrorBar(tt, noadapt_resp_tc_all(:,iC,im,it,1), noadapt_resp_tc_all(:,iC,im,it,2));
        hold on
        shadedErrorBar(tt, singadapt_resp_tc_all(:,iC,im,it,1), singadapt_resp_tc_all(:,iC,im,it,2), {'r-','markerfacecolor','r'});
        ii = ii+1;
        ylim([-0.5 2])
        title(['T=' num2str(testCons(it)) '; M=' num2str(testCons(im))]) 
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date '- Cell #' num2str(iC)])
end
print(fullfile(summaryDir_F6, 'Figure6_ExampleCell_posMod.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ExampleCell_posMod.fig'))


iC = ind_l(27); %13 is also nice
iexp = expt_ind(iC);
figure;
movegui('center')
ii = 1;
for im = 1:nTest
    for it = 1:nTest
        subplot(nTest,nTest,ii)
        shadedErrorBar(tt, noadapt_resp_tc_all(:,iC,im,it,1), noadapt_resp_tc_all(:,iC,im,it,2));
        hold on
        shadedErrorBar(tt, singadapt_resp_tc_all(:,iC,im,it,1), singadapt_resp_tc_all(:,iC,im,it,2), {'r-','markerfacecolor','r'});
        ii = ii+1;
        ylim([-0.5 2])
        title(['T=' num2str(testCons(it)) '; M=' num2str(testCons(im))]) 
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date '- Cell #' num2str(iC)])
print(fullfile(summaryDir_F6, 'Figure6_ExampleCell_negMod.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ExampleCell_negMod.fig'))

noadapt_preftest_resp_avg = zeros(nTest,nTest,2);
noadapt_prefmask_resp_avg = zeros(nTest,nTest,2);
noadapt_nopref_resp_avg = zeros(nTest,nTest,2);
singadapt_preftest_resp_avg = zeros(nTest,nTest,2);
singadapt_prefmask_resp_avg = zeros(nTest,nTest,2);
singadapt_nopref_resp_avg = zeros(nTest,nTest,2);

noadapt_prefmask_posMod_resp_avg = zeros(nTest,nTest,2);
noadapt_prefmask_negMod_resp_avg = zeros(nTest,nTest,2);
singadapt_prefmask_posMod_resp_avg = zeros(nTest,nTest,2);
singadapt_prefmask_negMod_resp_avg = zeros(nTest,nTest,2);
noadapt_prefmask_sigPosMod_resp_avg = zeros(nTest,nTest,2);
noadapt_prefmask_sigNegMod_resp_avg = zeros(nTest,nTest,2);
singadapt_prefmask_sigPosMod_resp_avg = zeros(nTest,nTest,2);
singadapt_prefmask_sigNegMod_resp_avg = zeros(nTest,nTest,2);

nopref_ind_all = setdiff(unique([resptest_ind_all; respmask_ind_all]),unique([preftest_ind_all; prefmask_ind_all]));
preftest_resp_all = [];
prefmask_resp_all = [];
posMod_resp_all = [];
negMod_resp_all = [];
stimCond_test = [];
stimCond_mask = [];
stimCond_pos = [];
stimCond_neg = [];
% preftest_resp_nomask_all = [];
% prefmask_resp_nomask_all = [];
% preftest_resp_mask_all = [];
% prefmask_resp_mask_all = [];
% stimCond_nomask_test = [];
% stimCond_nomask_mask = [];
% stimCond_mask_test = [];
% stimCond_mask_mask = [];

for im = 1:nTest
    for it = 1:nTest
        noadapt_preftest_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(preftest_ind_all,:),1);
        noadapt_preftest_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(preftest_ind_all,:),[],1)./sqrt(size(preftest_ind_all,1));
        noadapt_prefmask_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(prefmask_ind_all,:),1);
        noadapt_prefmask_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(prefmask_ind_all,:),[],1)./sqrt(size(prefmask_ind_all,1));
        noadapt_nopref_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(nopref_ind_all,:),1);
        noadapt_nopref_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(nopref_ind_all,:),[],1)./sqrt(size(nopref_ind_all,1));
        singadapt_preftest_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(preftest_ind_all,:),1);
        singadapt_preftest_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(preftest_ind_all,:),[],1)./sqrt(size(preftest_ind_all,1));
        singadapt_prefmask_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(prefmask_ind_all,:),1);
        singadapt_prefmask_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(prefmask_ind_all,:),[],1)./sqrt(size(prefmask_ind_all,1));
        singadapt_nopref_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(nopref_ind_all,:),1);
        singadapt_nopref_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(nopref_ind_all,:),[],1)./sqrt(size(nopref_ind_all,1));
        noadapt_prefmask_posMod_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(ind_h,:),1);
        noadapt_prefmask_posMod_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(ind_h,:),[],1)./sqrt(length(ind_h));
        noadapt_prefmask_negMod_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(ind_l,:),1);
        noadapt_prefmask_negMod_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(ind_l,:),[],1)./sqrt(length(ind_l));
        singadapt_prefmask_posMod_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(ind_h,:),1);
        singadapt_prefmask_posMod_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(ind_h,:),[],1)./sqrt(length(ind_h));
        singadapt_prefmask_negMod_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(ind_l,:),1);
        singadapt_prefmask_negMod_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(ind_l,:),[],1)./sqrt(length(ind_l));
        noadapt_prefmask_sigPosMod_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(ind_h_sig,:),1);
        noadapt_prefmask_sigPosMod_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(ind_h_sig,:),[],1)./sqrt(length(ind_h_sig));
        noadapt_prefmask_sigNegMod_resp_avg(im,it,1) = mean(noadapt_resp_cell_all{im,it}(ind_l_sig,:),1);
        noadapt_prefmask_sigNegMod_resp_avg(im,it,2) = std(noadapt_resp_cell_all{im,it}(ind_l_sig,:),[],1)./sqrt(length(ind_l));
        singadapt_prefmask_sigPosMod_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(ind_h_sig,:),1);
        singadapt_prefmask_sigPosMod_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(ind_h_sig,:),[],1)./sqrt(length(ind_h_sig));
        singadapt_prefmask_sigNegMod_resp_avg(im,it,1) = mean(singadapt_resp_cell_all{im,it}(ind_l_sig,:),1);
        singadapt_prefmask_sigNegMod_resp_avg(im,it,2) = std(singadapt_resp_cell_all{im,it}(ind_l_sig,:),[],1)./sqrt(length(ind_l_sig));
        
        n_test = length(preftest_ind_all);
        n_mask = length(prefmask_ind_all);
        preftest_resp_all = [preftest_resp_all noadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
        prefmask_resp_all = [prefmask_resp_all noadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
        preftest_resp_all = [preftest_resp_all singadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
        prefmask_resp_all = [prefmask_resp_all singadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
        
%         if it == 1
%             preftest_resp_nomask_all = [preftest_resp_nomask_all noadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
%             prefmask_resp_nomask_all = [prefmask_resp_nomask_all noadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
%             preftest_resp_nomask_all = [preftest_resp_nomask_all singadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
%             prefmask_resp_nomask_all = [prefmask_resp_nomask_all singadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
%             stimCond_nomask_test = [stimCond_nomask_test; it.*ones(n_test,1) im.*ones(n_test,1) zeros(n_test,1)];
%             stimCond_nomask_mask = [stimCond_nomask_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) zeros(n_mask,1)];
%             stimCond_nomask_test = [stimCond_nomask_test; it.*ones(n_test,1) im.*ones(n_test,1) ones(n_test,1)];
%             stimCond_nomask_mask = [stimCond_nomask_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) ones(n_mask,1)];
%         end
%         
%         if it == 5
%             preftest_resp_mask_all = [preftest_resp_mask_all noadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
%             prefmask_resp_mask_all = [prefmask_resp_mask_all noadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
%             preftest_resp_mask_all = [preftest_resp_mask_all singadapt_resp_cell_all{im,it}(preftest_ind_all,:)'];
%             prefmask_resp_mask_all = [prefmask_resp_mask_all singadapt_resp_cell_all{im,it}(prefmask_ind_all,:)'];
%             stimCond_mask_test = [stimCond_mask_test; it.*ones(n_test,1) im.*ones(n_test,1) zeros(n_test,1)];
%             stimCond_mask_mask = [stimCond_mask_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) zeros(n_mask,1)];
%             stimCond_mask_test = [stimCond_mask_test; it.*ones(n_test,1) im.*ones(n_test,1) ones(n_test,1)];
%             stimCond_mask_mask = [stimCond_mask_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) ones(n_mask,1)];
%         end
        
        stimCond_test = [stimCond_test; it.*ones(n_test,1) im.*ones(n_test,1) zeros(n_test,1)];
        stimCond_mask = [stimCond_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) zeros(n_mask,1)];
        stimCond_test = [stimCond_test; it.*ones(n_test,1) im.*ones(n_test,1) ones(n_test,1)];
        stimCond_mask = [stimCond_mask; it.*ones(n_mask,1) im.*ones(n_mask,1) ones(n_mask,1)];
        
        n_pos = length(ind_h);
        n_neg = length(ind_l);
        posMod_resp_all = [posMod_resp_all noadapt_resp_cell_all{im,it}(ind_h,:)'];
        negMod_resp_all = [negMod_resp_all noadapt_resp_cell_all{im,it}(ind_l,:)'];
        posMod_resp_all = [posMod_resp_all singadapt_resp_cell_all{im,it}(ind_h,:)'];
        negMod_resp_all = [negMod_resp_all singadapt_resp_cell_all{im,it}(ind_l,:)'];
        stimCond_pos = [stimCond_pos; it.*ones(n_pos,1) im.*ones(n_pos,1) zeros(n_pos,1)];
        stimCond_neg = [stimCond_neg; it.*ones(n_neg,1) im.*ones(n_neg,1) zeros(n_neg,1)];
        stimCond_pos = [stimCond_pos; it.*ones(n_pos,1) im.*ones(n_pos,1) ones(n_pos,1)];
        stimCond_neg = [stimCond_neg; it.*ones(n_neg,1) im.*ones(n_neg,1) ones(n_neg,1)];
    end
end

%compare test preferring in control and adapt with no mask
ind_test = find(stimCond_test(:,2)==1);
test_mat = stimCond_test(ind_test,1);
adapt_mat = stimCond_test(ind_test,3);
[h_test_nomask, p_test_nomask, stats_test] = anovan(preftest_resp_all(ind_test)',{test_mat,adapt_mat});

%compare mask preferring in control and adapt with no mask
ind_mask = find(stimCond_mask(:,1)==1);
mask_mat = stimCond_mask(ind_mask,2);
adapt_mat = stimCond_mask(ind_mask,3);
[h_mask_nomask, p_mask_nomask, stats_mask] = anovan(prefmask_resp_all(ind_mask)',{mask_mat,adapt_mat});

%compare test preferring in control and adapt with mask
ind_test = find(stimCond_test(:,2)==5);
test_mat = stimCond_test(ind_test,1);
adapt_mat = stimCond_test(ind_test,3);
[h_test_mask, p_test_mask, stats_test] = anovan(preftest_resp_all(ind_test)',{test_mat,adapt_mat});

%compare mask preferring in control and adapt with mask
ind_mask = find(stimCond_mask(:,1)==5);
mask_mat = stimCond_mask(ind_mask,2);
adapt_mat = stimCond_mask(ind_mask,3);
[h_mask_mask, p_mask_mask, stats_mask] = anovan(prefmask_resp_all(ind_mask)',{mask_mat,adapt_mat});

ind_pos = find(stimCond_pos(:,1)==5);
pos_mat = stimCond_pos(ind_pos,2);
adapt_mat = stimCond_pos(ind_pos,3);
[h_pos, p_pos, stats_pos] = anovan(posMod_resp_all(ind_pos)',{pos_mat,adapt_mat});

ind_neg = find(stimCond_neg(:,1)==5);
neg_mat = stimCond_neg(ind_neg,2);
adapt_mat = stimCond_neg(ind_neg,3);
[h_neg, p_neg, stats_neg] = anovan(negMod_resp_all(ind_neg)',{neg_mat,adapt_mat});

%main effect of mask
mask_mat = stimCond_pos(:,1);
pos_mat = stimCond_pos(:,2);
adapt_mat = stimCond_pos(:,3);
[h_pos, p_pos, stats_pos] = anovan(posMod_resp_all',{mask_mat,pos_mat,adapt_mat});

n = find(isnan(noadapt_MI_all+singadapt_MI_all));

[h_MI_pos, p_MI_pos_paired] = ttest(noadapt_MI_all(ind_h),singadapt_MI_all(ind_h));
[h_MI_neg, p_MI_neg_paired] = ttest(noadapt_MI_all(ind_l),singadapt_MI_all(ind_l));

n_pos = length(find(sum(~isnan([noadapt_MI_all(ind_h) singadapt_MI_all(ind_h)]),2)==2));
n_neg = length(find(sum(~isnan([noadapt_MI_all(ind_l) singadapt_MI_all(ind_l)]),2)==2));

MI_no_pos_avg = mean(noadapt_MI_all(ind_h));
MI_no_pos_sem = std(noadapt_MI_all(ind_h))./sqrt(length(ind_h));
MI_no_neg_avg = mean(noadapt_MI_all(ind_l));
MI_no_neg_sem = std(noadapt_MI_all(ind_l))./sqrt(length(ind_l));
MI_sing_pos_avg = nanmean(singadapt_MI_all(ind_h));
MI_sing_pos_sem = nanstd(singadapt_MI_all(ind_h))./sqrt(length(ind_h));
MI_sing_neg_avg = nanmean(singadapt_MI_all(ind_l));
MI_sing_neg_sem = nanstd(singadapt_MI_all(ind_l))./sqrt(length(ind_l));

[h_MI_pos_sig, p_MI_pos_sig_paired] = ttest(noadapt_MI_all(ind_h_sig),singadapt_MI_all(ind_h_sig));
[h_MI_neg_sig, p_MI_neg_sig_paired] = ttest(noadapt_MI_all(ind_l_sig),singadapt_MI_all(ind_l_sig));
n_pos_sig = length(find(sum(~isnan([noadapt_MI_all(ind_h_sig) singadapt_MI_all(ind_h_sig)]),2)==2));
n_neg_sig = length(find(sum(~isnan([noadapt_MI_all(ind_l_sig) singadapt_MI_all(ind_l_sig)]),2)==2));

MI_no_pos_avg_sig = mean(noadapt_MI_all(ind_h_sig));
MI_no_pos_sem_sig = std(noadapt_MI_all(ind_h_sig))./sqrt(length(ind_h_sig));
MI_no_neg_avg_sig = mean(noadapt_MI_all(ind_l_sig));
MI_no_neg_sem_sig = std(noadapt_MI_all(ind_l_sig))./sqrt(length(ind_l_sig));
MI_sing_pos_avg_sig = nanmean(singadapt_MI_all(ind_h_sig));
MI_sing_pos_sem_sig = nanstd(singadapt_MI_all(ind_h_sig))./sqrt(length(ind_h_sig));
MI_sing_neg_avg_sig = nanmean(singadapt_MI_all(ind_l_sig));
MI_sing_neg_sem_sig = nanstd(singadapt_MI_all(ind_l_sig))./sqrt(length(ind_l_sig));

r_MI_pos = triu2vec(corrcoef(noadapt_MI_all(setdiff(ind_h,n)),singadapt_MI_all(setdiff(ind_h,n))));
r_MI_pos_sig = triu2vec(corrcoef(noadapt_MI_all(intersect(find(h_MI_all),setdiff(ind_h,n))),singadapt_MI_all(intersect(find(h_MI_all),setdiff(ind_h,n)))));

r_MI_neg = triu2vec(corrcoef(noadapt_MI_all(setdiff(ind_l,n)),singadapt_MI_all(setdiff(ind_l,n))));
r_MI_neg_sig = triu2vec(corrcoef(noadapt_MI_all(setdiff(ind_l_sig,n)),singadapt_MI_all(setdiff(ind_l_sig,n))));

r_MI = triu2vec(corrcoef(noadapt_MI_all(setdiff([ind_h; ind_l],n)),singadapt_MI_all(setdiff([ind_h; ind_l],n))));
r_MI_sig = triu2vec(corrcoef(noadapt_MI_all(setdiff([ind_h_sig; ind_l_sig],n)),singadapt_MI_all(setdiff([ind_h_sig; ind_l_sig],n))));

[h_MI_pos, p_MI_pos_students] = ttest(singadapt_MI_all(ind_h),0,'tail','right');
[h_MI_neg, p_MI_neg_students] = ttest(singadapt_MI_all(ind_l),0,'tail','left');

[h_MI_pos_sig, p_MI_pos_sig_students] = ttest(singadapt_MI_all(ind_h_sig),0,'tail','right');
[h_MI_neg_sig, p_MI_neg_sig_students] = ttest(singadapt_MI_all(ind_l_sig),0,'tail','left');


figure;
subplot(1,3,1)
errorbar(testCons,noadapt_preftest_resp_avg(1,:,1),noadapt_preftest_resp_avg(1,:,2),'-ok')
hold on
errorbar(testCons,singadapt_preftest_resp_avg(1,:,1),singadapt_preftest_resp_avg(1,:,2),'-or')
errorbar(testCons,noadapt_preftest_resp_avg(:,1,1),noadapt_preftest_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_preftest_resp_avg(:,1,1),singadapt_preftest_resp_avg(:,1,2),'-om')
ylim([-0.1 0.5])
xlabel('Contrast')
ylabel('dF/F')
title(['Adapter preferring- n = ' num2str(size(preftest_ind_all,1))])
legend({'No mask- adapter','Adapt- adapter','No mask- ortho','Adapt- ortho'},'location','northwest')
subplot(1,3,2)
errorbar(testCons,noadapt_nopref_resp_avg(1,:,1),noadapt_nopref_resp_avg(1,:,2),'-ok')
hold on
errorbar(testCons,singadapt_nopref_resp_avg(1,:,1),singadapt_nopref_resp_avg(1,:,2),'-or')
errorbar(testCons,noadapt_nopref_resp_avg(:,1,1),noadapt_nopref_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_nopref_resp_avg(:,1,1),singadapt_nopref_resp_avg(:,1,2),'-om')
ylim([-0.1 0.5])
xlabel('Contrast')
ylabel('dF/F')
title(['No preference- n = ' num2str(size(nopref_ind_all,1))])
subplot(1,3,3)
errorbar(testCons,noadapt_prefmask_resp_avg(1,:,1),noadapt_prefmask_resp_avg(1,:,2),'-ok')
hold on
errorbar(testCons,singadapt_prefmask_resp_avg(1,:,1),singadapt_prefmask_resp_avg(1,:,2),'-or')
errorbar(testCons,noadapt_prefmask_resp_avg(:,1,1),noadapt_prefmask_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_prefmask_resp_avg(:,1,1),singadapt_prefmask_resp_avg(:,1,2),'-om')
ylim([-0.1 0.5])
xlabel('Contrast')
ylabel('dF/F')
title(['Ortho preferring- n = ' num2str(size(prefmask_ind_all,1))])
suptitle(['Single-stim adapt- n = ' num2str(nmice) ' mice'])
print(fullfile(summaryDir_F6, 'Figure6_avgResp_AllAdaptStimCombo.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_avgResp_AllAdaptStimCombo.fig'))

figure;
subplot(2,2,1)
errorbar(testCons,noadapt_preftest_resp_avg(1,:,1),noadapt_preftest_resp_avg(1,:,2),'-ok')
hold on
errorbar(testCons,singadapt_preftest_resp_avg(1,:,1),singadapt_preftest_resp_avg(1,:,2),'-or')
ylim([-0.1 0.5])
xlabel('Adapter Contrast')
ylabel('dF/F')
title(['Adapter preferring- n = ' num2str(size(preftest_ind_all,1))])
legend({'No mask','Adapt'},'location','northwest')
subplot(2,2,2)
errorbar(testCons,noadapt_prefmask_resp_avg(:,1,1),noadapt_prefmask_resp_avg(:,1,2),'-ok')
hold on
errorbar(testCons,singadapt_prefmask_resp_avg(:,1,1),singadapt_prefmask_resp_avg(:,1,2),'-or')
ylim([-0.1 0.5])
xlabel('Ortho Contrast')
ylabel('dF/F')
legend({'No mask','Adapt'},'location','northwest')
title(['Ortho preferring- n = ' num2str(size(prefmask_ind_all,1))])
subplot(2,2,3)
errorbar(testCons,noadapt_preftest_resp_avg(end,:,1),noadapt_preftest_resp_avg(end,:,2),'-ok')
hold on
errorbar(testCons,singadapt_preftest_resp_avg(end,:,1),singadapt_preftest_resp_avg(end,:,2),'-or')
ylim([-0.1 0.5])
xlabel('Adapter Contrast')
ylabel('dF/F')
legend({['Mask (ortho) = ' num2str(testCons(end))],'Adapt'},'location','northwest')
subplot(2,2,4)
errorbar(testCons,noadapt_prefmask_resp_avg(:,end,1),noadapt_prefmask_resp_avg(:,end,2),'-ok')
hold on
errorbar(testCons,singadapt_prefmask_resp_avg(:,end,1),singadapt_prefmask_resp_avg(:,end,2),'-or')
ylim([-0.1 0.5])
xlabel('Ortho Contrast')
ylabel('dF/F')
legend({['Mask (adapter) = ' num2str(testCons(end))],'Adapt'},'location','northwest')
suptitle(['Single-stim adapt- n = ' num2str(nmice) ' mice'])
print(fullfile(summaryDir_F6, 'Figure6_avgResp_AdaptVOrthogPref.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_avgResp_AdaptVOrthogPref.fig'))

figure;
subplot(2,2,1)
errorbar(testCons,noadapt_prefmask_posMod_resp_avg(:,end,1),noadapt_prefmask_posMod_resp_avg(:,end,2),'-ok')
hold on
errorbar(testCons,noadapt_prefmask_posMod_resp_avg(:,1,1),noadapt_prefmask_posMod_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_prefmask_posMod_resp_avg(:,end,1),singadapt_prefmask_posMod_resp_avg(:,end,2),'-or')
ylim([-0.1 1])
xlabel('Contrast')
ylabel('dF/F')
legend({['Mask = ' num2str(testCons(end))],['Mask = ' num2str(testCons(1))],'Adapt'},'location','northwest')
title(['Pos mod- n = ' num2str(length(ind_h))])
subplot(2,2,2)
errorbar(testCons,noadapt_prefmask_negMod_resp_avg(:,end,1),noadapt_prefmask_negMod_resp_avg(:,end,2),'-ok')
hold on
errorbar(testCons,noadapt_prefmask_negMod_resp_avg(:,1,1),noadapt_prefmask_negMod_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_prefmask_negMod_resp_avg(:,end,1),singadapt_prefmask_negMod_resp_avg(:,end,2),'-or')
ylim([-0.1 0.5])
xlabel('Contrast')
ylabel('dF/F')
legend({['Mask = ' num2str(testCons(end))],['Mask = ' num2str(testCons(1))],'Adapt'},'location','northwest')
title(['Neg mod- n = ' num2str(length(ind_l))])

MI_all_diff = -noadapt_MI_all+singadapt_MI_all;
MI_all_diff_zero = -noadapt_MI_all;

subplot(2,2,3)
cdfplot(MI_all_diff(ind_h))
hold on
cdfplot(MI_all_diff_zero(ind_h))
legend({'Actual','If adapt->MI=0'},'location','southeast')
xlabel('Adapt-Control MI')
ylabel('Fraction of cells')
xlim([-2 2])
title('')
subplot(2,2,4)
cdfplot(MI_all_diff(ind_l))
hold on
cdfplot(MI_all_diff_zero(ind_l))
legend({'Actual','If adapt->MI=0'},'location','northwest')
xlabel('Adapt-Control MI')
ylabel('Fraction of cells')
xlim([-2 2])
title('')
suptitle(['Single-stim adapt- n = ' num2str(nmice) ' mice'])

print(fullfile(summaryDir_F6, 'Figure6_avgResp_OrthogPref_PosVNegMod.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_avgResp_OrthogPref_PosVNegMod.fig'))

figure;
subplot(2,2,1)
errorbar(testCons,noadapt_prefmask_sigPosMod_resp_avg(:,end,1),noadapt_prefmask_sigPosMod_resp_avg(:,end,2),'-ok')
hold on
errorbar(testCons,noadapt_prefmask_sigPosMod_resp_avg(:,1,1),noadapt_prefmask_sigPosMod_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_prefmask_sigPosMod_resp_avg(:,end,1),singadapt_prefmask_sigPosMod_resp_avg(:,end,2),'-or')
ylim([-0.1 1])
xlabel('Contrast')
ylabel('dF/F')
legend({['Mask = ' num2str(testCons(end))],['Mask = ' num2str(testCons(1))],'Adapt'},'location','northwest')
title(['Sig pos mod- n = ' num2str(length(ind_h_sig))])
subplot(2,2,2)
errorbar(testCons,noadapt_prefmask_sigNegMod_resp_avg(:,end,1),noadapt_prefmask_sigNegMod_resp_avg(:,end,2),'-ok')
hold on
errorbar(testCons,noadapt_prefmask_sigNegMod_resp_avg(:,1,1),noadapt_prefmask_sigNegMod_resp_avg(:,1,2),'-ob')
errorbar(testCons,singadapt_prefmask_sigNegMod_resp_avg(:,end,1),singadapt_prefmask_sigNegMod_resp_avg(:,end,2),'-or')
ylim([-0.1 0.5])
xlabel('Contrast')
ylabel('dF/F')
legend({['Mask = ' num2str(testCons(end))],['Mask = ' num2str(testCons(1))],'Adapt'},'location','northwest')
title(['Sig neg mod- n = ' num2str(length(ind_l_sig))])



figure;
subplot(2,2,1)
scatter(noadapt_MI_all(ind_l), singadapt_MI_all(ind_l))
hold on
scatter(noadapt_MI_all(ind_h), singadapt_MI_all(ind_h))
vline(0)
hline(0)
xlim([-1 1])
ylim([-1 1])
xlabel('MI- Control')
ylabel('MI- Adapt')

[n_nol edges bin_nol] = histcounts(noadapt_MI_all(ind_l),[-1:0.2:1]);%,'Normalization','pdf');
[n_noh edges bin_noh] = histcounts(noadapt_MI_all(ind_h),[-1:0.2:1]);%,'Normalization','pdf');
subplot(2,3,3)
bar(edges(1:end-1),n_nol./length(bin_nol))
hold on
bar(edges(1:end-1),n_noh./length(bin_noh))
xlabel('MI- Control')
ylabel('Fraction N')
vline(nanmean(noadapt_MI_all(ind_l)),'b')
vline(nanmean(noadapt_MI_all(ind_h)),'r')
ylim([0 1])

[n_adl edges bin_adl] = histcounts(singadapt_MI_all(ind_l),[-1:0.2:1]);%,'Normalization','pdf');
[n_adh edges bin_adh] = histcounts(singadapt_MI_all(ind_h),[-1:0.2:1]);%,'Normalization','pdf');
subplot(2,2,4)
bar(edges(1:end-1),n_adl./length(bin_adl))
hold on
bar(edges(1:end-1),n_adh./length(bin_adh))
xlabel('MI- Adapt')
ylabel('Fraction N')
ylim([0 1])
vline(nanmean(singadapt_MI_all(ind_l)),'b')
vline(nanmean(singadapt_MI_all(ind_h)),'r')

suptitle(['Orthogonal responsive- n = ' num2str(length([ind_l; ind_h]))])

print(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter&PDFs.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter&PDFs.fig'))

figure;
subplot(2,2,1)
scatter(singadapt_SI_all(ind_l)-noadapt_SI_all(ind_l),singadapt_MI_all(ind_l)-noadapt_MI_all(ind_l))
xlim([-1 1])
ylim([-1 1])
xlabel('SI- Adapt-Control')
ylabel('MI- Adapt-Control')
title('Suppressed')
subplot(2,2,2)
scatter(singadapt_SI_all(ind_h)-noadapt_SI_all(ind_h),singadapt_MI_all(ind_h)-noadapt_MI_all(ind_h))
xlim([-1 1])
ylim([-1 1])
title('Facilitated')
xlabel('SI- Adapt-Control')
ylabel('MI- Adapt-Control')
subplot(2,2,3)
ind0 = find(singadapt_SI_all(ind_l)-noadapt_SI_all(ind_l));
scatter(singadapt_SI_all(ind_l(ind0))-noadapt_SI_all(ind_l(ind0)),singadapt_MI_all(ind_l(ind0))-noadapt_MI_all(ind_l(ind0)))
xlim([-1 1])
ylim([-1 1])
xlabel('SI- Adapt-Control')
ylabel('MI- Adapt-Control')
title('Suppressed- SI=0 removed')
print(fullfile(summaryDir_F6, 'Figure6_diffSIMI_scatters.pdf'),'-dpdf','-bestfit')

%sig only
figure;
subplot(2,2,1)
scatter(noadapt_MI_all(ind_l_sig), singadapt_MI_all(ind_l_sig))
hold on
scatter(noadapt_MI_all(ind_h_sig), singadapt_MI_all(ind_h_sig))
vline(0)
hline(0)
xlim([-1 1])
ylim([-1 1])
xlabel('MI- Control')
ylabel('MI- Adapt')

[n_nol edges bin_nol] = histcounts(noadapt_MI_all(ind_l_sig),[-1:0.2:1]);%,'Normalization','pdf');
[n_noh edges bin_noh] = histcounts(noadapt_MI_all(ind_h_sig),[-1:0.2:1]);%,'Normalization','pdf');
subplot(2,2,3)
bar(edges(1:end-1),n_nol./length(bin_nol))
hold on
bar(edges(1:end-1),n_noh./length(bin_noh))
xlabel('MI- Control')
ylabel('Fraction N')
vline(nanmean(noadapt_MI_all(ind_l_sig)),'b')
vline(nanmean(noadapt_MI_all(ind_h_sig)),'r')
ylim([0 1])

[n_adl edges bin_adl] = histcounts(singadapt_MI_all(ind_l_sig),[-1:0.2:1]);%,'Normalization','pdf');
[n_adh edges bin_adh] = histcounts(singadapt_MI_all(ind_h_sig),[-1:0.2:1]);%,'Normalization','pdf');
subplot(2,2,4)
bar(edges(1:end-1),n_adl./length(bin_adl))
hold on
bar(edges(1:end-1),n_adh./length(bin_adh))
xlabel('MI- Adapt')
ylabel('Fraction N')
ylim([0 1])
vline(nanmean(singadapt_MI_all(ind_l_sig)),'b')
vline(nanmean(singadapt_MI_all(ind_h_sig)),'r')

suptitle(['Orthogonal responsive- n = ' num2str(length([ind_l_sig; ind_h_sig]))])

print(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter&PDFs_sig.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter&PDFs_sig.fig'))


figure;
subplot(2,2,1)
scatter(noadapt_MI_all(ind_l), singadapt_MI_all(ind_l))
hold on
scatter(noadapt_MI_all(ind_h), singadapt_MI_all(ind_h))
vline(0)
hline(0)
xlim([-1 1])
ylim([-1 1])
xlabel('MI- Control')
ylabel('MI- Adapt')
title(['All- n = ' num2str(length([ind_l; ind_h]))])

subplot(2,2,2)
scatter(noadapt_MI_all(ind_l_sig), singadapt_MI_all(ind_l_sig))
hold on
scatter(noadapt_MI_all(ind_h_sig), singadapt_MI_all(ind_h_sig))
vline(0)
hline(0)
xlim([-1 1])
ylim([-1 1])
xlabel('MI- Control')
ylabel('MI- Adapt')
title(['Resp after adapt- n = ' num2str(length([ind_l_sig; ind_h_sig]))])

print(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter_sig.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ControlAdaptScatter_sig.fig'))

[n_nol edges bin_nol] = histcounts(noadapt_MI_all(ind_l),[-1:0.2:1]);%,'Normalization','pdf');
[n_noh edges bin_noh] = histcounts(noadapt_MI_all(ind_h),[-1:0.2:1]);%,'Normalization','pdf');
[n_adl edges bin_adl] = histcounts(singadapt_MI_all(ind_l),[-1:0.2:1]);%,'Normalization','pdf');
[n_adh edges bin_adh] = histcounts(singadapt_MI_all(ind_h),[-1:0.2:1]);%,'Normalization','pdf');
figure
subplot(2,2,1)
bar(-0.9:0.2:0.9,n_nol./length(bin_nol))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Suppressed- pre adapt')
vline(nanmean(noadapt_MI_all(ind_l)))
subplot(2,2,3)
bar(-0.9:0.2:0.9,n_noh./length(bin_noh))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Facilitated- pre adapt')
vline(nanmean(noadapt_MI_all(ind_h)))

subplot(2,2,2)
bar(-0.9:0.2:0.9,n_adl./length(bin_adl))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Suppressed- post adapt')
vline(nanmean(singadapt_MI_all(ind_l)))
subplot(2,2,4)
bar(-0.9:0.2:0.9,n_adh./length(bin_adh))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Facilitated- post adapt')
vline(nanmean(singadapt_MI_all(ind_h)))
suptitle('All cells')
print(fullfile(summaryDir_F6, 'Figure6_ControlAdaptHistograms.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ControlAdaptHistograms.fig'))

[n_sig_nol edges bin_sig_nol] = histcounts(noadapt_MI_all(ind_l_sig),[-1:0.2:1]);%,'Normalization','pdf');
[n_sig_noh edges bin_sig_noh] = histcounts(noadapt_MI_all(ind_h_sig),[-1:0.2:1]);%,'Normalization','pdf');
[n_sig_adl edges bin_sig_adl] = histcounts(singadapt_MI_all(ind_l_sig),[-1:0.2:1]);%,'Normalization','pdf');
[n_sig_adh edges bin_sig_adh] = histcounts(singadapt_MI_all(ind_h_sig),[-1:0.2:1]);%,'Normalization','pdf');
figure
subplot(2,2,1)
bar(-0.9:0.2:0.9,n_sig_nol./length(bin_sig_nol))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Suppressed- pre adapt')
vline(nanmean(noadapt_MI_all(ind_l_sig)))
subplot(2,2,3)
bar(-0.9:0.2:0.9,n_sig_noh./length(bin_sig_noh))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Facilitated- pre adapt')
vline(nanmean(noadapt_MI_all(ind_h_sig)))

subplot(2,2,2)
bar(-0.9:0.2:0.9,n_sig_adl./length(bin_sig_adl))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Suppressed- post adapt')
vline(nanmean(singadapt_MI_all(ind_l_sig)))
subplot(2,2,4)
bar(-0.9:0.2:0.9,n_sig_adh./length(bin_sig_adh))
ylim([0 1])
xlim([-1 1])
xlabel('MI- Control')
ylabel('Fraction N')
title('Facilitated- post adapt')
vline(nanmean(singadapt_MI_all(ind_h_sig)))
suptitle('Resp after adapt')
print(fullfile(summaryDir_F6, 'Figure6_ControlAdaptHistograms_sig.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F6, 'Figure6_ControlAdaptHistograms_sig.fig'))
