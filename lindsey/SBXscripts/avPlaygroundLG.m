nexp = size(dc_attn,2);
attn_rsc_weight_corr = zeros(2,2,nexp);
attn_rsc_all = [];
attn_vwT_all = [];
attn_awT_all = [];
attn_vwD_all = [];
attn_awD_all = [];
for iexp = 1:nexp
    VIx = strcmp(dc_attn(iexp).av(1).trOut,'cr');
    VIx = VIx+strcmp(dc_attn(iexp).av(1).trOut,'fa');
    AIx = strcmp(dc_attn(iexp).av(2).trOut,'cr');
    AIx = AIx+strcmp(dc_attn(iexp).av(2).trOut,'fa');
    resp_temp_v = dc_attn(iexp).av(1).respAllCells(find(VIx),:);
    resp_temp_a = dc_attn(iexp).av(2).respAllCells(find(AIx),:);
    resp_temp_vz = (resp_temp_v-mean(resp_temp_v,1))./std(resp_temp_v,[],1);
    resp_temp_az = (resp_temp_a-mean(resp_temp_a,1))./std(resp_temp_a,[],1);
    resp_corr_v = corrcoef(resp_temp_vz);
    resp_corr_a = corrcoef(resp_temp_az);
    resp_corr_avg = zeros(2,size(resp_temp_v,2));
    for ic = 1:size(resp_temp,2)
        temp = resp_corr_v;
        temp(:,ic) = [];
        resp_corr_avg(1,ic) = mean(temp(ic,:),2);
        temp = resp_corr_a;
        temp(:,ic) = [];
        resp_corr_avg(2,ic) = mean(temp(ic,:),2);
    end
    attn_rsc_weight_corr(1,1,iexp) = triu2vec(corrcoef(resp_corr_avg(1,:),dc_attn(iexp).av(1).weightTarget)');
    attn_rsc_weight_corr(2,1,iexp) = triu2vec(corrcoef(resp_corr_avg(2,:),dc_attn(iexp).av(2).weightTarget)');
    attn_rsc_weight_corr(1,2,iexp) = triu2vec(corrcoef(resp_corr_avg(1,:),dc_attn(iexp).av(1).weightDetect)');
    attn_rsc_weight_corr(2,2,iexp) = triu2vec(corrcoef(resp_corr_avg(2,:),dc_attn(iexp).av(2).weightDetect)');
    attn_rsc_all = [attn_rsc_all resp_corr_avg];
    attn_vwT_all = [attn_vwT_all dc_attn(iexp).av(1).weightTarget'];
    attn_awT_all = [attn_awT_all dc_attn(iexp).av(2).weightTarget'];
    attn_vwD_all = [attn_vwD_all dc_attn(iexp).av(1).weightDetect'];
    attn_awD_all = [attn_awD_all dc_attn(iexp).av(2).weightDetect'];
end

figure;
subplot(2,2,1)
scatter(attn_rsc_all(1,:),attn_vwT_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Vis target- ' num2str(chop(triu2vec(corrcoef(attn_rsc_all(1,:),attn_vwT_all)),2))])
subplot(2,2,2)
scatter(attn_rsc_all(2,:),attn_awT_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Aud target- ' num2str(chop(triu2vec(corrcoef(attn_rsc_all(2,:),attn_awT_all)),2))])
subplot(2,2,3)
scatter(attn_rsc_all(1,:),attn_vwD_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Vis detect- ' num2str(chop(triu2vec(corrcoef(attn_rsc_all(1,:),attn_vwD_all)),2))])
subplot(2,2,4)
scatter(attn_rsc_all(2,:),attn_awD_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Aud detect- ' num2str(chop(triu2vec(corrcoef(attn_rsc_all(2,:),attn_awD_all)),2))])

nexp = size(dc_noAttn,2);
noAttn_rsc_weight_corr = zeros(2,2,nexp);
noAttn_rsc_all = [];
noAttn_vwT_all = [];
noAttn_awT_all = [];
noAttn_vwD_all = [];
noAttn_awD_all = [];
for iexp = 1:nexp
    VIx = strcmp(dc_noAttn(iexp).av(1).trOut,'cr');
    VIx = VIx+strcmp(dc_noAttn(iexp).av(1).trOut,'fa');
    AIx = strcmp(dc_noAttn(iexp).av(2).trOut,'cr');
    AIx = AIx+strcmp(dc_noAttn(iexp).av(2).trOut,'fa');
    resp_temp_v = dc_noAttn(iexp).av(1).respAllCells(find(VIx),:);
    resp_temp_a = dc_noAttn(iexp).av(2).respAllCells(find(AIx),:);
    resp_temp_vz = (resp_temp_v-mean(resp_temp_v,1))./std(resp_temp_v,[],1);
    resp_temp_az = (resp_temp_a-mean(resp_temp_a,1))./std(resp_temp_a,[],1);
    resp_corr_v = corrcoef(resp_temp_vz);
    resp_corr_a = corrcoef(resp_temp_az);
    resp_corr_avg = zeros(2,size(resp_temp_v,2));
    for ic = 1:size(resp_temp,2)
        temp = resp_corr_v;
        temp(:,ic) = [];
        resp_corr_avg(1,ic) = mean(temp(ic,:),2);
        temp = resp_corr_a;
        temp(:,ic) = [];
        resp_corr_avg(2,ic) = mean(temp(ic,:),2);
    end
    noAttn_rsc_weight_corr(1,1,iexp) = triu2vec(corrcoef(resp_corr_avg(1,:),dc_noAttn(iexp).av(1).weightTarget)');
    noAttn_rsc_weight_corr(2,1,iexp) = triu2vec(corrcoef(resp_corr_avg(2,:),dc_noAttn(iexp).av(2).weightTarget)');
    noAttn_rsc_weight_corr(1,2,iexp) = triu2vec(corrcoef(resp_corr_avg(1,:),dc_noAttn(iexp).av(1).weightDetect)');
    noAttn_rsc_weight_corr(2,2,iexp) = triu2vec(corrcoef(resp_corr_avg(2,:),dc_noAttn(iexp).av(2).weightDetect)');
    noAttn_rsc_all = [noAttn_rsc_all resp_corr_avg];
    noAttn_vwT_all = [noAttn_vwT_all dc_noAttn(iexp).av(1).weightTarget'];
    noAttn_awT_all = [noAttn_awT_all dc_noAttn(iexp).av(2).weightTarget'];
    noAttn_vwD_all = [noAttn_vwD_all dc_noAttn(iexp).av(1).weightDetect'];
    noAttn_awD_all = [noAttn_awD_all dc_noAttn(iexp).av(2).weightDetect'];
end

figure;
subplot(2,2,1)
scatter(noAttn_rsc_all(1,:),noAttn_vwT_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Vis target- ' num2str(chop(triu2vec(corrcoef(noAttn_rsc_all(1,:),noAttn_vwT_all)),2))])
subplot(2,2,2)
scatter(noAttn_rsc_all(2,:),noAttn_awT_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Aud target- ' num2str(chop(triu2vec(corrcoef(noAttn_rsc_all(2,:),noAttn_awT_all)),2))])
subplot(2,2,3)
scatter(noAttn_rsc_all(1,:),noAttn_vwD_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Vis detect- ' num2str(chop(triu2vec(corrcoef(noAttn_rsc_all(1,:),noAttn_vwD_all)),2))])
subplot(2,2,4)
scatter(noAttn_rsc_all(2,:),noAttn_awD_all)
xlim([-0.4 0.4])
ylim([-3 3])
axis square
title(['Aud detect- ' num2str(chop(triu2vec(corrcoef(noAttn_rsc_all(2,:),noAttn_awD_all)),2))])

ind_h = [];
ind_m = [];
ind_cr = [];
ind_fa = [];
for iexp = 1:nexp
resp_temp_v = dc_attn(iexp).av(1).respAllCells;
resp_temp_a = dc_attn(iexp).av(2).respAllCells;
resp_temp_av = [resp_temp_v; resp_temp_a];
resp_temp_avz = (resp_temp_av-mean(resp_temp_av,1))./std(resp_temp_av,[],1);
[coeff,score,latent,tsquared,explained] = pca(resp_temp_avz);
trOut = [dc_attn(iexp).av(1).trOut dc_attn(iexp).av(2).trOut];
av = [zeros(1,length(dc_attn(iexp).av(1).trOut)) ones(1,length(dc_attn(iexp).av(2).trOut))];
ind_v_cr = intersect(find(av==0),find(strcmp(trOut,'cr')));
ind_v_fa = intersect(find(av==0),find(strcmp(trOut,'fa')));
ind_v_h = intersect(find(av==0),find(strcmp(trOut,'h')));
ind_v_m = intersect(find(av==0),find(strcmp(trOut,'m')));
ind_a_cr = intersect(find(av),find(strcmp(trOut,'cr')));
ind_a_fa = intersect(find(av),find(strcmp(trOut,'fa')));
ind_a_h = intersect(find(av),find(strcmp(trOut,'h')));
ind_a_m = intersect(find(av),find(strcmp(trOut,'m')));
ind_h = [ind_h; length(ind_v_h) length(ind_a_h)];
ind_m = [ind_m; length(ind_v_m) length(ind_a_m)];
ind_cr = [ind_cr; length(ind_v_cr) length(ind_a_cr)];
ind_fa = [ind_fa; length(ind_v_fa) length(ind_a_fa)];
figure;
subplot(2,2,1)
scatter(score(ind_v_cr,1),score(ind_v_cr,2))
hold on
scatter(score(ind_v_fa,1),score(ind_v_fa,2))
scatter(score(ind_v_h,1),score(ind_v_h,2))
scatter(score(ind_v_m,1),score(ind_v_m,2))
ylim([-6 6])
xlim([-6 6])

subplot(2,2,2)
scatter(score(ind_a_cr,1),score(ind_a_cr,2))
hold on
scatter(score(ind_a_fa,1),score(ind_a_fa,2))
scatter(score(ind_a_h,1),score(ind_a_h,2))
scatter(score(ind_a_m,1),score(ind_a_m,2))
ylim([-6 6])
xlim([-6 6])

subplot(2,2,3)
errorbar(mean(score(ind_v_cr,1),1),mean(score(ind_v_cr,2),1),std(score(ind_v_cr,2),[],1)./length(ind_v_cr), std(score(ind_v_cr,2),[],1)./length(ind_v_cr),std(score(ind_v_cr,1),[],1)./length(ind_v_cr),std(score(ind_v_cr,1),[],1)./length(ind_v_cr))
hold on
errorbar(mean(score(ind_v_fa,1),1),mean(score(ind_v_fa,2),1),std(score(ind_v_fa,2),[],1)./length(ind_v_fa), std(score(ind_v_fa,2),[],1)./length(ind_v_fa),std(score(ind_v_fa,1),[],1)./length(ind_v_fa),std(score(ind_v_fa,1),[],1)./length(ind_v_fa))
errorbar(mean(score(ind_v_h,1),1),mean(score(ind_v_h,2),1),std(score(ind_v_h,2),[],1)./length(ind_v_h), std(score(ind_v_h,2),[],1)./length(ind_v_h),std(score(ind_v_h,1),[],1)./length(ind_v_h),std(score(ind_v_h,1),[],1)./length(ind_v_h))
errorbar(mean(score(ind_v_m,1),1),mean(score(ind_v_m,2),1),std(score(ind_v_m,2),[],1)./length(ind_v_m), std(score(ind_v_m,2),[],1)./length(ind_v_m),std(score(ind_v_m,1),[],1)./length(ind_v_m),std(score(ind_v_m,1),[],1)./length(ind_v_m))
ylim([-2 2])
xlim([-2 2])

subplot(2,2,4)
errorbar(mean(score(ind_a_cr,1),1),mean(score(ind_a_cr,2),1),std(score(ind_a_cr,2),[],1)./length(ind_a_cr), std(score(ind_a_cr,2),[],1)./length(ind_a_cr),std(score(ind_a_cr,1),[],1)./length(ind_a_cr),std(score(ind_a_cr,1),[],1)./length(ind_a_cr))
hold on
errorbar(mean(score(ind_a_fa,1),1),mean(score(ind_a_fa,2),1),std(score(ind_a_fa,2),[],1)./length(ind_a_fa), std(score(ind_a_fa,2),[],1)./length(ind_a_fa),std(score(ind_a_fa,1),[],1)./length(ind_a_fa),std(score(ind_a_fa,1),[],1)./length(ind_a_fa))
errorbar(mean(score(ind_a_h,1),1),mean(score(ind_a_h,2),1),std(score(ind_a_h,2),[],1)./length(ind_a_h), std(score(ind_a_h,2),[],1)./length(ind_a_h),std(score(ind_a_h,1),[],1)./length(ind_a_h),std(score(ind_a_h,1),[],1)./length(ind_a_h))
errorbar(mean(score(ind_a_m,1),1),mean(score(ind_a_m,2),1),std(score(ind_a_m,2),[],1)./length(ind_a_m), std(score(ind_a_m,2),[],1)./length(ind_a_m),std(score(ind_a_m,1),[],1)./length(ind_a_m),std(score(ind_a_m,1),[],1)./length(ind_a_m))
ylim([-2 2])
xlim([-2 2])
end

%%
close all
clear all
ds = 'FSAV_attentionV1'; % 'FSAV_V1_100ms_naive'  'FSAV_V1_naive_GCaMP6m'  'FSAV_attentionV1'   'FSAV_attentionV1_noAttn'
analysisDate = '201016'; % attn= '201016'  noAttn= '201016'  naive= '201021'
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];
load(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));

start = 1;
nmouse = size(mouse,2);
for imouse = 1:nmouse
    nexp = size(mouse(imouse).expt,2);
    for iexp = 1:nexp
        for iav = 1:2
            trResp = [];
            trOut_act = [];
            trStim_act = [];
            trOut_dig = [];
            trStim_dig = [];
            de = mouse(imouse).expt(iexp).av(iav);
            for i = 2:4
                if i<4
                    trResp = [trResp squeeze(mean(de.align(i).respTC(respwin,:,:),1)-mean(de.align(i).respTC(basewin_0,:,:),1))];
                else
                    trResp = [trResp squeeze( mean(de.align(i).respTC(respwin_target,:,:),1)-mean(de.align(i).respTC(basewin_0_target,:,:),1))];
                end
                if i == 2
                    trOut_act = [trOut_act de.align(i).outcome];
                    trStim_act = [trStim_act zeros(size(de.align(i).outcome))];
                    trOut_dig = [trOut_dig ones(size(de.align(i).outcome))];
                    trStim_dig = [trStim_dig zeros(size(de.align(i).outcome))];
                end
                if i == 3
                    temp = cell(size(de.align(i).outcome));
                    temp(:) = {'cr'};
                    trOut_act = [trOut_act temp];
                    trStim_act = [trStim_act zeros(size(de.align(i).outcome))];
                    trOut_dig = [trOut_dig zeros(size(de.align(i).outcome))];
                    trStim_dig = [trStim_dig zeros(size(de.align(i).outcome))];
                end
                if i == 4
                    trOut_act = [trOut_act de.align(i).outcome];
                    if iav == 1
                        trStim_act = [trStim_act de.align(i).ori];
                    else
                        trStim_act = [trStim_act de.align(i).amp];
                    end
                    temp = zeros(size(de.align(i).outcome));
                    temp(find(strcmp(de.align(i).outcome, 'success'))) = 1;
                    trOut_dig = [trOut_dig temp];
                    trStim_dig = [trStim_dig ones(size(de.align(i).outcome))];
                end
            end

            nCells = size(trResp,1);
            stims = unique(trStim_act);
            nstim = length(stims);
            trResp_zall = zscore(trResp);
            trResp_z = zeros(size(trResp));
            for i = 1:nstim
                ind= find(trStim_act == stims(i));
                trResp_z(:,ind) = zscore(trResp(:,ind));
            end
            trResp_corr = corrcoef(trResp_z');
            trResp_corravg = zeros(1,nCells);
            for i = 1:nCells
                temp = trResp_corr(i,:);
                temp(i) = [];
                trResp_corravg(i) = mean(temp);
            end
            
            trStim_dig_z = trStim_dig-mean(trStim_dig);
            trOut_dig_z = trOut_dig-mean(trOut_dig);
            [mdl_stim pctcorr_stim] = fitrsvm_withcrossval(trResp_zall',trStim_dig_z);
            [mdl_choice pctcorr_choice] = fitrsvm_withcrossval(trResp_zall',trOut_dig_z');
            [mdl_choice_dist pctcorr_choice_dist] = fitrsvm_withcrossval(trResp_zall(:,find(trStim_dig==0))',trOut_dig_z(find(trStim_dig==0))');
            
            lm_stimchoice = fitlm(mdl_stim.Beta,mdl_choice.Beta);
            lm_stimcorr = fitlm(trResp_corravg,mdl_stim.Beta);
            lm_choicecorr = fitlm(trResp_corravg,mdl_choice.Beta);
            
            attn(start).av(iav).lm_stimchoice = lm_stimchoice;
            attn(start).av(iav).lm_stimcorr = lm_stimcorr;
            attn(start).av(iav).lm_choicecorr = lm_choicecorr;
            attn(start).av(iav).mdl_stim = mdl_stim;
            attn(start).av(iav).mdl_choice = mdl_choice;
            attn(start).av(iav).mdl_choice_dist = mdl_choice_dist;
            attn(start).av(iav).pctcorr_stim = pctcorr_stim;
            attn(start).av(iav).pctcorr_choice = pctcorr_choice;
            attn(start).av(iav).pctcorr_choice_dist = pctcorr_choice_dist;
            attn(start).av(iav).trResp = trResp;
            attn(start).av(iav).trStim_act = trStim_act;
            attn(start).av(iav).trStim_dig = trStim_dig;
            attn(start).av(iav).trOut_dig = trOut_dig;
            attn(start).av(iav).trOut_act = trOut_act;
            attn(start).av(iav).trResp_corravg = trResp_corravg;
        end
        attn(start).lm_audvis_stim = fitlm(attn(start).av(1).mdl_stim.Beta,attn(start).av(2).mdl_stim.Beta); 
        attn(start).lm_audvis_choice = fitlm(attn(start).av(1).mdl_choice.Beta,attn(start).av(2).mdl_choice.Beta); 
        attn(start).lm_audvis_choice_dist = fitlm(attn(start).av(1).mdl_choice_dist.Beta,attn(start).av(2).mdl_choice_dist.Beta); 
        attn(start).pred_AV_stim = predict(attn(start).av(2).mdl_stim,attn(start).av(1).trResp')+mean(attn(start).av(2).trStim_dig);
        attn(start).pred_AV_choice = predict(attn(start).av(2).mdl_choice,attn(start).av(1).trResp')+mean(attn(start).av(2).trOut_dig);
        attn(start).pred_AV_choice_dist = predict(attn(start).av(2).mdl_choice_dist,attn(start).av(1).trResp(:,find(attn(start).av(1).trOut_dig==0))')+mean(attn(start).av(2).trOut_dig(find(attn(start).av(2).trOut_dig==0)));
        attn(start).pred_VA_stim = predict(attn(start).av(1).mdl_stim,attn(start).av(2).trResp')+mean(attn(start).av(1).trStim_dig);
        attn(start).pred_VA_choice = predict(attn(start).av(1).mdl_choice,attn(start).av(2).trResp')+mean(attn(start).av(1).trOut_dig);
        attn(start).pred_VA_choice_dist = predict(attn(start).av(1).mdl_choice_dist,attn(start).av(2).trResp(:,find(attn(start).av(2).trOut_dig==0))')+mean(attn(start).av(1).trOut_dig(find(attn(start).av(1).trOut_dig==0)));
        start = start + 1;
    end
end
save('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\attn_struct.mat','attn')
%%
clear all
ds = 'FSAV_attentionV1_noAttn'; % 'FSAV_V1_100ms_naive'  'FSAV_V1_naive_GCaMP6m'  'FSAV_attentionV1'   'FSAV_attentionV1_noAttn'
analysisDate = '201016'; % attn= '201016'  noAttn= '201016'  naive= '201021'
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];
load(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));

start = 1;
nmouse = size(mouse,2);
for imouse = 1:nmouse
    nexp = size(mouse(imouse).expt,2);
    for iexp = 1:nexp
        for iav = 1:2
            trResp = [];
            trOut_act = [];
            trStim_act = [];
            trOut_dig = [];
            trStim_dig = [];
            de = mouse(imouse).expt(iexp).av(iav);
            for i = 2:4
                if i<4
                    trResp = [trResp squeeze(mean(de.align(i).respTC(respwin,:,:),1)-mean(de.align(i).respTC(basewin_0,:,:),1))];
                else
                    trResp = [trResp squeeze( mean(de.align(i).respTC(respwin_target,:,:),1)-mean(de.align(i).respTC(basewin_0_target,:,:),1))];
                end
                if i == 2
                    trOut_act = [trOut_act de.align(i).outcome];
                    trStim_act = [trStim_act zeros(size(de.align(i).outcome))];
                    trOut_dig = [trOut_dig ones(size(de.align(i).outcome))];
                    trStim_dig = [trStim_dig zeros(size(de.align(i).outcome))];
                end
                if i == 3
                    temp = cell(size(de.align(i).outcome));
                    temp(:) = {'cr'};
                    trOut_act = [trOut_act temp];
                    trStim_act = [trStim_act zeros(size(de.align(i).outcome))];
                    trOut_dig = [trOut_dig zeros(size(de.align(i).outcome))];
                    trStim_dig = [trStim_dig zeros(size(de.align(i).outcome))];
                end
                if i == 4
                    trOut_act = [trOut_act de.align(i).outcome];
                    if iav == 1
                        trStim_act = [trStim_act de.align(i).ori];
                    else
                        trStim_act = [trStim_act de.align(i).amp];
                    end
                    temp = zeros(size(de.align(i).outcome));
                    temp(find(strcmp(de.align(i).outcome, 'success'))) = 1;
                    trOut_dig = [trOut_dig temp];
                    trStim_dig = [trStim_dig ones(size(de.align(i).outcome))];
                end
            end

            nCells = size(trResp,1);
            stims = unique(trStim_act);
            nstim = length(stims);
            trResp_z = zeros(size(trResp));
            for i = 1:nstim
                ind= find(trStim_act == stims(i));
                trResp_z(:,ind) = (trResp(:,ind)-mean(trResp(:,ind),2))./std(trResp(:,ind),[],2);
            end
            trResp_corr = corrcoef(trResp_z');
            trResp_corravg = zeros(1,nCells);
            for i = 1:nCells
                temp = trResp_corr(i,:);
                temp(i) = [];
                trResp_corravg(i) = mean(temp);
            end
            
            trStim_dig_z = trStim_dig-mean(trStim_dig);
            trOut_dig_z = trOut_dig-mean(trOut_dig);
            mdl_stim = fitrsvm(trResp',trStim_dig_z);
            mdl_choice = fitrsvm(trResp',trOut_dig_z);
            mdl_choice_dist = fitrsvm(trResp(:,find(trStim_dig==0))',trOut_dig_z(find(trStim_dig==0)));
            pred_stim = predict(mdl_stim,trResp')+mean(trStim_dig_z);
            pred_choice = predict(mdl_choice,trResp')+mean(trOut_dig_z);
            pred_choice_dist = predict(mdl_choice_dist,trResp(:,find(trStim_dig==0))')+mean(trOut_dig_z(trStim_dig==0));

            lm_stimchoice = fitlm(mdl_stim.Beta,mdl_choice.Beta);
            lm_stimcorr = fitlm(trResp_corravg,mdl_stim.Beta);
            lm_choicecorr = fitlm(trResp_corravg,mdl_choice.Beta);
            noattn(start).av(iav).lm_stimchoice = lm_stimchoice;
            noattn(start).av(iav).lm_stimcorr = lm_stimcorr;
            noattn(start).av(iav).lm_choicecorr = lm_choicecorr;
            noattn(start).av(iav).mdl_stim = mdl_stim;
            noattn(start).av(iav).mdl_choice = mdl_choice;
            noattn(start).av(iav).mdl_choice_dist = mdl_choice_dist;
            noattn(start).av(iav).pred_stim = pred_stim;
            noattn(start).av(iav).pred_choice = pred_choice;
            noattn(start).av(iav).pred_choice_dist = pred_choice_dist;
            noattn(start).av(iav).trResp = trResp;
            noattn(start).av(iav).trStim_act = trStim_act;
            noattn(start).av(iav).trStim_dig = trStim_dig;
            noattn(start).av(iav).trOut_dig = trOut_dig;
            noattn(start).av(iav).trOut_act = trOut_act;
            noattn(start).av(iav).trResp_corravg = trResp_corravg;
        end
        noattn(start).lm_audvis_stim = fitlm(noattn(start).av(1).mdl_stim.Beta,noattn(start).av(2).mdl_stim.Beta); 
        noattn(start).lm_audvis_choice = fitlm(noattn(start).av(1).mdl_choice.Beta,noattn(start).av(2).mdl_choice.Beta); 
        noattn(start).lm_audvis_choice_dist = fitlm(noattn(start).av(1).mdl_choice_dist.Beta,noattn(start).av(2).mdl_choice_dist.Beta); 
        noattn(start).pred_AV_stim = predict(noattn(start).av(2).mdl_stim,noattn(start).av(1).trResp')+mean(noattn(start).av(2).trStim_dig);
        noattn(start).pred_AV_choice = predict(noattn(start).av(2).mdl_choice,noattn(start).av(1).trResp')+mean(noattn(start).av(2).trOut_dig);
        noattn(start).pred_AV_choice_dist = predict(noattn(start).av(2).mdl_choice_dist,noattn(start).av(1).trResp(:,find(noattn(start).av(1).trOut_dig==0))')+mean(noattn(start).av(2).trOut_dig(find(noattn(start).av(2).trOut_dig==0)));
        noattn(start).pred_VA_stim = predict(noattn(start).av(1).mdl_stim,noattn(start).av(2).trResp')+mean(noattn(start).av(1).trStim_dig);
        noattn(start).pred_VA_choice = predict(noattn(start).av(1).mdl_choice,noattn(start).av(2).trResp')+mean(noattn(start).av(1).trOut_dig);
        noattn(start).pred_VA_choice_dist = predict(noattn(start).av(1).mdl_choice_dist,noattn(start).av(2).trResp(:,find(noattn(start).av(2).trOut_dig==0))')+mean(noattn(start).av(1).trOut_dig(find(noattn(start).av(1).trOut_dig==0)));
        start = start + 1;
    end
end
save('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\noattn_struct.mat','noattn')

%%
clear all
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\noattn_struct.mat')
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\attn_struct.mat')

figure;
subplot(1,3,1)
hist(attn(1).av(1).pred_stim,[-0.5:.1:1.5])
xlim([-0.5 1.5])
title('Stim')
subplot(1,3,2)
hist(attn(1).av(1).pred_choice,[-0.5:.1:1.5])
xlim([-0.5 1.5])
title('Choice')
subplot(1,3,3)
hist(attn(1).av(1).pred_choice_dist,[-0.5:.1:1.5])
xlim([-0.5 1.5])
title('Dist')
n = size(attn,2);
for i = 1:n 
    for ii = 1:2
        pctCorr_stim(i,ii) = mean(attn(i).av(ii).pred_stim>0.5 == attn(i).av(ii).trStim_dig');
        pctCorr_choice(i,ii) = mean(attn(i).av(ii).pred_choice>0.5 == attn(i).av(ii).trOut_dig');
        pctCorr_choice_dist(i,ii) = mean(attn(i).av(ii).pred_choice_dist>0.5 == attn(i).av(ii).trOut_dig(find(attn(i).av(ii).trStim_dig==0))');
    end
end
figure;
subplot(2,3,1)
plot(pctCorr_stim','k')
title('Stim')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)
subplot(2,3,2)
plot(pctCorr_choice','k')
title('Choice')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)
subplot(2,3,3)
plot(pctCorr_choice_dist','k')
title('Distractors')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)

n2 = size(noattn,2);
for i = 1:n2
    for ii = 1:2
        pctCorr_stim_n(i,ii) = mean(noattn(i).av(ii).pred_stim>0.5 == noattn(i).av(ii).trStim_dig');
        pctCorr_choice_n(i,ii) = mean(noattn(i).av(ii).pred_choice>0.5 == noattn(i).av(ii).trOut_dig');
        pctCorr_choice_dist_n(i,ii) = mean(noattn(i).av(ii).pred_choice_dist>0.5 == noattn(i).av(ii).trOut_dig(find(noattn(i).av(ii).trStim_dig==0))');
    end
end

subplot(2,3,4)
plot(pctCorr_stim_n','k')
title('Stim')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)
subplot(2,3,5)
plot(pctCorr_choice_n','k')
title('Choice')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)
subplot(2,3,6)
plot(pctCorr_choice_dist_n','k')
title('Distractors')
set(gca,'XTick',1:2,'XTickLabel',{'Vis','Aud'})
xlim([0 3])
ylim([0.5 1])
hline(0.55)


n = size(attn,2);
for i = 1:n
    for ii = 1:2
        SCh(i,ii,1) = attn(i).av(ii).lm_stimchoice.Coefficients.Estimate(end);
        SCo(i,ii,1) = attn(i).av(ii).lm_stimcorr.Coefficients.Estimate(end);
        ChCo(i,ii,1) = attn(i).av(ii).lm_choicecorr.Coefficients.Estimate(end);
        SCh(i,ii,2) = attn(i).av(ii).lm_stimchoice.Rsquared.Ordinary;
        SCo(i,ii,2) = attn(i).av(ii).lm_stimcorr.Rsquared.Ordinary;
        ChCo(i,ii,2) = attn(i).av(ii).lm_choicecorr.Rsquared.Ordinary;
    end
end
figure; 
subplot(2,3,1)
plot(SCh(:,:,1)','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
ylabel('Slope')
title('Stimulus-Choice')
[h p] = ttest(SCh(:,1,1),SCh(:,2,1));
hold on
text(1.5,.1, num2str(chop(p,2)))
text(.1,.9, 'Attn')
subplot(2,3,2)
plot(SCo(:,:,2)','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
[h p] = ttest(SCo(:,1,2),SCo(:,2,2));
title('Stimulus-Corr')
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')
subplot(2,3,3)
plot(ChCo(:,:,2)','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
[h p] = ttest(ChCo(:,1,2),ChCo(:,2,2));
title('Choice-Corr')
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')

n = size(noattn,2);
for i = 1:n
    for ii = 1:2
        SCh_n(i,ii,1) = noattn(i).av(ii).lm_stimchoice.Coefficients.Estimate(end);
        SCo_n(i,ii,1) = noattn(i).av(ii).lm_stimcorr.Coefficients.Estimate(end);
        ChCo_n(i,ii,1) = noattn(i).av(ii).lm_choicecorr.Coefficients.Estimate(end);
        SCh_n(i,ii,2) = noattn(i).av(ii).lm_stimchoice.Rsquared.Ordinary;
        SCo_n(i,ii,2) = noattn(i).av(ii).lm_stimcorr.Rsquared.Ordinary;
        ChCo_n(i,ii,2) = noattn(i).av(ii).lm_choicecorr.Rsquared.Ordinary;
    end
end

subplot(2,3,4)
plot(SCh_n(:,:,1)','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
ylabel('Slope')
[h p] = ttest(SCh_n(:,1,1),SCh_n(:,2,1));
hold on
text(1.5,.1, num2str(chop(p,2)))
text(.1,.9, 'No Attn')
subplot(2,3,5)
plot(SCo_n(:,:,2)','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
[h p] = ttest(SCo_n(:,1,2),SCo_n(:,2,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')
subplot(2,3,6)
plot(ChCo_n(:,:,2)','k')
ylim([0 1])
xlim([0 3])

set(gca, 'Xtick',1:2,'XTickLabel',{'Vis', 'Aud'})
[h p] = ttest(ChCo_n(:,1,2),ChCo_n(:,2,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\VisAud_weightComps.pdf','-dpdf');

figure;
subplot(2,2,1)
plot([SCo(:,1,2) ChCo(:,1,2)]','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Stim', 'Choice'})
[h p] = ttest(SCo(:,1,2),ChCo(:,1,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
text(.1,.9, 'Attn')
ylabel('R2')
title('Visual')
subplot(2,2,2)
plot([SCo(:,2,2) ChCo(:,2,2)]','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Stim', 'Choice'})
[h p] = ttest(SCo(:,2,2),ChCo(:,2,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')
title('Auditory')

subplot(2,2,3)
plot([SCo_n(:,1,2) ChCo_n(:,1,2)]','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Stim', 'Choice'})
[h p] = ttest(SCo_n(:,1,2),ChCo_n(:,1,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
text(.1,.9, 'No Attn')
ylabel('R2')
subplot(2,2,4)
plot([SCo_n(:,2,2) ChCo_n(:,2,2)]','k')
ylim([0 1])
xlim([0 3])
set(gca, 'Xtick',1:2,'XTickLabel',{'Stim', 'Choice'})
[h p] = ttest(SCo_n(:,2,2),ChCo_n(:,2,2));
hold on
text(1.5,.8, num2str(chop(p,2)))
ylabel('R2')
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\StimChoice_weightComps.pdf','-dpdf');

%%
figure
n = size(attn,2);
tit_str = {'Visual','Auditory'};
for i = 1:n
    for ii = 1:2
        subplot(2,2,ii)
        scatter(attn(i).av(ii).mdl_stim.Beta,attn(i).av(ii).mdl_choice.Beta,'ok')
        hold on
        if i == n
            xlabel('Stim Weight')
            ylabel('Choice Weight')
            title(tit_str{ii})
            ylim([-3 3])
            xlim([-4 4])
            refline(1)
            text(-3,2,'Attn')
        end
    end
end
n = size(noattn,2);
for i = 1:n
    for ii = 1:2
        subplot(2,2,ii+2)
        scatter(noattn(i).av(ii).mdl_stim.Beta,noattn(i).av(ii).mdl_choice.Beta,'ok')
        hold on
        if i == n
            xlabel('Stim Weight')
            ylabel('Choice Weight')
            ylim([-3 3])
            xlim([-4 4])
            refline(1)
            text(-3,2,'No Attn')
        end
    end
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\StimChoice_weightScatters.pdf','-dpdf');
%% auditory-visual
tit_str = {'Stim','Choice','Choice-D'};
n = size(attn,2);
for i = 1:n
    AVS(i,1) = attn(i).lm_audvis_stim.Coefficients.Estimate(end);
    AVC(i,1) = attn(i).lm_audvis_choice.Coefficients.Estimate(end);
    AVCD(i,1) = attn(i).lm_audvis_choice_dist.Coefficients.Estimate(end);
    AVS(i,2) = attn(i).lm_audvis_stim.Rsquared.Ordinary;
    AVC(i,2) = attn(i).lm_audvis_choice.Rsquared.Ordinary;
    AVCD(i,2) = attn(i).lm_audvis_choice_dist.Rsquared.Ordinary;
    AVS(i,3) = attn(i).lm_audvis_stim.Coefficients.pValue(end)<0.05;
    AVC(i,3) = attn(i).lm_audvis_choice.Coefficients.pValue(end)<0.05;
    AVCD(i,3) = attn(i).lm_audvis_choice_dist.Coefficients.pValue(end)<0.05;
end
n2 = size(noattn,2);
for i = 1:n2
    AVS_n(i,1) = noattn(i).lm_audvis_stim.Coefficients.Estimate(end);
    AVC_n(i,1) = noattn(i).lm_audvis_choice.Coefficients.Estimate(end);
    AVCD_n(i,1) = noattn(i).lm_audvis_choice_dist.Coefficients.Estimate(end);
    AVS_n(i,2) = noattn(i).lm_audvis_stim.Rsquared.Ordinary;
    AVC_n(i,2) = noattn(i).lm_audvis_choice.Rsquared.Ordinary;
    AVCD_n(i,2) = noattn(i).lm_audvis_choice_dist.Rsquared.Ordinary;
    AVS_n(i,3) = noattn(i).lm_audvis_stim.Coefficients.pValue(end)<0.05;
    AVC_n(i,3) = noattn(i).lm_audvis_choice.Coefficients.pValue(end)<0.05;
    AVCD_n(i,3) = noattn(i).lm_audvis_choice_dist.Coefficients.pValue(end)<0.05;
end
subplot(2,2,1)
scatter(1.*ones(n,1),AVS(:,1))
hold on
text(1,0.9,num2str(sum(AVS(:,3))))
scatter(2.*ones(n,1),AVC(:,1))
text(2,0.9,num2str(sum(AVC(:,3))))
scatter(3.*ones(n,1),AVCD(:,1))
text(3,0.9,num2str(sum(AVCD(:,3))))
xlim([0 4])
ylim([-.1 1])
ylabel('Slope')
set(gca,'XTick',1:3,'XTicklabel',tit_str)
subplot(2,2,2)
scatter(1.*ones(n,1),AVS(:,2))
hold on
scatter(2.*ones(n,1),AVC(:,2))
scatter(3.*ones(n,1),AVCD(:,2))
xlim([0 4])
ylim([0 1])
ylabel('R2')
set(gca,'XTick',1:3,'XTicklabel',tit_str)
subplot(2,2,3)
scatter(1.*ones(n2,1),AVS_n(:,1))
hold on
text(1,0.9,num2str(sum(AVS_n(:,3))))
scatter(2.*ones(n2,1),AVC_n(:,1))
text(2,0.9,num2str(sum(AVC_n(:,3))))
scatter(3.*ones(n2,1),AVCD_n(:,1))
text(3,0.9,num2str(sum(AVCD_n(:,3))))
xlim([0 4])
ylim([-.1 1])
ylabel('Slope')
set(gca,'XTick',1:3,'XTicklabel',tit_str)
subplot(2,2,4)
scatter(1.*ones(n2,1),AVS_n(:,2))
hold on
scatter(2.*ones(n2,1),AVC_n(:,2))
scatter(3.*ones(n2,1),AVCD_n(:,2))
xlim([0 4])
ylim([0 1])
ylabel('Slope')
set(gca,'XTick',1:3,'XTicklabel',tit_str)


figure;
n = size(attn,2);

for i = 1:n
    subplot(2,3,1)
    scatter(attn(i).av(1).mdl_stim.Beta,attn(i).av(2).mdl_stim.Beta,'ok')
    hold on
    subplot(2,3,2)
    scatter(attn(i).av(1).mdl_choice.Beta,attn(i).av(2).mdl_choice.Beta,'ok')
    hold on
    subplot(2,3,3)
    scatter(attn(i).av(1).mdl_choice_dist.Beta,attn(i).av(2).mdl_choice_dist.Beta,'ok')
    hold on
    if i == n
        for ii = 1:3
            subplot(2,3,ii)
            xlabel('Visual Weight')
            ylabel('Auditory Weight')
            ylim([-3 3])
            xlim([-4 4])
            refline(1)
            text(-3,2,'Attn')
            title(tit_str{ii})
        end
    end
end
n = size(noattn,2);
for i = 1:n
    subplot(2,3,4)
    scatter(noattn(i).av(1).mdl_stim.Beta,noattn(i).av(2).mdl_stim.Beta,'ok')
    hold on
    subplot(2,3,5)
    scatter(noattn(i).av(1).mdl_choice.Beta,noattn(i).av(2).mdl_choice.Beta,'ok')
    hold on
    subplot(2,3,6)
    scatter(noattn(i).av(1).mdl_choice_dist.Beta,noattn(i).av(2).mdl_choice_dist.Beta,'ok')
    hold on
    if i == n
        for ii = 1:3
            subplot(2,3,3+ii)
            xlabel('Visual Weight')
            ylabel('Auditory Weight')
            ylim([-3 3])
            xlim([-4 4])
            refline(1)
            text(-3,2,'No Attn')
            title(tit_str{ii})
        end
    end
end
print('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\FromAshley\VisAud_weightScatters.pdf','-dpdf');

