data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\';
pupil_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Analysis\2P\';
mouse = strvcat('i1380','i1381','i1386','i1374','i1387','i1375');
area = 'V1';
date = strvcat('230330','230404','230406','230411','230418','230425');
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg\Grat6_Img8';

nexp = length(mouse);

min_resp = 0.02;
doEyeDist = 0;
min_dist = 4;

R1_avg_resp_all = [];
R2_avg_resp_all = [];
R1_snr_resp_all = [];
Adapt_avg_resp_all = [];
pref_sf_all = [];
pref_nat_all = [];
max_val_sf_all = [];
max_val_nat_all = [];
max_snr_sf_all = [];
max_snr_nat_all = [];
h_stim_all = [];
ntrialperstim = zeros(nexp,14);

for iexp = 1:nexp
    fprintf([mouse(iexp,:) ' ' date(iexp,:) '\n'])
    load(fullfile(data_pn,[area '_' mouse(iexp,:) '_' date(iexp,:) '_cellpose'],'trace_trial_stim.mat'))
    load(fullfile(pupil_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004'],[date(iexp,:) '_' mouse(iexp,:) '_pupil.mat']))
   
    if doEyeDist
        ind_dist = find(centroid.dist<=min_dist);
    else
        ind_dist = 1:length(stim_seq);
    end

    stims = unique(stim_seq);
    nStim = length(stims);
    nCells = size(R1_cell_trial,1);
    
    h_stim = zeros(nStim,nCells);
    p_stim = zeros(nStim,nCells);
    R1_avg = zeros(nStim,nCells);
    R2_avg = zeros(nStim,nCells);
    R1_snr = zeros(nStim,nCells);
    for istim = 1:nStim
        ind = intersect(ind_dist,find(stim_seq == stims(istim)));
        ntrialperstim(iexp,istim) = length(ind);
        if length(ind)>40
            R1_avg(istim,:) = mean(R1_cell_trial(:,ind),2,'omitnan');
            R2_avg(istim,:) = mean(R2_cell_trial(:,ind),2,'omitnan');
            R1_snr(istim,:) =  R1_avg(istim,:)./std(R1_cell_trial(:,ind),[],2,'omitnan')';
            [h_stim(istim,:) p_stim(istim,:)] = ttest(R1_cell_trial(:,ind),0,'tail','right','alpha',0.05/(nStim-1),'Dim',2);
        else
            R1_avg(istim,:) = nan(1,nCells);
            R2_avg(istim,:) = nan(1,nCells);
            R1_snr(istim,:) = nan(1,nCells);
            h_stim(istim,:) = nan(1,nCells);
            p_stim(istim,:) = nan(1,nCells);
        end
    end

    h_stim(find(R1_avg<min_resp)) = 0; %minimum response amp
    nonsigind = find(h_stim==0);
    R1_avg_resp = R1_avg;
    R2_avg_resp = R2_avg;
    R1_avg_resp(nonsigind) = nan;
    R2_avg_resp(nonsigind) = nan;
    R1_avg_resp = reshape(R1_avg_resp,size(R1_avg));
    R2_avg_resp = reshape(R2_avg_resp,size(R2_avg));
    R1_snr_resp = R1_snr;
    R1_snr_resp(nonsigind) = nan;
    R1_snr_resp = reshape(R1_snr_resp,size(R1_snr));

    Adapt_avg_resp = (R2_avg_resp-R1_avg_resp)./R1_avg_resp;
    
    [max_val_sf pref_sf] = max(R1_avg_resp(1:6,:),[],1);
    [max_val_nat pref_nat] = max(R1_avg_resp(7:end,:),[],1);
    max_snr_sf = indOnly(R1_snr_resp(1:6,:)',pref_sf')';
    max_snr_nat = indOnly(R1_snr_resp(7:end,:)',pref_nat')';

    R1_avg_resp_all = [R1_avg_resp_all R1_avg_resp];
    R2_avg_resp_all = [R2_avg_resp_all R2_avg_resp];
    R1_snr_resp_all = [R1_snr_resp_all R1_snr_resp];
    Adapt_avg_resp_all = [Adapt_avg_resp_all Adapt_avg_resp];
    pref_sf_all = [pref_sf_all pref_sf];
    pref_nat_all = [pref_nat_all pref_nat];
    max_val_sf_all = [max_val_sf_all max_val_sf];
    max_val_nat_all = [max_val_nat_all max_val_nat];
    max_snr_sf_all = [max_snr_sf_all max_snr_sf];
    max_snr_nat_all = [max_snr_nat_all max_snr_nat];
    h_stim_all = [h_stim_all h_stim];
end

Adapt_avg_resp_grating_mean = cell(1,7);
Adapt_avg_resp_natimg_mean = cell(1,7);
Adapt_avg_resp_grating_max = cell(1,7);
Adapt_avg_resp_natimg_max = cell(1,7);
Adapt_avg_resp_grating_mean_all = [];
Adapt_avg_resp_natimg_mean_all = [];
Adapt_avg_resp_grating_max_all = [];
Adapt_avg_resp_natimg_max_all = [];

figure;
for i = 1:6
    ind_use = intersect(find(h_stim_all(i,:)),intersect(find(pref_sf_all==i),find(sum(h_stim_all(7:end,:),1))));
    Adapt_avg_resp_grating_max{i} = mean(Adapt_avg_resp_all(i,ind_use),1,'omitnan');
    Adapt_avg_resp_grating_mean{i} = mean(Adapt_avg_resp_all(1:6,ind_use),1,'omitnan');
    for iCell = 1:length(ind_use)
        iC = ind_use(iCell);
        Adapt_avg_resp_natimg_max{i} = [Adapt_avg_resp_natimg_max{i} Adapt_avg_resp_all(pref_nat_all(iC)+6,iC)];
    end
    Adapt_avg_resp_natimg_mean{i} = mean(Adapt_avg_resp_all(7:end,ind_use),1,'omitnan');
    Adapt_avg_resp_grating_max_all = [Adapt_avg_resp_grating_max_all; Adapt_avg_resp_all(i,ind_use)' i.*ones(length(ind_use),1)];
    Adapt_avg_resp_grating_mean_all = [Adapt_avg_resp_grating_mean_all; mean(Adapt_avg_resp_all(1:6,ind_use),1,'omitnan')' i.*ones(length(ind_use),1)];
    Adapt_avg_resp_natimg_mean_all = [Adapt_avg_resp_natimg_mean_all; mean(Adapt_avg_resp_all(7:end,ind_use),1,'omitnan')' i.*ones(length(ind_use),1)];
    for iCell = 1:length(ind_use)
        iC = ind_use(iCell);
        Adapt_avg_resp_natimg_max_all = [Adapt_avg_resp_natimg_max_all; Adapt_avg_resp_all(pref_nat_all(iC)+6,iC) i];
    end
    subplot(3,3,i)
    swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{i},'k')
    hold on
    swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{i},'k')
    errorbar(1,mean(Adapt_avg_resp_grating_max{i},2,'omitnan'),std(Adapt_avg_resp_grating_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_max{i}))),'or')
    errorbar(2,mean(Adapt_avg_resp_natimg_max{i},2,'omitnan'),std(Adapt_avg_resp_natimg_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_max{i}))),'or')
    [h p] = ttest(Adapt_avg_resp_natimg_max{i},Adapt_avg_resp_grating_max{i});
    title(['Pref SF' num2str(i) ';- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{i}))) '; p = ' num2str(chop(p,2))])
    ylim([-2 2])
    ylabel('Adapt Ix')
    set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
end

ind_use = intersect(find(sum(h_stim_all(1:6,:),1)),find(sum(h_stim_all(7:end,:),1)));
Adapt_avg_resp_grating_mean{7} = mean(Adapt_avg_resp_all(1:6,ind_use),1,'omitnan');
Adapt_avg_resp_natimg_mean{7} = mean(Adapt_avg_resp_all(7:end,ind_use),1,'omitnan');
subplot(3,3,7)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_mean{7},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_mean{7},'k')
errorbar(1,mean(Adapt_avg_resp_grating_mean{7},2,'omitnan'),std(Adapt_avg_resp_grating_mean{7},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_mean{7}))),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_mean{7},2,'omitnan'),std(Adapt_avg_resp_natimg_mean{7},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_mean{7}))),'or')
[h p] = ttest(Adapt_avg_resp_natimg_mean{7},Adapt_avg_resp_grating_mean{7});
title(['All cells avg- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_mean{7}))) '; p = ' num2str(chop(p,2))])
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

Adapt_avg_resp_grating_max{7} = [];
Adapt_avg_resp_natimg_max{7} = [];
for iCell = 1:length(ind_use)
    iC = ind_use(iCell);
    Adapt_avg_resp_grating_max{7} = [Adapt_avg_resp_grating_max{7} Adapt_avg_resp_all(pref_sf_all(iC),iC)];
    Adapt_avg_resp_natimg_max{7} = [Adapt_avg_resp_natimg_max{7} Adapt_avg_resp_all(pref_nat_all(iC)+6,iC)];
end
subplot(3,3,8)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{7},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{7},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{7},2,'omitnan'),std(Adapt_avg_resp_grating_max{7},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_max{7}))),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{7},2,'omitnan'),std(Adapt_avg_resp_natimg_max{7},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_max{7}))),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{7},Adapt_avg_resp_grating_max{7});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{7}))) '; p = ' num2str(chop(p,2))])
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

subplot(3,3,9)
swarmchart(ones(1,length(ind_use)),max_val_sf_all(ind_use),'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),max_val_nat_all(ind_use),'k')
errorbar(1,mean(max_val_sf_all(ind_use)),std(max_val_sf_all(ind_use))./sqrt(sum(~isnan(max_val_sf_all(ind_use)))),'or')
errorbar(2,mean(max_val_nat_all(ind_use)),std(max_val_nat_all(ind_use))./sqrt(sum(~isnan(max_val_nat_all(ind_use)))),'or')
[h p] = ttest(max_val_sf_all(ind_use),max_val_nat_all(ind_use));
title(['All cells max resp- n = ' num2str(sum(~isnan(max_val_sf_all(ind_use)))) '; p = ' num2str(chop(p,2))])
ylim([0 .5])
ylabel('dF/F')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
suptitle(['Min resp = ' num2str(min_resp) 'dF/F; MinEyeDist = ' num2str(min_dist) ' deg'])
save(fullfile(fn_out,'Grat6_Img8_allData.mat'), 'Adapt_avg_resp_grating_mean','Adapt_avg_resp_natimg_mean','Adapt_avg_resp_grating_max', ...
'Adapt_avg_resp_natimg_max','Adapt_avg_resp_grating_mean_all','Adapt_avg_resp_natimg_mean_all','Adapt_avg_resp_grating_max_all', ...
'Adapt_avg_resp_natimg_max_all', 'R1_avg_resp_all','R2_avg_resp_all','Adapt_avg_resp_all','pref_sf_all','pref_nat_all','max_val_sf_all', ...
'max_val_nat_all','h_stim_all','min_resp','min_dist')
print(fullfile(fn_out,'Grat6_Img8_GratImgComp.pdf'),'-dpdf','-fillpage')

%violin plots
figure;
subplot(2,2,1)
for i = 1:6
    plot([1 2],[mean(Adapt_avg_resp_grating_max{i},2,'omitnan') mean(Adapt_avg_resp_natimg_max{i},2,'omitnan')])
    hold on
end
errorbar(1,mean(Adapt_avg_resp_grating_max{7},2,'omitnan'),std(Adapt_avg_resp_grating_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_mean{7}))),'ok' )
errorbar(2,mean(Adapt_avg_resp_natimg_max{7},2,'omitnan'),std(Adapt_avg_resp_natimg_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_mean{7}))),'ok' )
xlim([0 3])
set(gca,'Xtick',[1 2],'XTickLabels',{'grating','image'})
ylim([-0.5 0.5])
axis square
print(fullfile(fn_out,'Grat6_Img8_GratImgCompBySF.pdf'),'-dpdf','-fillpage')


figure;
subplot(3,2,1)
swarmchart(Adapt_avg_resp_grating_max_all(:,2),Adapt_avg_resp_grating_max_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_grating_max_all(:,1),Adapt_avg_resp_grating_max_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:6
    errorbar(i,mean(Adapt_avg_resp_grating_max{i},2,'omitnan'),std(Adapt_avg_resp_grating_max{i},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_max{i})),'or')
end
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Gratings- max')
ylabel('Adaptation Ix')
xlabel('SF')
ylim([-2 3])
xlim([0 7])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,2)
swarmchart(Adapt_avg_resp_natimg_max_all(:,2),Adapt_avg_resp_natimg_max_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_natimg_max_all(:,1),Adapt_avg_resp_natimg_max_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:6
    errorbar(i,mean(Adapt_avg_resp_natimg_max{i},2,'omitnan'),std(Adapt_avg_resp_natimg_max{i},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_max{i})),'or')
end
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Natural images- max')
ylabel('Adaptation Ix')
xlabel('SF')
ylim([-2 3])
xlim([0 7])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,3)
swarmchart(Adapt_avg_resp_grating_mean_all(:,2),Adapt_avg_resp_grating_mean_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_grating_mean_all(:,1),Adapt_avg_resp_grating_mean_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:6
    errorbar(i,mean(Adapt_avg_resp_grating_mean{i},2,'omitnan'),std(Adapt_avg_resp_grating_mean{i},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_mean{i})),'or')
end
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Gratings- mean')
ylabel('Adaptation Ix')
xlabel('SF')
ylim([-2 3])
xlim([0 7])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,4)
swarmchart(Adapt_avg_resp_natimg_mean_all(:,2),Adapt_avg_resp_natimg_mean_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_natimg_mean_all(:,1),Adapt_avg_resp_natimg_mean_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:6
    errorbar(i,mean(Adapt_avg_resp_natimg_mean{i},2,'omitnan'),std(Adapt_avg_resp_natimg_mean{i},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_mean{i})),'or')
end
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Natural images- mean')
ylabel('Adaptation Ix')
xlabel('SF')
ylim([-2 3])
xlim([0 7])
sigstar(starmat,out(sigind,end))
hline(0)

Adapt_avg_resp_grating_all = [];
subplot(3,2,5)
for i = 1:6
    ind_use = find(h_stim_all(i,:));
    Adapt_avg_resp_grating_all = [Adapt_avg_resp_grating_all; Adapt_avg_resp_all(i,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart(i.*ones(length(ind_use),1),Adapt_avg_resp_all(i,ind_use)','k')
    hold on
    errorbar(i,mean(Adapt_avg_resp_all(i,ind_use),2), std(Adapt_avg_resp_all(i,ind_use),[],2)./length(ind_use),'or');
end
[p tbl stats] = anovan(Adapt_avg_resp_grating_all(:,1),Adapt_avg_resp_grating_all(:,2),'display','off');
out = multcompare(stats,'display','off');
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Gratings- all')
ylabel('Adaptation Ix')
xlabel('SF')
ylim([-2 3])
xlim([0 7])
sigstar(starmat,out(sigind,end))
hline(0)
print(fullfile(fn_out,'Grat6_Img8_PrefSF.pdf'),'-dpdf','-fillpage')

% amplitude match
ind_use = intersect(find(sum(h_stim_all(1:6,:),1)),find(sum(h_stim_all(7:end,:),1)));
figure;
subplot(2,2,1)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{7},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{7},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{7},2,'omitnan'),std(Adapt_avg_resp_grating_max{7},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_max{7})),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{7},2,'omitnan'),std(Adapt_avg_resp_natimg_max{7},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_max{7})),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{7},Adapt_avg_resp_grating_max{7});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{7}))) '; p = ' num2str(chop(p,2))])
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

subplot(2,2,2)
swarmchart(ones(1,length(ind_use)),max_val_sf_all(ind_use),'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),max_val_nat_all(ind_use),'k')
errorbar(1,mean(max_val_sf_all(ind_use)),std(max_val_sf_all(ind_use))./sum(~isnan(max_val_sf_all(ind_use))),'or')
errorbar(2,mean(max_val_nat_all(ind_use)),std(max_val_nat_all(ind_use))./sum(~isnan(max_val_nat_all(ind_use))),'or')
[h p] = ttest(max_val_sf_all(ind_use),max_val_nat_all(ind_use));
title(['p = ' num2str(chop(p,2))])
ylim([0 .5])
ylabel('dF/F')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

ind_use = intersect(find(sum(h_stim_all(1:6,:),1)),find(sum(h_stim_all(7:end,:),1)));
ind_use_amp = intersect(find(max_val_sf_all(ind_use)>0.05), intersect(find(max_val_sf_all(ind_use)<0.15), intersect(find(max_val_nat_all(ind_use)>0.05), find(max_val_nat_all(ind_use)<0.15))));
[h_amp p_amp] = ttest(max_val_sf_all(ind_use(ind_use_amp)),max_val_nat_all(ind_use(ind_use_amp)));
[h_ad p_ad] = ttest(Adapt_avg_resp_grating_max{7}(ind_use_amp),Adapt_avg_resp_natimg_max{7}(ind_use_amp));
subplot(2,2,4)
swarmchart(ones(size(max_val_sf_all(ind_use(ind_use_amp)))), max_val_sf_all(ind_use(ind_use_amp)),'k')
hold on
errorbar(1,mean(max_val_sf_all(ind_use(ind_use_amp)),2),std(max_val_sf_all(ind_use(ind_use_amp)),[],2)./sqrt(length(ind_use_amp)),'or')
swarmchart(2.*ones(size(max_val_nat_all(ind_use(ind_use_amp)))), max_val_nat_all(ind_use(ind_use_amp)),'k')
errorbar(2,mean(max_val_nat_all(ind_use(ind_use_amp)),2),std(max_val_nat_all(ind_use(ind_use_amp)),[],2)./sqrt(length(ind_use_amp)),'or')
ylim([0 0.5])
ylabel('dF/F')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['p = ' num2str(chop(p_amp,2))])
subplot(2,2,3)
swarmchart(ones(size(max_val_sf_all(ind_use(ind_use_amp)))), Adapt_avg_resp_grating_max{7}(ind_use_amp),'k')
hold on
errorbar(1,mean(Adapt_avg_resp_grating_max{7}(ind_use_amp),2),std(Adapt_avg_resp_grating_max{7}(ind_use_amp),[],2)./sqrt(length(ind_use_amp)),'or')
swarmchart(2.*ones(size(max_val_nat_all(ind_use(ind_use_amp)))), Adapt_avg_resp_natimg_max{7}(ind_use_amp),'k')
errorbar(2,mean(Adapt_avg_resp_natimg_max{7}(ind_use_amp),2),std(Adapt_avg_resp_natimg_max{7}(ind_use_amp),[],2)./sqrt(length(ind_use_amp)),'or')
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['Matched: n = ' num2str(length(ind_use_amp)) '; p = ' num2str(chop(p_ad,2))])
suptitle('Amplitude match')
print(fullfile(fn_out,'Grat6_Img8_AmpMatch.pdf'),'-dpdf','-fillpage')

% snr match
ind_use = intersect(find(sum(h_stim_all(1:6,:),1)),find(sum(h_stim_all(7:end,:),1)));
figure;
subplot(2,2,1)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{7},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{7},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{7},2,'omitnan'),std(Adapt_avg_resp_grating_max{7},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_max{7})),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{7},2,'omitnan'),std(Adapt_avg_resp_natimg_max{7},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_max{7})),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{7},Adapt_avg_resp_grating_max{7});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{7}))) '; p = ' num2str(chop(p,2))])
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

subplot(2,2,2)
swarmchart(ones(1,length(ind_use)),max_snr_sf_all(ind_use),'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),max_snr_nat_all(ind_use),'k')
errorbar(1,mean(max_snr_sf_all(ind_use)),std(max_snr_sf_all(ind_use))./sum(~isnan(max_snr_sf_all(ind_use))),'or')
errorbar(2,mean(max_snr_nat_all(ind_use)),std(max_snr_nat_all(ind_use))./sum(~isnan(max_snr_nat_all(ind_use))),'or')
[h p] = ttest(max_snr_sf_all(ind_use),max_snr_nat_all(ind_use));
title(['p = ' num2str(chop(p,2))])
ylim([0 2])
ylabel('SNR')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

ind_use = intersect(find(sum(h_stim_all(1:6,:),1)),find(sum(h_stim_all(7:end,:),1)));
ind_use_snr = intersect(find(max_snr_sf_all(ind_use)>0.4), intersect(find(max_snr_sf_all(ind_use)<0.9), intersect(find(max_snr_nat_all(ind_use)>0.4), find(max_snr_nat_all(ind_use)<0.9))));
[h_amp p_snr] = ttest(max_snr_sf_all(ind_use(ind_use_snr)),max_snr_nat_all(ind_use(ind_use_snr)));
[h_ad p_ad] = ttest(Adapt_avg_resp_grating_max{7}(ind_use_snr),Adapt_avg_resp_natimg_max{7}(ind_use_snr));
subplot(2,2,4)
swarmchart(ones(size(max_snr_sf_all(ind_use(ind_use_snr)))), max_snr_sf_all(ind_use(ind_use_snr)),'k')
hold on
errorbar(1,mean(max_snr_sf_all(ind_use(ind_use_snr)),2),std(max_snr_sf_all(ind_use(ind_use_snr)),[],2)./sqrt(length(ind_use_snr)),'or')
swarmchart(2.*ones(size(max_snr_nat_all(ind_use(ind_use_snr)))), max_snr_nat_all(ind_use(ind_use_snr)),'k')
errorbar(2,mean(max_snr_nat_all(ind_use(ind_use_snr)),2),std(max_snr_nat_all(ind_use(ind_use_snr)),[],2)./sqrt(length(ind_use_snr)),'or')
ylim([0 2])
ylabel('SNR')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['p = ' num2str(chop(p_snr,2))])
subplot(2,2,3)
swarmchart(ones(size(max_val_sf_all(ind_use(ind_use_snr)))), Adapt_avg_resp_grating_max{7}(ind_use_snr),'k')
hold on
errorbar(1,mean(Adapt_avg_resp_grating_max{7}(ind_use_snr),2),std(Adapt_avg_resp_grating_max{7}(ind_use_snr),[],2)./sqrt(length(ind_use_snr)),'or')
swarmchart(2.*ones(size(max_val_nat_all(ind_use(ind_use_snr)))), Adapt_avg_resp_natimg_max{7}(ind_use_snr),'k')
errorbar(2,mean(Adapt_avg_resp_natimg_max{7}(ind_use_snr),2),std(Adapt_avg_resp_natimg_max{7}(ind_use_snr),[],2)./sqrt(length(ind_use_snr)),'or')
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['Matched: n = ' num2str(length(ind_use_snr)) '; p = ' num2str(chop(p_ad,2))])
suptitle('SNR match')
print(fullfile(fn_out,'Grat6_Img8_SNRMatch.pdf'),'-dpdf','-fillpage')

%sf pref of nat image response
figure;
Adapt_avg_resp_natimg_all = [];
subplot(3,1,1)
for i = 7:14
    ind_use = find(h_stim_all(i,:));
    Adapt_avg_resp_natimg_all = [Adapt_avg_resp_natimg_all; Adapt_avg_resp_all(i,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart((i-6).*ones(length(ind_use),1),Adapt_avg_resp_all(i,ind_use)','k')
    hold on
    errorbar(i-6,mean(Adapt_avg_resp_all(i,ind_use),2), std(Adapt_avg_resp_all(i,ind_use),[],2)./length(ind_use),'or');
end
[p tbl stats] = anovan(Adapt_avg_resp_natimg_all(:,1),Adapt_avg_resp_natimg_all(:,2),'display','off');
out = multcompare(stats,'display','off');
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Natural images- all')
ylabel('Adaptation Ix')
xlabel('Nat image')
ylim([-2 3])
xlim([0 9])
sigstar(starmat,out(sigind,end))
hline(0)

Adapt_avg_resp_natimg_pref = [];
subplot(3,1,2)
for i = 7:14
    ind_use = intersect(find(pref_nat_all==i-6),find(h_stim_all(i,:)));
    Adapt_avg_resp_natimg_pref = [Adapt_avg_resp_natimg_pref; Adapt_avg_resp_all(i,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart((i-6).*ones(length(ind_use),1),Adapt_avg_resp_all(i,ind_use)','k')
    hold on
    errorbar(i-6,mean(Adapt_avg_resp_all(i,ind_use),2), std(Adapt_avg_resp_all(i,ind_use),[],2)./length(ind_use),'or');
end
[p tbl stats] = anovan(Adapt_avg_resp_natimg_pref(:,1),Adapt_avg_resp_natimg_pref(:,2),'display','off');
out = multcompare(stats,'display','off');
sigind = find(out(:,end)<0.05);
sign = length(sigind);
starmat = cell(1,sign);
for i = 1:sign
    starmat{i} = out(sigind(i),1:2);
end
title('Natural images- pref')
ylabel('Adaptation Ix')
xlabel('Nat image')
ylim([-2 3])
xlim([0 9])
sigstar(starmat,out(sigind,end))
hline(0)
subplot(3,2,5)
for i = 7:14
    ind_use = intersect(find(pref_nat_all==i-6),find(h_stim_all(i,:)));
    if length(ind_use)>0
        cdfplot(pref_sf_all(ind_use))
        hold on
    end
end
xlabel('Preferred SF')
ylabel('Fraction of cells')
title('Pref image')
subplot(3,2,6)
for i = 7:14
    ind_use = intersect(find(pref_nat_all==i-6),find(h_stim_all(i,:)));
    errorbar(mean(pref_sf_all(ind_use),2), mean(Adapt_avg_resp_all(i,ind_use),2),...
        std(Adapt_avg_resp_all(i,ind_use),[],2)./sqrt(length(ind_use)), std(Adapt_avg_resp_all(i,ind_use),[],2)./sqrt(length(ind_use)),...
        std(pref_sf_all(ind_use),[],2)./sqrt(length(ind_use)),std(pref_sf_all(ind_use),[],2)./sqrt(length(ind_use)),'o')
    hold on
end
xlabel('Preferred SF')
ylabel('Adapt Ix- Pref image')
ylim([-1 1.2])
print(fullfile(fn_out,'Grat6_Img8_PrefImg.pdf'),'-dpdf','-fillpage')

%% eye-tracking
eye_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\';
mworks_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Analysis\2P\';
mouse = strvcat('i1380','i1381','i1386','i1374','i1387','i1375');
area = 'V1';
date = strvcat('230330','230404','230406','230411','230418','230425');
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg\Grat6_Img8';


for iexp = 1:length(mouse)
    fprintf([mouse(iexp,:) ' ' date(iexp,:) '\n'])
    eyedata = [];
    clear input
    for irun = 1:3
        fprintf(['00' num2str(irun+1) '\n'])
        load(fullfile(eye_pn,mouse(iexp,:),date(iexp,:),['00' num2str(irun+1)],['00' num2str(irun+1) '_000_000_eye.mat']));
        data_temp = squeeze(data);
        load(fullfile(mworks_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-00' num2str(irun+1)], [date(iexp,:) '_' mouse(iexp,:) '_runs-00' num2str(irun+1) '_input.mat']));
        input(irun) = behav_input;
        nframes = behav_input.counterValues{end}(end);
        ntrials = size(behav_input.counterValues,2);
        if irun > 1
            for itrial = 1:ntrials
                input(irun).cStimOneOn{itrial} = input(irun).cStimOneOn{itrial} + nframes;
                input(irun).cStimTwoOn{itrial} = input(irun).cStimTwoOn{itrial} + nframes;
            end
        end
        
        eyedata = cat(3,eyedata,data_temp(:,:,1:nframes));
    end
    input = concatenateStructuresLG(input);
    mkdir(fullfile(mworks_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004']))
    save(fullfile(mworks_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004'],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004_input.mat']),'input')
    [data_crop rect] = cropEyeData(eyedata);
    rad_range = [3 15]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
    Eye_data = extractEyeData(data_crop,rad_range);
    while 1  % till broken out of
        % interactively get clicks
        [X Y selectionType] = getAPoint(gca);
        if isnan(X)
            key = lower(Y);
            switch key
              case char(13) % return
                break;  % out of loop, done
              case 'z' 
              end
            continue
        end
    end
    [rad centroid] = alignEyeData(Eye_data,input);
    ind_use = find(centroid.dist<4);
    save(fullfile(mworks_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004'],[date(iexp,:) '_' mouse(iexp,:) '_pupil.mat']),'rad','centroid','rect','rad_range','ind_use');
    close all
end

%% PC analysis

data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\';
pupil_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Analysis\2P\';
mouse = strvcat('i1380','i1381','i1386','i1374','i1387','i1375');
area = 'V1';
date = strvcat('230330','230404','230406','230411','230418','230425');
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg\Grat6_Img8';

nexp = length(mouse);
sfs = [0.02 0.04 0.08 0.16 0.32 0.64];
min_resp = 0.02;
doEyeDist = 0;
min_dist = 4;
PCbySF = nan(100,6,nexp);
PCbyAdG = nan(100,2,nexp);
PCbyAdN = nan(100,2,nexp);
prefSF_lowAd_G = [];
prefSF_highAd_G = [];
prefSF_lowAd_N = [];
prefSF_highAd_N = [];
figure;
for iexp = 1:nexp
    fprintf([mouse(iexp,:) ' ' date(iexp,:) '\n'])
    load(fullfile(data_pn,[area '_' mouse(iexp,:) '_' date(iexp,:) '_cellpose'],'trace_trial_stim.mat'))
    load(fullfile(pupil_pn,[date(iexp,:) '_' mouse(iexp,:)],[date(iexp,:) '_' mouse(iexp,:) '_runs-002-004'],[date(iexp,:) '_' mouse(iexp,:) '_pupil.mat']))
    
    if doEyeDist
        ind_dist = find(centroid.dist<=min_dist);
    else
        ind_dist = 1:length(stim_seq);
    end

    stims = unique(stim_seq);
    nStim = length(stims);
    nCells = size(R1_cell_trial,1);
    
    h_stim = zeros(nStim,nCells);
    p_stim = zeros(nStim,nCells);
    R1_avg = zeros(nStim,nCells);
    R2_avg = zeros(nStim,nCells);
    for istim = 1:nStim
        ind = intersect(ind_dist,find(stim_seq == stims(istim)));
        ntrialperstim(iexp,istim) = length(ind);
        if length(ind)>40
            R1_avg(istim,:) = mean(R1_cell_trial(:,ind),2,'omitnan');
            R2_avg(istim,:) = mean(R2_cell_trial(:,ind),2,'omitnan');
            [h_stim(istim,:) p_stim(istim,:)] = ttest(R1_cell_trial(:,ind),0,'tail','right','alpha',0.05/(nStim-1),'Dim',2);
        else
            R1_avg(istim,:) = nan(1,nCells);
            R2_avg(istim,:) = nan(1,nCells);
            h_stim(istim,:) = nan(1,nCells);
            p_stim(istim,:) = nan(1,nCells);
        end
    end
    
    if nCells > 50
        h_stim(find(R1_avg<min_resp)) = 0; %minimum response amp
        nonsigind = find(h_stim==0);
        R1_avg_resp = R1_avg;
        R2_avg_resp = R2_avg;
        R1_avg_resp(nonsigind) = nan;
        R2_avg_resp(nonsigind) = nan;
        R1_avg_resp = reshape(R1_avg_resp,size(R1_avg));
        R2_avg_resp = reshape(R2_avg_resp,size(R2_avg));
        
        Adapt_avg_resp = (R2_avg_resp-R1_avg_resp)./R1_avg_resp;
        
        [max_val_sf pref_sf] = max(R1_avg_resp(1:6,:),[],1);
        [max_val_nat pref_nat] = max(R1_avg_resp(7:end,:),[],1);
        
        R1_cell_trial_z = zscore(R1_cell_trial);
        natimg_tr = find(stim_seq>5);
        ind_use = intersect(find(sum(h_stim(1:6,:))),find(sum(h_stim(7:end,:))));
        % ind_low = intersect(ind_use,find(pref_sf<=3));
        % ind_high = intersect(ind_use,find(pref_sf>3));
        % [u s v_lowsf] = pca(R1_cell_trial_z(ind_low,natimg_tr)');
        % [u s v_highsf] = pca(R1_cell_trial_z(ind_high,natimg_tr)');
        % subplot(nexp,3,1+(3*(iexp-1)))
        % plot(cumsum(v_lowsf))
        % hold on
        % plot(cumsum(v_highsf))
        % title('Low vs High SF')
        subplot(nexp,3,1+(3*(iexp-1)))
        for i = 1:6
            ind = intersect(ind_use,find(pref_sf==i));
            [u s v] = pca(R1_cell_trial_z(ind,natimg_tr)');
            plot(cumsum(v))
            hold on
            PCbySF(1:size(v,1),i,iexp) = cumsum(v);
        end
        if iexp == 1
            legend(num2str(sfs'))
        end
        Adapt_maxSF = indOnly(Adapt_avg_resp',pref_sf');
        Adapt_maxNat = indOnly(Adapt_avg_resp',pref_nat'+6);
        medAdapt_SF = median(Adapt_maxSF,1,'omitnan');
        medAdapt_Nat = median(Adapt_maxNat,1,'omitnan');
        ind_low = intersect(ind_use,find(Adapt_maxSF<=medAdapt_SF));
        ind_high = intersect(ind_use,find(Adapt_maxSF>medAdapt_SF));
        prefSF_lowAd_G = [prefSF_lowAd_G pref_sf(ind_low)];
        prefSF_highAd_G = [prefSF_highAd_G pref_sf(ind_high)];
        [u s v_lowad] = pca(R1_cell_trial_z(ind_low,natimg_tr)');
        [u s v_highad] = pca(R1_cell_trial_z(ind_high,natimg_tr)');
        subplot(nexp,3,2+(3*(iexp-1)))
        plot(cumsum(v_lowad))
        hold on
        plot(cumsum(v_highad))
        PCbyAdG(1:size(v_lowad,1),1,iexp) = cumsum(v_lowad);
        PCbyAdG(1:size(v_highad,1),2,iexp) = cumsum(v_highad);
        title(['Strong vs Weak Adapt- Grating median ' num2str(chop(medAdapt_SF,2))])
        ind_low = intersect(ind_use,find(Adapt_maxNat<=medAdapt_Nat));
        ind_high = intersect(ind_use,find(Adapt_maxNat>medAdapt_Nat));
        prefSF_lowAd_N = [prefSF_lowAd_G pref_sf(ind_low)];
        prefSF_highAd_N = [prefSF_highAd_N pref_sf(ind_high)];
        [u s v_lowad] = pca(R1_cell_trial_z(ind_low,natimg_tr)');
        [u s v_highad] = pca(R1_cell_trial_z(ind_high,natimg_tr)');
        subplot(nexp,3,3+(3*(iexp-1)))
        plot(cumsum(v_lowad))
        hold on
        plot(cumsum(v_highad))
        PCbyAdN(1:size(v_lowad,1),1,iexp) = cumsum(v_lowad);
        PCbyAdN(1:size(v_highad,1),2,iexp) = cumsum(v_highad);
        title(['Strong vs Weak Adapt- Image median ' num2str(chop(medAdapt_Nat,2))])
    end
end
figure;
subplot(3,2,1)
PCbySF_norm = PCbySF./max(max(PCbySF,[],1),[],2);
for i = 1:6
    errorbar(1:100,mean(PCbySF_norm(:,i,:),3,'omitnan'),std(PCbySF_norm(:,i,:),[],3,'omitnan')./sqrt(5))
    hold on
end
subplot(3,2,3)
PCbyAd_norm = PCbyAdG./max(max(PCbyAdG,[],1),[],2);
for i = 1:2
    errorbar(1:100,mean(PCbyAd_norm(:,i,:),3,'omitnan'),std(PCbyAd_norm(:,i,:),[],3,'omitnan')./sqrt(5))
    hold on
end
subplot(3,2,5)
PCbyAdN_norm = PCbyAdN./max(max(PCbyAdN,[],1),[],2);
for i = 1:2
    errorbar(1:100,mean(PCbyAdN_norm(:,i,:),3,'omitnan'),std(PCbyAdN_norm(:,i,:),[],3,'omitnan')./sqrt(5))
    hold on
end
subplot(3,2,4)
cdfplot(prefSF_lowAd_G)
hold on
cdfplot(prefSF_highAd_G)
subplot(3,2,6)
cdfplot(prefSF_lowAd_N)
hold on
cdfplot(prefSF_highAd_N)

%% MNR
data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\';
pupil_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Analysis\2P\';
mouse = strvcat('i1380','i1381','i1386','i1374','i1387','i1375');
area = 'V1';
date = strvcat('230330','230404','230406','230411','230418','230425');
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg\Grat6_Img8';

nexp = length(mouse);

for iexp = 1:nexp
    fprintf([mouse(iexp,:) ' ' date(iexp,:) '\n'])
    load(fullfile(data_pn,[area '_' mouse(iexp,:) '_' date(iexp,:) '_cellpose'],'trace_trial_stim.mat'))
    
    natimg_tr = find(stim_seq>5);
    nCells = size(R1_cell_trial,1);
    cellname = cell(1,nCells);
    cellname{1} = 'Image';
    for i = 1:nCells
        cellname{i+1} = num2str(i);
    end
    R1_cell_trial_tab = array2table([stim_seq(natimg_tr)' R1_cell_trial(:,natimg_tr)'], 'VariableNames',cellname);
    MnrMdl = fitmnr(R1_cell_trial(:,natimg_tr)',stim_seq(natimg_tr)');
end