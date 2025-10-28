close all
clear all
clc

% mouse = strvcat('i1412','i2585','i1406');
% area = 'V1';
% date = strvcat('241129', '241202', '241202');
% ImgFolder = strvcat({'003'},{'002'},{'003'});

mouse = strvcat('i1414','i1423');
area = 'V1';
date = strvcat('251009', '251009');
ImgFolder = [{'003'},{'002'}];

nrun = size(ImgFolder,2);

stim_set = 'Grat5_Img5_300ms';

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey';
fn_out = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg_LG', stim_set);
if ~exist(fn_out)
    mkdir(fn_out)
end 
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
    run_str = catRunName(ImgFolder(iexp), nrun);
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

    stim_seq = celleqel2mat_padded(input.tstimOne);
    R1_cell_trial = resp_mat(:,:,1);
    R2_cell_trial = resp_mat(:,:,2);
    stims = unique(stim_seq);
    nStim = length(stims);
    [nCells nTrials] = size(R1_cell_trial);
    
    h_stim = zeros(nStim,nCells);
    p_stim = zeros(nStim,nCells);
    R1_avg = zeros(nStim,nCells);
    R2_avg = zeros(nStim,nCells);
    R1_snr = zeros(nStim,nCells);
    for istim = 1:nStim
        ind = find(stim_seq == stims(istim));
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
%     [h1 p1] = ttest(R1_cell_trial,0,'tail','right','alpha',0.05,'Dim',2);
%     good_ind = unique([find(h1)' find(sum(h_stim,1))]);
%     resp_mat = zeros(nCells,nTrials,2);
%     resp_mat(:,:,1) = R1_cell_trial;
%     resp_mat(:,:,2) = R2_cell_trial;
% 
%     if ~exist(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)]))
%         mkdir(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)]))
%     end
%     save(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)], [date(iexp,:) '_' mouse(iexp,:) '_respData.mat']),'resp_mat','good_ind','h1','h_stim','stim_seq')
% end

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
    
    [max_val_sf pref_sf] = max(R1_avg_resp(sf_stim,:),[],1);
    [max_val_nat pref_nat] = max(R1_avg_resp(nat_stim,:),[],1);
    max_snr_sf = indOnly(R1_snr_resp(sf_stim,:)',pref_sf')';
    max_snr_nat = indOnly(R1_snr_resp(nat_stim,:)',pref_nat')';

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

%% 241123 analysis
%outpn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Grants\Adaptation R01\AdaptationR01_Dec2024';
nCells = size(R1_avg_resp_all,2);
R1_avg_resp_all_nan = R1_avg_resp_all;
R1_avg_resp_all_nan(find(h_stim_all==0)) = nan;
R1_avg_resp_all_nan(find(R1_avg_resp_all_nan<min_resp)) = nan;
norm_resp_sf = mean(R2_avg_resp_all(4,:)./R1_avg_resp_all_nan(4,:),1,'omitnan'); 
R1_resp_sf = mean(R1_avg_resp_all_nan(4,:),1,'omitnan'); 
norm_resp_nat = zeros(1,nCells);
R1_resp_nat = zeros(1,nCells);
for i = 1:nCells
    norm_resp_nat(i) = R2_avg_resp_all(pref_nat_all(i)+nsf,i)./R1_avg_resp_all_nan(pref_nat_all(i)+nsf,i);
    R1_resp_nat(i) = R1_avg_resp_all_nan(pref_nat_all(i)+nsf,i);
end

ind_sf = find(~isnan(norm_resp_sf));
ind_nat = find(~isnan(norm_resp_nat));
n_sf = sum(~isnan(norm_resp_sf));
n_nat = sum(~isnan(norm_resp_nat));

figure;
subplot(2,2,1)
errorbar([1 2], [mean(norm_resp_sf,2,'omitnan') mean(norm_resp_nat,2,'omitnan')],[std(norm_resp_sf,[],2,'omitnan')./sqrt(n_sf) std(norm_resp_nat,[],2,'omitnan')./sqrt(n_nat)])
ylim([0 1])
xlim([0 3])
ylabel('Norm dF/F')
[h p] = ttest2(norm_resp_sf,norm_resp_nat);
title(['p = ' num2str(p)])
subplot(2,2,2)
errorbar([1 2], [mean(R1_resp_sf,2,'omitnan') mean(R1_resp_nat,2,'omitnan')],[std(R1_resp_sf,[],2,'omitnan')./sqrt(n_sf) std(R1_resp_nat,[],2,'omitnan')./sqrt(n_nat)])
ylim([0 0.1])
xlim([0 3])
ylabel('R1 dF/F')
[h p] = ttest2(R1_resp_sf,R1_resp_nat);
title(['p = ' num2str(p)])
ind_sub_sf = find(R1_resp_sf>0.05 & R1_resp_sf<0.5);
n_sub_sf = length(ind_sub_sf);
ind_sub_nat = find(R1_resp_nat>0.05 & R1_resp_nat<0.5);
n_sub_nat = length(ind_sub_nat);
subplot(2,2,3)
errorbar([1 2], [mean(norm_resp_sf(ind_sub_sf),2,'omitnan') mean(norm_resp_nat(ind_sub_nat),2,'omitnan')],[std(norm_resp_sf(ind_sub_sf),[],2,'omitnan')./sqrt(n_sub_sf) std(norm_resp_nat(ind_sub_nat),[],2,'omitnan')./sqrt(n_sub_nat)])
ylim([0 1])
xlim([0 3])
ylabel('Norm dF/F')
[h p] = ttest2(norm_resp_sf(ind_sub_sf),norm_resp_nat(ind_sub_nat));
title(['p = ' num2str(p)])
subplot(2,2,4)
errorbar([1 2], [mean(R1_resp_sf(ind_sub_sf),2,'omitnan') mean(R1_resp_nat(ind_sub_nat),2,'omitnan')],[std(R1_resp_sf(ind_sub_sf),[],2,'omitnan')./sqrt(n_sub_sf) std(R1_resp_nat(ind_sub_nat),[],2,'omitnan')./sqrt(n_sub_nat)])
ylim([0 0.1])
xlim([0 3])
ylabel('R1 dF/F')
[h p] = ttest2(R1_resp_sf(ind_sub_sf),R1_resp_nat(ind_sub_nat));
title(['p = ' num2str(p)])
sgtitle('All cells resp to 0.16 OR any nat image')
print(fullfile(fn_out,'NatImgVGratingAdapt_RespEither.pdf'),'-dpdf')

figure;
subplot(2,2,1)
ind_match = intersect(ind_sf,ind_nat);
ind_match_n = length(ind_match);
errorbar([1 2], [mean(norm_resp_sf(ind_match),2,'omitnan') mean(norm_resp_nat(ind_match),2,'omitnan')],[std(norm_resp_sf(ind_match),[],2,'omitnan')./sqrt(ind_match_n) std(norm_resp_nat(ind_match),[],2,'omitnan')./sqrt(ind_match_n)])
ylim([0 1])
xlim([0 3])
ylabel('Norm dF/F')
[h p] = ttest(norm_resp_sf(ind_match),norm_resp_nat(ind_match));
title(['p = ' num2str(p)])
subplot(2,2,2)
errorbar([1 2], [mean(R1_resp_sf(ind_match),2,'omitnan') mean(R1_resp_nat(ind_match),2,'omitnan')],[std(R1_resp_sf(ind_match),[],2,'omitnan')./sqrt(ind_match_n) std(R1_resp_nat(ind_match),[],2,'omitnan')./sqrt(ind_match_n)])
ylim([0 0.1])
xlim([0 3])
ylabel('R1 dF/F')
[h p] = ttest(R1_resp_sf(ind_match),R1_resp_nat(ind_match));
title(['p = ' num2str(p)])
ind_sub_match = intersect(ind_match,intersect(ind_sub_sf,ind_sub_nat));
n_sub_match = length(ind_sub_match);
subplot(2,2,3)
errorbar([1 2], [mean(norm_resp_sf(ind_sub_match),2,'omitnan') mean(norm_resp_nat(ind_sub_match),2,'omitnan')],[std(norm_resp_sf(ind_sub_match),[],2,'omitnan')./sqrt(n_sub_match) std(norm_resp_nat(ind_sub_match),[],2,'omitnan')./sqrt(n_sub_match)])
ylim([0 1])
xlim([0 3])
ylabel('Norm dF/F')
[h p] = ttest(norm_resp_sf(ind_sub_match),norm_resp_nat(ind_sub_match));
title(['p = ' num2str(p)])
subplot(2,2,4)
errorbar([1 2], [mean(R1_resp_sf(ind_sub_match),2,'omitnan') mean(R1_resp_nat(ind_sub_match),2,'omitnan')],[std(R1_resp_sf(ind_sub_match),[],2,'omitnan')./sqrt(n_sub_match) std(R1_resp_nat(ind_sub_match),[],2,'omitnan')./sqrt(n_sub_match)])
ylim([0 0.1])
xlim([0 3])
ylabel('R1 dF/F')
[h p] = ttest(R1_resp_sf(ind_sub_match),R1_resp_nat(ind_sub_match));
title(['p = ' num2str(p)])
sgtitle('All cells resp to 0.16 AND any nat image')

print(fullfile(fn_out,'NatImgVGratingAdapt_RespBoth.pdf'),'-dpdf')
%% 
Adapt_avg_resp_grating_mean = cell(1,nsf+1);
Adapt_avg_resp_natimg_mean = cell(1,nsf+1);
Adapt_avg_resp_grating_max = cell(1,nsf+1);
Adapt_avg_resp_natimg_max = cell(1,nsf+1);
Adapt_avg_resp_grating_mean_all = [];
Adapt_avg_resp_natimg_mean_all = [];
Adapt_avg_resp_grating_max_all = [];
Adapt_avg_resp_natimg_max_all = [];

figure;
for i = 1:nsf
    stim_ind = sf_stim(i);
    ind_use = intersect(find(h_stim_all(stim_ind,:)),intersect(find(pref_sf_all==i),find(sum(h_stim_all(nat_stim,:),1))));
    Adapt_avg_resp_grating_max{i} = mean(Adapt_avg_resp_all(stim_ind,ind_use),1,'omitnan');
    Adapt_avg_resp_grating_mean{i} = mean(Adapt_avg_resp_all(sf_stim,ind_use),1,'omitnan');
    for iCell = 1:length(ind_use)
        iC = ind_use(iCell);
        Adapt_avg_resp_natimg_max{i} = [Adapt_avg_resp_natimg_max{i} Adapt_avg_resp_all(nat_stim(pref_nat_all(iC)),iC)];
    end
    Adapt_avg_resp_natimg_mean{i} = mean(Adapt_avg_resp_all(nat_stim,ind_use),1,'omitnan');
    Adapt_avg_resp_grating_max_all = [Adapt_avg_resp_grating_max_all; Adapt_avg_resp_all(stim_ind,ind_use)' i.*ones(length(ind_use),1)];
    Adapt_avg_resp_grating_mean_all = [Adapt_avg_resp_grating_mean_all; mean(Adapt_avg_resp_all(sf_stim,ind_use),1,'omitnan')' i.*ones(length(ind_use),1)];
    Adapt_avg_resp_natimg_mean_all = [Adapt_avg_resp_natimg_mean_all; mean(Adapt_avg_resp_all(nat_stim,ind_use),1,'omitnan')' i.*ones(length(ind_use),1)];
    for iCell = 1:length(ind_use)
        iC = ind_use(iCell);
        Adapt_avg_resp_natimg_max_all = [Adapt_avg_resp_natimg_max_all; Adapt_avg_resp_all(nat_stim(pref_nat_all(iC)),iC) i];
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

ind_use = intersect(find(sum(h_stim_all(sf_stim,:),1)),find(sum(h_stim_all(nat_stim,:),1)));
Adapt_avg_resp_grating_mean{nsf+1} = mean(Adapt_avg_resp_all(sf_stim,ind_use),1,'omitnan');
Adapt_avg_resp_natimg_mean{nsf+1} = mean(Adapt_avg_resp_all(nat_stim,ind_use),1,'omitnan');
subplot(3,3,nsf+1)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_mean{nsf+1},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_mean{nsf+1},'k')
errorbar(1,mean(Adapt_avg_resp_grating_mean{nsf+1},2,'omitnan'),std(Adapt_avg_resp_grating_mean{nsf+1},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_mean{nsf+1}))),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_mean{nsf+1},2,'omitnan'),std(Adapt_avg_resp_natimg_mean{nsf+1},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_mean{nsf+1}))),'or')
[h p] = ttest(Adapt_avg_resp_natimg_mean{nsf+1},Adapt_avg_resp_grating_mean{nsf+1});
title(['All cells avg- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_mean{nsf+1}))) '; p = ' num2str(chop(p,2))])
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})

Adapt_avg_resp_grating_max{nsf+1} = [];
Adapt_avg_resp_natimg_max{nsf+1} = [];
for iCell = 1:length(ind_use)
    iC = ind_use(iCell);
    Adapt_avg_resp_grating_max{nsf+1} = [Adapt_avg_resp_grating_max{nsf+1} Adapt_avg_resp_all(sf_stim(pref_sf_all(iC)),iC)];
    Adapt_avg_resp_natimg_max{nsf+1} = [Adapt_avg_resp_natimg_max{nsf+1} Adapt_avg_resp_all(nat_stim(pref_nat_all(iC)),iC)];
end
subplot(3,3,8)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{nsf+1},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{nsf+1},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_grating_max{nsf+1},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_max{nsf+1}))),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_natimg_max{nsf+1},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1}))),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{nsf+1},Adapt_avg_resp_grating_max{nsf+1});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1}))) '; p = ' num2str(chop(p,2))])
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
save(fullfile(fn_out,[stim_set '_allData.mat']), 'Adapt_avg_resp_grating_mean','Adapt_avg_resp_natimg_mean','Adapt_avg_resp_grating_max', ...
'Adapt_avg_resp_natimg_max','Adapt_avg_resp_grating_mean_all','Adapt_avg_resp_natimg_mean_all','Adapt_avg_resp_grating_max_all', ...
'Adapt_avg_resp_natimg_max_all', 'R1_avg_resp_all','R2_avg_resp_all','Adapt_avg_resp_all','pref_sf_all','pref_nat_all','max_val_sf_all', ...
'max_val_nat_all','h_stim_all','min_resp','min_dist')
print(fullfile(fn_out, [stim_set '_GratImgComp.pdf']),'-dpdf','-fillpage')

%grating/nat pairs
Adapt_avg_resp_grating_nat_max = cell(nsf,nnat);
col_mat = defaultPlotColors(1:max([nnat nsf]));
figure;
[n n2] = subplotn(nsf+1);
for i = 1:nsf
    sf_ind = sf_stim(i);
    for ii = 1:nnat
        nat_ind = nat_stim(ii);
        ind_use = intersect(find(h_stim_all(sf_ind,:)),intersect(find(pref_sf_all==i),intersect(find(h_stim_all(nat_ind,:)),find(pref_nat_all==ii))));
        Adapt_avg_resp_grating_nat_max{i,ii} = [Adapt_avg_resp_all(sf_ind,ind_use); Adapt_avg_resp_all(nat_ind,ind_use)];
        subplot(n,n2,i)
        errorbar(nanmean(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),2),nanmean(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),2), ...
            nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)), ...
            nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),'o')
        hold on;
        subplot(n,n2,nsf+1)
        errorbar(nanmean(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),2),nanmean(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),2), ...
            nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)), ...
            nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_grating_nat_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),'o','Color',col_mat(i,:))
        hold on;
    end
    subplot(n,n2,i)
    xlim([-1 1])
    ylim([-1 1])
    axis square
    h = refline(1);
    set(h, 'Color', 'black')
    title(['SF ' num2str(i)])
    xlabel('Grating Adapt')
    ylabel('Nat img adapt')
end
subplot(n,n2,nsf+1)
xlim([-1 1])
ylim([-1 1])
axis square
h = refline(1);
set(h, 'Color', 'black')
title(['All pairs'])
xlabel('Grating Adapt')
ylabel('Nat img adapt')
print(fullfile(fn_out,[stim_set '_GratImgComp_Allpairs_bySF.pdf']),'-dpdf','-fillpage')

Adapt_avg_resp_nat_grating_max = cell(nnat,nsf);
figure;
[n n2] = subplotn(nnat+1);
for i = 1:nnat
    nat_ind = nat_stim(i);
    for ii = 1:nsf
        sf_ind = sf_stim(ii);
        ind_use = intersect(find(h_stim_all(sf_ind,:)),intersect(find(pref_sf_all==ii),intersect(find(h_stim_all(nat_ind,:)),find(pref_nat_all==i))));
        Adapt_avg_resp_nat_grating_max{i,ii} = [Adapt_avg_resp_all(sf_ind,ind_use); Adapt_avg_resp_all(nat_ind,ind_use)];
        subplot(n,n2,i)
        errorbar(nanmean(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),2),nanmean(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),2), ...
            nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)), ...
            nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),'o')
        hold on;
        subplot(n,n2,nnat+1)
        errorbar(nanmean(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),2),nanmean(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),2), ...
            nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(1,:),[],2)./sqrt(length(ind_use)), ...
            nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),nanstd(Adapt_avg_resp_nat_grating_max{i,ii}(2,:),[],2)./sqrt(length(ind_use)),'o','Color',col_mat(i,:))
        hold on;
    end
    subplot(n,n2,i)
    xlim([-1 1])
    ylim([-1 1])
    axis square
    h = refline(1);
    set(h, 'Color', 'black')
    title(['NatImg ' num2str(i)])
    xlabel('Grating Adapt')
    ylabel('Nat img adapt')
end
subplot(n,n2,nnat+1)
xlim([-1 1])
ylim([-1 1])
axis square
h = refline(1);
set(h, 'Color', 'black')
title(['All pairs'])
xlabel('Grating Adapt')
ylabel('Nat img adapt')
print(fullfile(fn_out,[stim_set '_GratImgComp_Allpairs_byNat.pdf']),'-dpdf','-fillpage')

figure;
Adapt_avg_resp_pair = [];
subplot(1,2,1)
for i = 1:nsf
    temp_avg = [];
    for ii = 1:nnat
        temp_avg = [temp_avg Adapt_avg_resp_grating_nat_max{i,ii}];
    end
    errorbar(nanmean(temp_avg(1,:),2),nanmean(temp_avg(2,:),2),...
        nanstd(temp_avg(1,:),[],2)./sqrt(size(temp_avg,2)),nanstd(temp_avg(1,:),[],2)./sqrt(size(temp_avg,2)),...
        nanstd(temp_avg(2,:),[],2)./sqrt(size(temp_avg,2)),nanstd(temp_avg(2,:),[],2)./sqrt(size(temp_avg,2)),'o','Color',col_mat(i,:))
    hold on
    Adapt_avg_resp_pair = [Adapt_avg_resp_pair temp_avg];
end
xlim([-0.6 0.4])
ylim([-0.6 0.4])
axis square
h = refline(1);
set(h, 'Color', 'black')
title(['All Nat by SF'])
xlabel('Grating Adapt')
ylabel('Nat img adapt')

subplot(1,2,2)
for i = 1:nnat
    temp_avg = [];
    for ii = 1:nsf
        temp_avg = [temp_avg Adapt_avg_resp_nat_grating_max{i,ii}];
    end
    errorbar(nanmean(temp_avg(1,:),2),nanmean(temp_avg(2,:),2),...
        nanstd(temp_avg(1,:),[],2)./sqrt(size(temp_avg,2)),nanstd(temp_avg(1,:),[],2)./sqrt(size(temp_avg,2)),...
        nanstd(temp_avg(2,:),[],2)./sqrt(size(temp_avg,2)),nanstd(temp_avg(2,:),[],2)./sqrt(size(temp_avg,2)),'o','Color',col_mat(i,:))
    hold on
end
xlim([-0.6 0.4])
ylim([-0.6 0.4])
axis square
h = refline(1);
set(h, 'Color', 'black')
title(['All SF by Nat'])
xlabel('Grating Adapt')
ylabel('Nat img adapt')
print(fullfile(fn_out,[stim_set '_GratImgComp_AllPairsAvg.pdf']),'-dpdf','-fillpage')

%violin plots
figure;
for i = 1:nsf
    plot([1 2],[mean(Adapt_avg_resp_grating_max{i},2,'omitnan') mean(Adapt_avg_resp_natimg_max{i},2,'omitnan')])
    hold on
end
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_grating_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_grating_mean{nsf+1}))),'ok' )
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_natimg_max{i},[],2,'omitnan')./sqrt(sum(~isnan(Adapt_avg_resp_natimg_mean{nsf+1}))),'ok' )
xlim([0 3])
set(gca,'Xtick',[1 2],'XTickLabels',{'grating','image'})
ylim([-0.5 0.5])
axis square
print(fullfile(fn_out,[stim_set '_GratImgCompBySF.pdf']),'-dpdf','-fillpage')

figure;
violin(Adapt_avg_resp_pair','facealpha',1)
ylim([-2 2])
set(gca,'XTick',1:2,'XTickLabel',{'Grating','Image'})
print(fullfile(fn_out,[stim_set '_GratImgPairViolin.pdf']),'-dpdf','-fillpage')

figure;
subplot(3,2,1)
swarmchart(Adapt_avg_resp_grating_max_all(:,2),Adapt_avg_resp_grating_max_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_grating_max_all(:,1),Adapt_avg_resp_grating_max_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:nsf
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
xlim([0 nsf+1])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,2)
swarmchart(Adapt_avg_resp_natimg_max_all(:,2),Adapt_avg_resp_natimg_max_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_natimg_max_all(:,1),Adapt_avg_resp_natimg_max_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:nsf
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
xlim([0 nsf+1])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,3)
swarmchart(Adapt_avg_resp_grating_mean_all(:,2),Adapt_avg_resp_grating_mean_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_grating_mean_all(:,1),Adapt_avg_resp_grating_mean_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:nsf
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
xlim([0 nsf+1])
sigstar(starmat,out(sigind,end))
hline(0)

subplot(3,2,4)
swarmchart(Adapt_avg_resp_natimg_mean_all(:,2),Adapt_avg_resp_natimg_mean_all(:,1),'k')
[p tbl stats] = anovan(Adapt_avg_resp_natimg_mean_all(:,1),Adapt_avg_resp_natimg_mean_all(:,2),'display','off');
out = multcompare(stats,'display','off');
hold on
for i = 1:nsf
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
xlim([0 nsf+1])
sigstar(starmat,out(sigind,end))
hline(0)

Adapt_avg_resp_grating_all = [];
subplot(3,2,5)
for i = 1:nsf
    stim_ind = sf_stim(i);
    ind_use = find(h_stim_all(stim_ind,:));
    Adapt_avg_resp_grating_all = [Adapt_avg_resp_grating_all; Adapt_avg_resp_all(stim_ind,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart(i.*ones(length(ind_use),1),Adapt_avg_resp_all(stim_ind,ind_use)','k')
    hold on
    errorbar(i,mean(Adapt_avg_resp_all(stim_ind,ind_use),2), std(Adapt_avg_resp_all(stim_ind,ind_use),[],2)./length(ind_use),'or');
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
xlim([0 nsf+1])
sigstar(starmat,out(sigind,end))
hline(0)
print(fullfile(fn_out,[stim_set '_PrefSF.pdf']),'-dpdf','-fillpage')

% amplitude match
ind_use = intersect(find(sum(h_stim_all(sf_stim,:),1)),find(sum(h_stim_all(nat_stim,:),1)));
figure;
subplot(2,2,1)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{nsf+1},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{nsf+1},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_grating_max{nsf+1},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_max{nsf+1})),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_natimg_max{nsf+1},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1})),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{nsf+1},Adapt_avg_resp_grating_max{nsf+1});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1}))) '; p = ' num2str(chop(p,2))])
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

ind_use = intersect(find(sum(h_stim_all(sf_stim,:),1)),find(sum(h_stim_all(nat_stim,:),1)));
ind_use_amp = intersect(find(max_val_sf_all(ind_use)>0.05), intersect(find(max_val_sf_all(ind_use)<0.15), intersect(find(max_val_nat_all(ind_use)>0.05), find(max_val_nat_all(ind_use)<0.15))));
[h_amp p_amp] = ttest(max_val_sf_all(ind_use(ind_use_amp)),max_val_nat_all(ind_use(ind_use_amp)));
[h_ad p_ad] = ttest(Adapt_avg_resp_grating_max{nsf+1}(ind_use_amp),Adapt_avg_resp_natimg_max{nsf+1}(ind_use_amp));
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
swarmchart(ones(size(max_val_sf_all(ind_use(ind_use_amp)))), Adapt_avg_resp_grating_max{nsf+1}(ind_use_amp),'k')
hold on
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1}(ind_use_amp),2),std(Adapt_avg_resp_grating_max{nsf+1}(ind_use_amp),[],2)./sqrt(length(ind_use_amp)),'or')
swarmchart(2.*ones(size(max_val_nat_all(ind_use(ind_use_amp)))), Adapt_avg_resp_natimg_max{nsf+1}(ind_use_amp),'k')
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1}(ind_use_amp),2),std(Adapt_avg_resp_natimg_max{nsf+1}(ind_use_amp),[],2)./sqrt(length(ind_use_amp)),'or')
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['Matched: n = ' num2str(length(ind_use_amp)) '; p = ' num2str(chop(p_ad,2))])
suptitle('Amplitude match')
print(fullfile(fn_out,[stim_set '_AmpMatch.pdf']),'-dpdf','-fillpage')

% snr match
ind_use = intersect(find(sum(h_stim_all(sf_stim,:),1)),find(sum(h_stim_all(nat_stim,:),1)));
figure;
subplot(2,2,1)
swarmchart(ones(1,length(ind_use)),Adapt_avg_resp_grating_max{nsf+1},'k')
hold on
swarmchart(2.*ones(1,length(ind_use)),Adapt_avg_resp_natimg_max{nsf+1},'k')
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_grating_max{nsf+1},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_grating_max{nsf+1})),'or')
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1},2,'omitnan'),std(Adapt_avg_resp_natimg_max{nsf+1},[],2,'omitnan')./sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1})),'or')
[h p] = ttest(Adapt_avg_resp_natimg_max{nsf+1},Adapt_avg_resp_grating_max{nsf+1});
title(['All cells max- n = ' num2str(sum(~isnan(Adapt_avg_resp_natimg_max{nsf+1}))) '; p = ' num2str(chop(p,2))])
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

ind_use = intersect(find(sum(h_stim_all(sf_stim,:),1)),find(sum(h_stim_all(nat_stim,:),1)));
ind_use_snr = intersect(find(max_snr_sf_all(ind_use)>0.4), intersect(find(max_snr_sf_all(ind_use)<0.9), intersect(find(max_snr_nat_all(ind_use)>0.4), find(max_snr_nat_all(ind_use)<0.9))));
[h_amp p_snr] = ttest(max_snr_sf_all(ind_use(ind_use_snr)),max_snr_nat_all(ind_use(ind_use_snr)));
[h_ad p_ad] = ttest(Adapt_avg_resp_grating_max{nsf+1}(ind_use_snr),Adapt_avg_resp_natimg_max{nsf+1}(ind_use_snr));
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
swarmchart(ones(size(max_val_sf_all(ind_use(ind_use_snr)))), Adapt_avg_resp_grating_max{nsf+1}(ind_use_snr),'k')
hold on
errorbar(1,mean(Adapt_avg_resp_grating_max{nsf+1}(ind_use_snr),2),std(Adapt_avg_resp_grating_max{nsf+1}(ind_use_snr),[],2)./sqrt(length(ind_use_snr)),'or')
swarmchart(2.*ones(size(max_val_nat_all(ind_use(ind_use_snr)))), Adapt_avg_resp_natimg_max{nsf+1}(ind_use_snr),'k')
errorbar(2,mean(Adapt_avg_resp_natimg_max{nsf+1}(ind_use_snr),2),std(Adapt_avg_resp_natimg_max{nsf+1}(ind_use_snr),[],2)./sqrt(length(ind_use_snr)),'or')
ylim([-2 2])
ylabel('Adapt Ix')
set(gca,'Xtick',[1 2], 'XtickLabel', {'grating', 'image'})
title(['Matched: n = ' num2str(length(ind_use_snr)) '; p = ' num2str(chop(p_ad,2))])
suptitle('SNR match')
print(fullfile(fn_out,[stim_set '_SNRMatch.pdf']),'-dpdf','-fillpage')

%sf pref of nat image response
figure;
Adapt_avg_resp_natimg_all = [];
subplot(3,1,1)
for i = 1:nnat
    stim_use = nat_stim(i);
    ind_use = find(h_stim_all(stim_use,:));
    Adapt_avg_resp_natimg_all = [Adapt_avg_resp_natimg_all; Adapt_avg_resp_all(stim_use,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart((i).*ones(length(ind_use),1),Adapt_avg_resp_all(stim_use,ind_use)','k')
    hold on
    errorbar(i,mean(Adapt_avg_resp_all(stim_use,ind_use),2), std(Adapt_avg_resp_all(stim_use,ind_use),[],2)./length(ind_use),'or');
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
for i = 1:nnat
    stim_use = nat_stim(i);
    ind_use = intersect(find(nat_stim(pref_nat_all)==stim_use),find(h_stim_all(stim_use,:)));
    Adapt_avg_resp_natimg_pref = [Adapt_avg_resp_natimg_pref; Adapt_avg_resp_all(stim_use,ind_use)' i.*ones(length(ind_use),1)];
    swarmchart(i.*ones(length(ind_use),1),Adapt_avg_resp_all(nat_ind,ind_use)','k')
    hold on
    errorbar(i,mean(Adapt_avg_resp_all(stim_use,ind_use),2), std(Adapt_avg_resp_all(stim_use,ind_use),[],2)./length(ind_use),'or');
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
for i = 1:nnat
    stim_use = nat_stim(i);
    ind_use = intersect(find(nat_stim(pref_nat_all)==stim_use),find(h_stim_all(stim_use,:)));
    if length(ind_use)>0
        cdfplot(pref_sf_all(ind_use))
        hold on
    end
end
xlabel('Preferred SF')
ylabel('Fraction of cells')
title('Pref image')
subplot(3,2,6)
for i = 1:nnat
    stim_use = nat_stim(i);
    ind_use = intersect(find(nat_stim(pref_nat_all)==stim_use),find(h_stim_all(stim_use,:)));
    errorbar(mean(pref_sf_all(ind_use),2), mean(Adapt_avg_resp_all(stim_use,ind_use),2),...
        std(Adapt_avg_resp_all(stim_use,ind_use),[],2)./sqrt(length(ind_use)), std(Adapt_avg_resp_all(stim_use,ind_use),[],2)./sqrt(length(ind_use)),...
        std(pref_sf_all(ind_use),[],2)./sqrt(length(ind_use)),std(pref_sf_all(ind_use),[],2)./sqrt(length(ind_use)),'o')
    hold on
end
xlabel('Preferred SF')
ylabel('Adapt Ix- Pref image')
ylim([-1 1.2])
print(fullfile(fn_out,[stim_set '_PrefImg.pdf']),'-dpdf','-fillpage')

%% eye-tracking
eye_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\';
mworks_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Analysis\2P\';
mouse = strvcat('i1380','i1381','i1386','i1374','i1387','i1375');
area = 'V1';
date = strvcat('230330','230404','230406','230411','230418','230425');
fn_out = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg', stim_set);


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
fn_out = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg', stim_set);

nexp = length(mouse);
sfs = [0.02 0.04 0.08 0.16 0.32 0.64];
min_resp = 0.02;
doEyeDist = 0;
min_dist = 4;
PCbySF = nan(100,nsf,nexp);
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
        
        [max_val_sf pref_sf] = max(R1_avg_resp(sf_stim,:),[],1);
        [max_val_nat pref_nat] = max(R1_avg_resp(nat_stim,:),[],1);
        
        R1_cell_trial_z = zscore(R1_cell_trial);
        natimg_tr = find(stim_seq>5);
        ind_use = intersect(find(sum(h_stim(sf_stim,:))),find(sum(h_stim(nat_stim,:))));
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
        for i = 1:nsf
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
        Adapt_maxNat = indOnly(Adapt_avg_resp',pref_nat'+nsf);
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
for i = 1:nsf
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
fn_out = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg',stim_set);

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