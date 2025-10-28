close all
clear all
clc

% mouse = strvcat('i1412','i2585','i1406');
% area = 'V1';
% date = strvcat('241129', '241202', '241202');
% ImgFolder = strvcat({'003'},{'002'},{'003'});
% 300 ms stim
% mouse = strvcat('i1414','i1423','i1425');
% area = 'V1';
% date = strvcat('251017', '251017', '251018');
% ImgFolder = [{'002'},{'002'},{'002'}];
% stim_set = 'Grat1_Img6_300ms';


mouse = strvcat('i1426','i1414','i1423');
area = 'V1';
date = strvcat('251018', '251019', '251019');
ImgFolder = [{'002'},{'002'},{'002'}];
stim_set = 'Grat1_Img6_1000ms';

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey';
outpn = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Adaptation\SFSummary\NatImg_LG', stim_set);
if ~exist(outpn)
    mkdir(outpn)
end 
nexp = size(mouse,1);



min_resp = 0.02;
doEyeDist = 0;
min_dist = 3;

R1_avg_resp_all = [];
R2_avg_resp_all = [];
R1_snr_resp_all = [];
Adapt_avg_resp_all = [];
R1_avg_eye_resp_all = [];
R2_avg_eye_resp_all = [];
Adapt_avg_eye_resp_all = [];
pref_sf_all = [];
pref_nat_all = [];
max_val_sf_all = [];
max_val_nat_all = [];
max_snr_sf_all = [];
max_snr_nat_all = [];
h_stim_all = [];

for iexp = 1:nexp
    fprintf([mouse(iexp,:) ' ' date(iexp,:) '\n'])
    run_str = catRunName(ImgFolder{iexp}, 1);
    load(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)], [date(iexp,:) '_' mouse(iexp,:) '_' run_str], [date(iexp,:) '_' mouse(iexp,:) '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)], [date(iexp,:) '_' mouse(iexp,:) '_' run_str], [date(iexp,:) '_' mouse(iexp,:) '_' run_str '_TCs.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)], [date(iexp,:) '_' mouse(iexp,:) '_' run_str], [date(iexp,:) '_' mouse(iexp,:) '_' run_str '_input.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date(iexp,:) '_' mouse(iexp,:)], [date(iexp,:) '_' mouse(iexp,:) '_' run_str], [date(iexp,:) '_' mouse(iexp,:) '_' run_str '_eyeAndWheel.mat']))
    
    cStimTwo = cell2mat(input.cStimTwoOn);
    [nFrames, nCells] = size(data_tc);
    nTrials = length(cStimOne);
    tc_one = nan(50,nCells,nTrials);
    tc_two = nan(50,nCells,nTrials);
    for itrial = 1:nTrials
        if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+29)<nFrames
            tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-20:cStimOne(itrial)+29,:);
        end
        if ~isnan(cStimTwo(itrial)) & (cStimTwo(itrial)+29)<nFrames
            tc_two(:,:,itrial) = npSub_tc(cStimTwo(itrial)-20:cStimTwo(itrial)+29,:);
        end
    end
    tc_one_f = mean(tc_one(1:20,:,:));
    tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;
    tc_two_dfof = (tc_two-tc_one_f)./tc_one_f;

    dfof_resp_one = squeeze(mean(tc_one_dfof(resp_win,:,:),1));
    dfof_base_one = squeeze(mean(tc_one_dfof(base_win,:,:),1));
    dfof_resp_two = squeeze(mean(tc_two_dfof(resp_win,:,:),1));
    dfof_base_two = squeeze(mean(tc_two_dfof(base_win,:,:),1));
    
    nStim = nImage + nGrating;

    h_stim = zeros(nStim,nCells);
    p_stim = zeros(nStim,nCells);
    R1_avg = zeros(nStim,nCells);
    R1_snr = zeros(nStim,nCells);
    R2_avg = zeros(nStim,nCells);
    R1_avg_eye = zeros(nStim,nCells);
    R2_avg_eye = zeros(nStim,nCells);

    nat_stim = [];
    start = 1;
    for i = 1:nImage
        ind = intersect(find(tGrating==0), find(tStimOne == stimOne(i)));
        R1_avg(start,:) = (mean(dfof_resp_one(:,ind)-dfof_base_one(:,ind),2,'omitnan'))';
        R1_snr(start,:) = (mean(dfof_resp_one(:,ind)-dfof_base_one(:,ind),2,'omitnan')./std(dfof_resp_one(:,ind)-dfof_base_one(:,ind),[],2,'omitnan'))';
        R2_avg(start,:) = (mean(dfof_resp_two(:,ind)-dfof_base_two(:,ind),2,'omitnan'))';
        [h_stim(start,:) p_stim(start,:)] = ttest(dfof_resp_one(:,ind),dfof_base_one(:,ind),'tail','right','alpha',0.05/(nStim-1),'Dim',2);
        nat_stim = [nat_stim start];
        ind_eye = intersect(find(centroid.dist<min_dist),ind);
        R1_avg_eye(start,:) = (mean(dfof_resp_one(:,ind_eye)-dfof_base_one(:,ind_eye),2,'omitnan'))';
        R2_avg_eye(start,:) = (mean(dfof_resp_two(:,ind_eye)-dfof_base_two(:,ind_eye),2,'omitnan'))';
        start = start+1;
    end
   
    sf_stim = [];
    for i = 1:nGrating
        ind = intersect(find(tGrating==1), find(tGratingSF == gratingSFs(i)));
        R1_avg(start,:) = (mean(dfof_resp_one(:,ind)-dfof_base_one(:,ind),2,'omitnan'))';
        R1_snr(start,:) = (mean(dfof_resp_one(:,ind)-dfof_base_one(:,ind),2,'omitnan')./std(dfof_resp_one(:,ind)-dfof_base_one(:,ind),[],2,'omitnan'))';
        R2_avg(start,:) = (mean(dfof_resp_two(:,ind)-dfof_base_two(:,ind),2,'omitnan'))';
        [h_stim(start,:) p_stim(start,:)] = ttest(dfof_resp_one(:,ind),dfof_base_one(:,ind),'tail','right','alpha',0.05/(nStim-1),'Dim',2);
        ind_eye = intersect(find(centroid.dist<min_dist),ind);
        R1_avg_eye(start,:) = (mean(dfof_resp_one(:,ind_eye)-dfof_base_one(:,ind_eye),2,'omitnan'))';
        R2_avg_eye(start,:) = (mean(dfof_resp_two(:,ind_eye)-dfof_base_two(:,ind_eye),2,'omitnan'))';
        sf_stim = [sf_stim start];
        start = start+1;
    end

    h_stim(find(R1_avg<min_resp)) = 0; %minimum response amp
    nonsigind = find(h_stim==0);
    R1_avg_resp = R1_avg;
    R1_snr_resp = R1_snr;
    R2_avg_resp = R2_avg;
    R1_avg_resp(nonsigind) = nan;
    R1_snr_resp(nonsigind) = nan;
    R2_avg_resp(nonsigind) = nan;
    R1_avg_resp = reshape(R1_avg_resp,size(R1_avg));
    R1_snr_resp = reshape(R1_snr_resp,size(R1_snr));
    R2_avg_resp = reshape(R2_avg_resp,size(R2_avg));
    R1_avg_eye_resp = R1_avg_eye;
    R2_avg_eye_resp = R2_avg_eye;
    R1_avg_eye_resp(nonsigind) = nan;
    R2_avg_eye_resp(nonsigind) = nan;
    R1_avg_eye_resp = reshape(R1_avg_eye_resp,size(R1_avg_eye));
    R2_avg_eye_resp = reshape(R2_avg_eye_resp,size(R2_avg_eye));

    Adapt_avg_resp = R2_avg_resp./R1_avg_resp;
    Adapt_avg_eye_resp = R2_avg_eye_resp./R1_avg_eye_resp;

    [max_val_sf pref_sf] = max(R1_avg_resp(sf_stim,:),[],1);
    [max_val_nat pref_nat] = max(R1_avg_resp(nat_stim,:),[],1);
    max_snr_sf = indOnly(R1_snr_resp(sf_stim,:)',pref_sf')';
    max_snr_nat = indOnly(R1_snr_resp(nat_stim,:)',pref_nat')';

    R1_avg_resp_all = [R1_avg_resp_all R1_avg_resp];
    R2_avg_resp_all = [R2_avg_resp_all R2_avg_resp];
    R1_snr_resp_all = [R1_snr_resp_all R1_snr_resp];
    Adapt_avg_resp_all = [Adapt_avg_resp_all Adapt_avg_resp];
    R1_avg_eye_resp_all = [R1_avg_eye_resp_all R1_avg_eye_resp];
    R2_avg_eye_resp_all = [R2_avg_eye_resp_all R2_avg_eye_resp];
    Adapt_avg_eye_resp_all = [Adapt_avg_eye_resp_all Adapt_avg_eye_resp];
    pref_sf_all = [pref_sf_all pref_sf];
    pref_nat_all = [pref_nat_all pref_nat];
    max_val_sf_all = [max_val_sf_all max_val_sf];
    max_val_nat_all = [max_val_nat_all max_val_nat];
    max_snr_sf_all = [max_snr_sf_all max_snr_sf];
    max_snr_nat_all = [max_snr_nat_all max_snr_nat];
    h_stim_all = [h_stim_all h_stim];
end


%% Grat vs Nat img - preferred
nCells = size(h_stim_all,2);

ind_sf = find(sum(h_stim_all(sf_stim,:),1));
ind_nat = find(sum(h_stim_all(nat_stim,:),1));
ind_both = intersect(ind_sf,ind_nat);

Adapt_resp_pref_sf = zeros(1,nCells);
Adapt_eye_resp_pref_sf = zeros(1,nCells);
Adapt_resp_pref_nat = zeros(1,nCells);
Adapt_eye_resp_pref_nat = zeros(1,nCells);
for iCell = 1:nCells
    Adapt_resp_pref_sf(iCell) = Adapt_avg_resp_all(nImage+pref_sf_all(iCell),iCell);
    Adapt_eye_resp_pref_sf(iCell) = Adapt_avg_eye_resp_all(nImage+pref_sf_all(iCell),iCell);
    Adapt_resp_pref_nat(iCell) = Adapt_avg_resp_all(pref_nat_all(iCell),iCell);
    Adapt_eye_resp_pref_nat(iCell) =  Adapt_avg_eye_resp_all(pref_nat_all(iCell),iCell);
end

figure; 
subplot(2,2,1)
errorbar([1 2],[mean(Adapt_resp_pref_sf(ind_sf)) mean(Adapt_resp_pref_nat(ind_nat))],[std(Adapt_resp_pref_sf(ind_sf))./sqrt(length(ind_sf)) std(Adapt_resp_pref_nat(ind_nat))./sqrt(length(ind_nat))])
[h p] = ttest2(Adapt_resp_pref_sf(ind_sf), Adapt_resp_pref_nat(ind_nat));
title(['Pref - p = ' num2str(chop(p,3))])
xlim([0 3])
ylim([0 1])
subplot(2,2,2)
errorbar([1 2],[mean(Adapt_resp_pref_sf(ind_both)) mean(Adapt_resp_pref_nat(ind_both))],[std(Adapt_resp_pref_sf(ind_both))./sqrt(length(ind_both)) std(Adapt_resp_pref_nat(ind_both))./sqrt(length(ind_both))])
[h p] = ttest(Adapt_resp_pref_sf(ind_both), Adapt_resp_pref_nat(ind_both));
title(['Match - p = ' num2str(chop(p,3))])
xlim([0 3])
ylim([0 1])
subplot(2,2,3)
errorbar([1 2],[mean(Adapt_eye_resp_pref_sf(ind_sf)) mean(Adapt_eye_resp_pref_nat(ind_nat))],[std(Adapt_eye_resp_pref_sf(ind_sf))./sqrt(length(ind_sf)) std(Adapt_eye_resp_pref_nat(ind_nat))./sqrt(length(ind_nat))])
[h p] = ttest2(Adapt_eye_resp_pref_sf(ind_sf), Adapt_eye_resp_pref_nat(ind_nat));
title(['Eye < ' num2str(min_dist) ' deg; Pref - p = ' num2str(chop(p,3))])
xlim([0 3])
ylim([0 1])
subplot(2,2,4)
errorbar([1 2],[mean(Adapt_eye_resp_pref_sf(ind_both)) mean(Adapt_eye_resp_pref_nat(ind_both))],[std(Adapt_eye_resp_pref_sf(ind_both))./sqrt(length(ind_both)) std(Adapt_eye_resp_pref_nat(ind_both))./sqrt(length(ind_both))])
[h p] = ttest(Adapt_eye_resp_pref_sf(ind_both), Adapt_eye_resp_pref_nat(ind_both));
title(['Eye < ' num2str(min_dist) ' deg; Match - p = ' num2str(chop(p,3))])
xlim([0 3])
ylim([0 1])
print(fullfile(outpn,'NatImgVGratingAdapt_RespBoth_eye.pdf'),'-dpdf')
