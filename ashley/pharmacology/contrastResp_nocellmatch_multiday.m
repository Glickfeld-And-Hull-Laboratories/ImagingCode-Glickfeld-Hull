%% which dataset and parameters
clear all
clear global
close all
ds = 'DART_V1_PV_contrast'; %dataset info
dataStructLabels = 'contrastxori';

rc = behavConstsAV; %directories
eval(ds)


%% which days to use?
day1_id = 8;
day2_id = 9;
day3_id = nan;

%% mouse and saving info
mouse = expt(day1_id).mouse;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
fnout = fullfile(fn,'multiday');
mkdir(fnout)

%% load time-courses
data_tag = cell(1,2);
data_notag = cell(1,2);
% day 1
expDate_day1 = expt(day1_id).date;
expTime_day1 = eval(['expt(day1_id).' dataStructLabels '_time']);
expRuns_day1 = eval(['expt(day1_id).' dataStructLabels '_runs']);
for irun = 1:length(expRuns_day1)
    load(fullfile(fn,expDate_day1,expRuns_day1{irun},'timecourses_cells'))
    data_tag{1}{irun} = tc_subnp_red;
    data_notag{1}{irun} = tc_subnp_green;
end
% day 2
expDate_day2 = expt(day2_id).date;
expTime_day2 = eval(['expt(day2_id).' dataStructLabels '_time']);
load(fullfile(fn,expDate_day2,eval(['expt(day2_id).' dataStructLabels '_runs{1}']),'timecourses_cells'))
data_tag{2} = tc_subnp_red;
data_notag{2} = tc_subnp_green;

if ~isnan(day3_id)
    % day 3
    expDate_day3 = expt(day3_id).date;
    expTime_day3 = eval(['expt(day3_id).' dataStructLabels '_time']);
    load(fullfile(fn,expDate_day3,eval(['expt(day3_id).' dataStructLabels '_runs{1}']),'timecourses_cells'))
    data_tag{3} = tc_subnp_red;
    data_notag{3} = tc_subnp_green;
end
%% load mworks data
mw_day1 = [];
for irun = 1:length(expRuns_day1)
    mw_day1{irun} = loadMworksFile(mouse,expDate_day1,expTime_day1{irun});
end
mw_day2 = loadMworksFile(mouse,expDate_day2,expTime_day2{1});
if ~isnan(day3_id)
    mw_day3 = loadMworksFile(mouse,expDate_day3,expTime_day3{1});
    mw = {mw_day1,mw_day2,mw_day3};
else
    mw = {mw_day1,mw_day2};
end
%% params
frameRateHz = expt(day1_id).frame_rate;
nBaseFr = frameRateHz;
nRespFr = frameRateHz;
winlengthFr = round(0.5*frameRateHz);
basewin = (nBaseFr-winlengthFr):nBaseFr;
respwin = nBaseFr+3:(nBaseFr+winlengthFr);
%% day 1 analysis
%tagged
[dff_eaTrial,avgResp_eaTrial] = getEaTrialResp_visstimret(data_tag{1},mw_day1,frameRateHz);

conID = cell2mat(cellfun(@(x) celleqel2mat_padded(x.tGratingContrast),mw_day1,'unif',0));
dirID = cell2mat(cellfun(@(x) celleqel2mat_padded(x.tGratingDirectionDeg),mw_day1,'unif',0));
oriID = nan(size(dirID));
oriID(dirID<180) = dirID(dirID<180);
oriID(dirID>=180) = dirID(dirID>=180)-180;
contrasts = unique(conID);
orientations = unique(oriID);
ncon = length(contrasts);
nori = length(orientations);
[nt,nc] = size(avgResp_eaTrial);

resp_oriXcon = nan(ncon,nori,nc);
respErr_oriXcon = nan(ncon,nori,nc);
respCell_oriXcon = nan(ncon,nori,nc);
for icon = 1:ncon
    for iori = 1:nori
        ind = conID == contrasts(icon) & oriID == orientations(iori);
        resp_oriXcon(icon,iori,:) = mean(avgResp_eaTrial(ind,:),1);
        respErr_oriXcon(icon,iori,:) = ste(avgResp_eaTrial(ind,:),1);
%         if icon == ncon
            b = cell2mat(cellfun(@(x) mean(x(basewin,:)),dff_eaTrial(ind),'unif',0)');
            r = cell2mat(cellfun(@(x) mean(x(respwin,:)),dff_eaTrial(ind),'unif',0)');
            respCell_oriXcon(icon,iori,:) = ttest(r,b,'dim',1,'tail','right','alpha',0.01./nori);
%         end
    end
end


[ori_fit,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(end,:,:))',...
    squeeze(resp_oriXcon(end,:,:))',double(orientations),1);
ori_fit = squeeze(ori_fit(1:180,1,:));
[ori_fit_lowcon,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(3,:,:))',...
    squeeze(resp_oriXcon(3,:,:))',double(orientations),1);
ori_fit_lowcon = squeeze(ori_fit_lowcon(1:180,1,:));

[~,pref_ori] = max(squeeze(resp_oriXcon(end,:,:)));

nexamplecells = 8;
if nexamplecells > nc
    nexamplecells = nc;
end
ind = randsample(find(sum(squeeze(respCell_oriXcon(end,:,:)),1)>0),nexamplecells);

figure
for i = 1:nexamplecells
    y = resp_oriXcon(end,:,ind(i));
    yerr = respErr_oriXcon(end,:,ind(i));
    f = ori_fit(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'ko')
    plot(0:179,f,'k-')
    
    y = resp_oriXcon(3,:,ind(i));
    yerr = respErr_oriXcon(3,:,ind(i));
    f = ori_fit_lowcon(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'mo')
    plot(0:179,f,'m-')
    
    figXAxis([],'Orienation (deg)',[0 180],orientations,orientations)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell #%s',num2str(ind(i))))
end
subplot 339
fractionResponsive = sum(squeeze(sum(respCell_oriXcon,2))>0,2)./nc;
plot(contrasts,fractionResponsive,'.','MarkerSize',20)
figXAxis([],'Contrast',[0 1],contrasts,contrasts)
set(gca, 'XScale', 'log')
figYAxis([],'Fraction Responsive Cells',[])
figAxForm
title(sprintf('DART %s HRs', num2str(expt(day1_id).multiday_timesincedrug_hours)))

print(fullfile(fnout,'oriTuning_tagged_exCells_day1'),'-dpdf','-fillpage')

%not tagged
[dff_eaTrial,avgResp_eaTrial] = getEaTrialResp_visstimret(data_notag{1},mw_day1,frameRateHz);

[nt,nc] = size(avgResp_eaTrial);

resp_oriXcon = nan(ncon,nori,nc);
respErr_oriXcon = nan(ncon,nori,nc);
respCell_oriXcon = nan(ncon,nori,nc);
for icon = 1:ncon
    for iori = 1:nori
        ind = conID == contrasts(icon) & oriID == orientations(iori);
        resp_oriXcon(icon,iori,:) = mean(avgResp_eaTrial(ind,:),1);
        respErr_oriXcon(icon,iori,:) = ste(avgResp_eaTrial(ind,:),1);
%         if icon == ncon
            b = cell2mat(cellfun(@(x) mean(x(basewin,:)),dff_eaTrial(ind),'unif',0)');
            r = cell2mat(cellfun(@(x) mean(x(respwin,:)),dff_eaTrial(ind),'unif',0)');
            respCell_oriXcon(icon,iori,:) = ttest(r,b,'dim',1,'tail','right','alpha',0.01./nori);
%         end
    end
end


[ori_fit,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(end,:,:))',...
    squeeze(resp_oriXcon(end,:,:))',double(orientations),1);
ori_fit = squeeze(ori_fit(1:180,1,:));
[ori_fit_lowcon,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(3,:,:))',...
    squeeze(resp_oriXcon(3,:,:))',double(orientations),1);
ori_fit_lowcon = squeeze(ori_fit_lowcon(1:180,1,:));

[~,pref_ori] = max(squeeze(resp_oriXcon(end,:,:)));

nexamplecells = 8;
if nexamplecells > nc
    nexamplecells = nc;
end
ind = randsample(find(sum(squeeze(respCell_oriXcon(end,:,:)),1)>0),nexamplecells);

figure
for i = 1:nexamplecells
    y = resp_oriXcon(end,:,ind(i));
    yerr = respErr_oriXcon(end,:,ind(i));
    f = ori_fit(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'ko')
    plot(0:179,f,'k-')
    
    y = resp_oriXcon(3,:,ind(i));
    yerr = respErr_oriXcon(3,:,ind(i));
    f = ori_fit_lowcon(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'mo')
    plot(0:179,f,'m-')
    
    figXAxis([],'Orienation (deg)',[0 180],orientations,orientations)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell #%s',num2str(ind(i))))
end
subplot 339
fractionResponsive = sum(squeeze(sum(respCell_oriXcon,2))>0,2)./nc;
plot(contrasts,fractionResponsive,'.','MarkerSize',20)
figXAxis([],'Contrast',[0 1],contrasts,contrasts)
set(gca, 'XScale', 'log')
figYAxis([],'Fraction Responsive Cells',[])
figAxForm
title(sprintf('DART %s HRs', num2str(expt(day1_id).multiday_timesincedrug_hours)))

print(fullfile(fnout,'oriTuning_nottag_exCells_day1'),'-dpdf','-fillpage')
%% day 2
%tagged
[dff_eaTrial,avgResp_eaTrial] = getEaTrialResp_visstimret({data_tag{2}},{mw_day2},frameRateHz);

conID = cell2mat(cellfun(@(x) celleqel2mat_padded(x.tGratingContrast),{mw_day2},'unif',0));
dirID = cell2mat(cellfun(@(x) celleqel2mat_padded(x.tGratingDirectionDeg),{mw_day2},'unif',0));
oriID = nan(size(dirID));
oriID(dirID<180) = dirID(dirID<180);
oriID(dirID>=180) = dirID(dirID>=180)-180;
contrasts = unique(conID);
orientations = unique(oriID);
ncon = length(contrasts);
nori = length(orientations);
[nt,nc] = size(avgResp_eaTrial);

resp_oriXcon = nan(ncon,nori,nc);
respErr_oriXcon = nan(ncon,nori,nc);
respCell_oriXcon = nan(ncon,nori,nc);
for icon = 1:ncon
    for iori = 1:nori
        ind = conID == contrasts(icon) & oriID == orientations(iori);
        resp_oriXcon(icon,iori,:) = mean(avgResp_eaTrial(ind,:),1);
        respErr_oriXcon(icon,iori,:) = ste(avgResp_eaTrial(ind,:),1);
%         if icon == ncon
            b = cell2mat(cellfun(@(x) mean(x(basewin,:)),dff_eaTrial(ind),'unif',0)');
            r = cell2mat(cellfun(@(x) mean(x(respwin,:)),dff_eaTrial(ind),'unif',0)');
            respCell_oriXcon(icon,iori,:) = ttest(r,b,'dim',1,'tail','right','alpha',0.01./nori);
%         end
    end
end


[ori_fit,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(end,:,:))',...
    squeeze(resp_oriXcon(end,:,:))',double(orientations),1);
ori_fit = squeeze(ori_fit(1:180,1,:));
[ori_fit_lowcon,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(3,:,:))',...
    squeeze(resp_oriXcon(3,:,:))',double(orientations),1);
ori_fit_lowcon = squeeze(ori_fit_lowcon(1:180,1,:));

[~,pref_ori] = max(squeeze(resp_oriXcon(end,:,:)));

nexamplecells = 8;
if nexamplecells > nc
    nexamplecells = nc;
end
ind = find(sum(squeeze(respCell_oriXcon(end,:,:)),1)>0);
if nexamplecells > length(ind)
    nexamplecells = length(ind);
end
ind = randsample(ind,nexamplecells);

figure
for i = 1:nexamplecells
    y = resp_oriXcon(end,:,ind(i));
    yerr = respErr_oriXcon(end,:,ind(i));
    f = ori_fit(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'ko')
    plot(0:179,f,'k-')
    
    y = resp_oriXcon(3,:,ind(i));
    yerr = respErr_oriXcon(3,:,ind(i));
    f = ori_fit_lowcon(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'mo')
    plot(0:179,f,'m-')
    
    figXAxis([],'Orienation (deg)',[0 180],orientations,orientations)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell #%s',num2str(ind(i))))
end
subplot 339
fractionResponsive = sum(squeeze(sum(respCell_oriXcon,2))>0,2)./nc;
plot(contrasts,fractionResponsive,'.','MarkerSize',20)
figXAxis([],'Contrast',[0 1],contrasts,contrasts)
set(gca, 'XScale', 'log')
figYAxis([],'Fraction Responsive Cells',[])
figAxForm
title(sprintf('DART %s HRs', num2str(expt(day2_id).multiday_timesincedrug_hours)))

print(fullfile(fnout,'oriTuning_tagged_exCells_day2'),'-dpdf','-fillpage')

%no tag
[dff_eaTrial,avgResp_eaTrial] = getEaTrialResp_visstimret({data_notag{2}},{mw_day2},frameRateHz);

[nt,nc] = size(avgResp_eaTrial);

resp_oriXcon = nan(ncon,nori,nc);
respErr_oriXcon = nan(ncon,nori,nc);
respCell_oriXcon = nan(ncon,nori,nc);
for icon = 1:ncon
    for iori = 1:nori
        ind = conID == contrasts(icon) & oriID == orientations(iori);
        resp_oriXcon(icon,iori,:) = mean(avgResp_eaTrial(ind,:),1);
        respErr_oriXcon(icon,iori,:) = ste(avgResp_eaTrial(ind,:),1);
%         if icon == ncon
            b = cell2mat(cellfun(@(x) mean(x(basewin,:)),dff_eaTrial(ind),'unif',0)');
            r = cell2mat(cellfun(@(x) mean(x(respwin,:)),dff_eaTrial(ind),'unif',0)');
            respCell_oriXcon(icon,iori,:) = ttest(r,b,'dim',1,'tail','right','alpha',0.01./nori);
%         end
    end
end


[ori_fit,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(end,:,:))',...
    squeeze(resp_oriXcon(end,:,:))',double(orientations),1);
ori_fit = squeeze(ori_fit(1:180,1,:));
[ori_fit_lowcon,~,~,R2] = vonmisesReliableFit(squeeze(resp_oriXcon(3,:,:))',...
    squeeze(resp_oriXcon(3,:,:))',double(orientations),1);
ori_fit_lowcon = squeeze(ori_fit_lowcon(1:180,1,:));

[~,pref_ori] = max(squeeze(resp_oriXcon(end,:,:)));

nexamplecells = 8;
if nexamplecells > nc
    nexamplecells = nc;
end
ind = find(sum(squeeze(respCell_oriXcon(end,:,:)),1)>0);
if nexamplecells > length(ind)
    nexamplecells = length(ind);
end
ind = randsample(ind,nexamplecells);

figure
for i = 1:nexamplecells
    y = resp_oriXcon(end,:,ind(i));
    yerr = respErr_oriXcon(end,:,ind(i));
    f = ori_fit(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'ko')
    plot(0:179,f,'k-')
    
    y = resp_oriXcon(3,:,ind(i));
    yerr = respErr_oriXcon(3,:,ind(i));
    f = ori_fit_lowcon(:,ind(i));
    subplot(3,3,i);hold on
    errorbar(orientations,y,yerr,'mo')
    plot(0:179,f,'m-')
    
    figXAxis([],'Orienation (deg)',[0 180],orientations,orientations)
    figYAxis([],'dF/F',[])
    figAxForm
    title(sprintf('Cell #%s',num2str(ind(i))))
end
subplot 339
fractionResponsive = sum(squeeze(sum(respCell_oriXcon,2))>0,2)./nc;
plot(contrasts,fractionResponsive,'.','MarkerSize',20)
figXAxis([],'Contrast',[0 1],contrasts,contrasts)
set(gca, 'XScale', 'log')
figYAxis([],'Fraction Responsive Cells',[])
figAxForm
title(sprintf('DART %s HRs', num2str(expt(day2_id).multiday_timesincedrug_hours)))

print(fullfile(fnout,'oriTuning_notag_exCells_day2'),'-dpdf','-fillpage')