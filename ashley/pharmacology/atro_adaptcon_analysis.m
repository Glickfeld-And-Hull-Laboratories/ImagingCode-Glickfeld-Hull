clear all
clear global
close all
ds = 'atro_adaptcon_V1'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doPreviousReg = false;
doGreenOnly = false;
%%
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');

%% load time-courses 
load(fullfile(fnout,'timecourses_cells'))

%% parameters

frame_rate = expt(slct_expt).frame_rate;
nBaseFr = frame_rate;
nRespFr = frame_rate;
winlengthFr = round(0.2*frame_rate);
basewin = (nBaseFr-winlengthFr):nBaseFr;
respwin = nBaseFr+4:(nBaseFr+4+winlengthFr);

tt_ms = ((1:(nBaseFr+nRespFr))-nBaseFr)./frame_rate.*1000;

nc = [size(tc_tagneg_adapt_subnp{1},2),size(tc_tagpos_adapt_subnp{1},2)];

%% labels
tagName = {'no tag',expt(slct_expt).redChannelLabel};
drugs = expt(slct_expt).adapt_condition;
%% spontaneous activity

%% contrast responses
%load mworks info
nruns = length(expt(slct_expt).adapt_runs);
data_mw = cell(1,nruns);
for irun = 1:nruns    
    mw_temp = loadMworksFile(mouse,expDate,expt(slct_expt).adapt_time{irun});
    data_mw{irun} = mw_temp;
end

% grab trial info
tc_stim1 = cell(1,nruns);
tc_stim2 = cell(1,nruns);
tCons = cell(1,nruns);
tISI_ms = cell(1,nruns);
stimOn_ms = double(data_mw{1}.stimOneGratingOnTimeMs);
for irun = 1:nruns
    
    tCons{irun} = celleqel2mat_padded(data_mw{irun}.tStimOneGratingContrast);
    cStim1 = celleqel2mat_padded(data_mw{irun}.cStimOneOn);
    cStim2 = celleqel2mat_padded(data_mw{irun}.cStimTwoOn);
    tISI_ms{irun} = (round((cStim2-cStim1)./frame_rate,2,'decimals').*1000)-stimOn_ms;

    nt = length(cStim1);
    d1 = nan(nt,nBaseFr+nRespFr,nc(1));
    d2 = nan(nt,nBaseFr+nRespFr,nc(1));
    for itrial = 1:nt
        ind = (cStim1(itrial)-nBaseFr):(cStim1(itrial)+nRespFr-1);
        d1(itrial,:,:) = tc_tagneg_adapt_subnp{irun}(ind,:);
        
        ind = (cStim2(itrial)-nBaseFr):(cStim2(itrial)+nRespFr-1);
        d2(itrial,:,:) = tc_tagneg_adapt_subnp{irun}(ind,:);
    end
    tc_stim1{irun} = d1;
    tc_stim2{irun} = d2;
end

% get dff ea trial
dff_stim1 = cellfun(@(x) (x-mean(x(:,basewin,:),2))./mean(x(:,basewin,:),2),...
    tc_stim1,'unif',0);
dff_stim2 = cellfun(@(x) (x-mean(x(:,basewin,:),2))./mean(x(:,basewin,:),2),...
    tc_stim2,'unif',0);

% get responses ea trial
resp_stim1 = cellfun(@(x) squeeze(mean(x(:,respwin,:),2)),dff_stim1,'unif',0);
resp_stim2 = cellfun(@(x) squeeze(mean(x(:,respwin,:),2)),dff_stim2,'unif',0);

% get adaptation ea trial 
adapt = cellfun(@(x,y) y./x,resp_stim1,resp_stim2,'unif',0);

% get responses ea contrast
cons = unique(cell2mat(tCons));
isis = unique(cell2mat(tISI_ms));
tMaxAdaptISI = cellfun(@(x) x == isis(1)|x==isis(2),tISI_ms,'unif',0);
conResp_stim1 = cell(1,nruns);
conResp_stim2 = cell(1,nruns);
conResp_err_stim1 = cell(1,nruns);
conAdapt = cell(1,nruns);
conRespCells = cell(1,nruns);
for icon = 1:length(cons)
    ind = cellfun(@(x) x == cons(icon),tCons,'unif',0);
    conResp_stim1 = cellfun(@(a,b) cat(1,a,b),conResp_stim1,...
        cellfun(@(x,y) mean(x(y,:),1),resp_stim1,ind,'unif',0),'unif',0);
    conResp_err_stim1 = cellfun(@(a,b) cat(1,a,b),conResp_stim1,...
        cellfun(@(x,y) ste(x(y,:),1),resp_stim1,ind,'unif',0),'unif',0);
    conResp_stim2 = cellfun(@(a,b) cat(1,a,b),conResp_stim2,...
        cellfun(@(x,y) mean(x(y,:),1),resp_stim2,ind,'unif',0),'unif',0);
    
%     conAdapt = cellfun(@(a,b) cat(1,a,b),conAdapt,...
%         cellfun(@(x,y) mean(x(y,:),1),adapt,ind,'unif',0),'unif',0);
    conAdapt = cellfun(@(a,b) cat(1,a,b),conAdapt,...
        cellfun(@(x,y,z,xx) mean(y(z & xx,:),1)./mean(x(z & xx,:),1),...
        resp_stim1,resp_stim2,ind,tMaxAdaptISI,'unif',0),'unif',0);
    
    conRespCells = cellfun(@(a,b) cat(1,a,b),conRespCells,...
        cellfun(@(x,y) ttest(squeeze(mean(x(y,respwin,:),2)),...
        squeeze(mean(x(y,basewin,:),2)),'dim',1,'tail','right'),dff_stim1,ind,...
        'unif',0),'unif',0);
end

%% more params
conColors = brewermap(length(cons)+3,'Greys'); conColors = conColors(4:end,:);

%% look at timecourses of responsive cells across conditions
resp_lim = [-0.02, 0.2];
respCells = sum(conRespCells{1},1)>0 | ...
    ttest(squeeze(mean(dff_stim1{1}(:,respwin,:),2)),...
    squeeze(mean(dff_stim1{1}(:,basewin,:),2)),'dim',1,'tail','right');

for icell = 1:nc
    if respCells(icell) ~= 1
        continue
    end
    figure
    suptitle(sprintf('Atropine, cell # %s',num2str(icell)))
    for idrug = 1:2
        subplot(3,2,idrug)
        hold on
        y = nan(1,length(cons));
        yerr = nan(1,length(cons));
        for icon = 1:length(cons)
            ind = tCons{idrug} == cons(icon) & tMaxAdaptISI{idrug}==1;
            tc = squeeze(mean(dff_stim1{idrug}(ind,:,icell)));
            h = plot(tt_ms,tc,'-','LineWidth',1);
            h.Color = conColors(icon,:);
            if sum(ind) > 0
                y(icon) = mean(resp_stim1{idrug}(ind,icell),1);
                yerr(icon) = ste(resp_stim1{idrug}(ind,icell),1);
            end
        end
        figXAxis([],'time from stim (ms)',[min(tt_ms),max(tt_ms)])
        figYAxis([],'dF/F',resp_lim)
        figAxForm
        title(sprintf('stim 1 + %s',drugs{idrug}))

        subplot(3,2,idrug+2)
        hold on
        for icon = 1:length(cons)
            ind = tCons{idrug} == cons(icon) & tMaxAdaptISI{idrug}==1;
            tc = squeeze(mean(dff_stim2{idrug}(ind,:,icell)));
            h = plot(tt_ms,tc,'-','LineWidth',1);
            h.Color = conColors(icon,:);
        end
        figXAxis([],'time from stim (ms)',[min(tt_ms),max(tt_ms)])
        figYAxis([],'dF/F',resp_lim)
        figAxForm
        title(sprintf('stim 2 + %s',drugs{idrug}))
        
        subplot(3,2,5)
        hold on
        x = cons(~isnan(y));
        if idrug == 1
            errorbar(x,y(~isnan(y)),yerr(~isnan(y)),'k.-','MarkerSize',20)
        else
            errorbar(x,y(~isnan(y)),yerr(~isnan(y)),'r.-','MarkerSize',20)
        end
        figXAxis([],'contrast',[0 1.1],cons,cons)
        figYAxis([],'dF/F',resp_lim)
        figAxForm
        if icell == 13
            print(fullfile(fnout,'exCellTC_notag'),'-dpdf','-fillpage')
        end
    end    
end

%% plot responsees across conditions
resp_lim = [-0.05 .25];
resp_lim_low = [-0.05 .1];
respCells = sum(conRespCells{1},1)>0 | ...
    ttest(squeeze(mean(dff_stim1{1}(:,respwin,:),2)),...
    squeeze(mean(dff_stim1{1}(:,basewin,:),2)),'dim',1,'tail','right');

% normResp_stim1 = cellfun(@(x) x(2:4,:)./max(conResp_stim1{1}(2:4,:),[],1),conResp_stim1,'unif',0);

figure
subplot 131
hold on
ind = cellfun(@(x) x==0.2,tCons,'unif',0);
r = cellfun(@(x,y) mean(x(y,respCells),1), resp_stim1,ind,'unif',0);
plot(r{1},r{2},'.','MarkerSize',20)
figXAxis([],'saline (dF/F)',resp_lim_low)
figYAxis([],'drug (dF/F)',resp_lim_low)
figAxForm
plot(resp_lim,resp_lim,'k--')
title('con=0.2 trials')

subplot 132
hold on
ind = cellfun(@(x) x==0.8,tCons,'unif',0);
r = cellfun(@(x,y) mean(x(y,respCells),1), resp_stim1,ind,'unif',0);
plot(r{1},r{2},'.','MarkerSize',20)
figXAxis([],'saline (dF/F)',resp_lim)
figYAxis([],'drug (dF/F)',resp_lim)
figAxForm
plot(resp_lim,resp_lim,'k--')
title('con=0.8 trials')

subplot 133
hold on
y = mean(conResp_stim1{1}(:,respCells),2);
yerr = ste(conResp_stim1{1}(:,respCells),2);
x = cons(~isnan(y));
errorbar(x,y(~isnan(y)),yerr(~isnan(y)),'k.-','MarkerSize',20)
y = mean(conResp_stim1{2}(:,respCells),2);
yerr = ste(conResp_stim1{2}(:,respCells),2);
x = cons(~isnan(y));
errorbar(x,y(~isnan(y)),yerr(~isnan(y)),'r.-','MarkerSize',20)
figXAxis([],'contrast',[0.05 1.1],cons,cons)
figYAxis([],'dF/F',resp_lim_low)
figAxForm
ax = gca; ax.XScale = 'log';

print(fullfile(fnout,'contrastResp_notag'),'-dpdf','-fillpage')

%% plot adaptation
resp_lim = [-0.05 0.15];
adapt_lim = [-0.5 2];
con_edges = -2:0.2:2;


adapt_fig = figure;
suptitle(sprintf('%s expt',expt(slct_expt).drug{1}))
for idrug = 1:2
    if idrug == 1
        iplot = 0;
    else
        iplot = iplot+2;
    end
    
    subplot(3,2,1+iplot)
    hold on
    for icon = 1:length(cons)
        x = conResp_stim1{idrug}(icon,respCells);
        y = conResp_stim2{idrug}(icon,respCells);
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = conColors(icon,:);
    end
    figXAxis([],'stim 1 dF/F',resp_lim)
    figYAxis([],'stim 2 dF/F',resp_lim)
    figAxForm
    plot(resp_lim,resp_lim,'k--')
    title(drugs{idrug})
    
    subplot(3,2,2+iplot)
    hold on
    for icon = 1:length(cons)
        ind = conRespCells{idrug}(icon,:)==1 & ...
            conResp_stim1{idrug}(icon,:)>0.015;
        a = conAdapt{idrug}(icon,ind);
        if ~isempty(a) 
            h = plot(ones(1,sum(ind)).*icon,a,'.','MarkerSize',10);
            h.Color = conColors(icon,:);
            h = plot(icon,mean(a),'.','MarkerSize',20);
            h.Color = conColors(icon,:);
        end
    end
    hline(1,'k--')
    figXAxis([],'contrast',[0 length(cons)+1],1:length(cons),cons)
    figYAxis([],'stim2/stim1',adapt_lim)
    figAxForm
    title(drugs{idrug})
end
subplot(3,2,5)
hold on
for icon = 1:length(cons)
    ind = conRespCells{1}(icon,:)==1 & conRespCells{2}(icon,:)==1 & conResp_stim1{idrug}(icon,:)>0.015;
    x = conAdapt{1}(icon,ind);
    y = conAdapt{2}(icon,ind);
    if isempty(x) || isempty(y)
        continue
    end
    h = plot(x,y,'.','MarkerSize',20);
    h.Color = conColors(icon,:);
end
figXAxis([],'saline - stim2/stim1',adapt_lim)
figYAxis([],'drug - stim2/stim1',adapt_lim)
figAxForm
plot(adapt_lim,adapt_lim,'k--')
print(fullfile(fnout,'adaptation_notag'),'-dpdf','-fillpage')

%% noise correlation
corr_lim = [-1 1];
respCells = sum(conRespCells{1},1)>0 | ...
    ttest(squeeze(mean(dff_stim1{1}(:,respwin,:),2)),...
    squeeze(mean(dff_stim1{1}(:,basewin,:),2)),'dim',1,'tail','right');
ind = cellfun(@(x) x==0.8,tCons,'unif',0);
r = cellfun(@(x,y) x(y,respCells),resp_stim1,ind,'unif',0);
rsc = cellfun(@corrcoef,r,'unif',0);   
[~,sortInd] = sort(mean(rsc{1}));

figure; colormap(brewermap([],'*RdBu'))
for idrug = 1:nruns
    subplot(1,nruns+1,idrug)
    imagesc(rsc{idrug}(sortInd,sortInd))
    figAxForm
    clim(corr_lim)
    colorbar
    title(drugs{idrug})
end
subplot(1,nruns+1,nruns+1)
imagesc(rsc{2}(sortInd,sortInd)-rsc{1}(sortInd,sortInd))
figAxForm
clim(corr_lim)
colorbar
title('subtraction')

print(fullfile(fnout,'noiseCorr_notag'),'-dpdf','-fillpage')

% %% running
% runningThreshold_cps = 5;
% ws_cps = cellfun(@(x) wheelSpeedCalc(x,32,'red'),data_mw,'unif',0);
% running_tc = cell(1,nruns);
% for irun = 1:nruns
%     
%     cStim1 = celleqel2mat_padded(data_mw{irun}.cStimOneOn);
%     
%     nt = length(cStim1);
%     d1 = nan(nt,nBaseFr+nRespFr);
%     for itrial = 1:nt
%         ind = (cStim1(itrial)-nBaseFr):(cStim1(itrial)+nRespFr-1);
%         d1(itrial,:) = ws_cps{irun}(ind);
%     end
%     running_tc{irun} = d1;
% end
% 
% running_prob = cellfun(@(x) ...
%     sum(x(:,respwin)>runningThreshold_cps,2)./length(respwin),...
%     running_tc,'unif',0);
% 
% figure; 
% subplot 121