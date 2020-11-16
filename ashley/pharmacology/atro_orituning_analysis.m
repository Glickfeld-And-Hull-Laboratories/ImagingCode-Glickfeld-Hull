clear all
clear global
close all
ds = 'nbqx_oriTuning_V1'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doPreviousReg = false;
doGreenOnly = false;
dsFactor = 10;
%%
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');

%% load time-courses and ori tuning
load(fullfile(fnout,'timecourses_cells'))

%% spontaneous

Fwin_s = 120;
Fwin_fr = Fwin_s.*expt(slct_expt).frame_rate./dsFactor;
nframes = cellfun(@(x) size(x,1),tc_tagneg_nostim_subnp);

% calc n events as 2std gt mean, calc df/f
events_nostim = cell(1,2);
dff_nostim = cell(1,2);
for idrug = 1:2
    d = tc_tagneg_nostim_subnp{idrug};
    nfr = size(d,1);
    halfwin = floor(Fwin_fr/2);
    e = nan(nfr,size(d,2));
    dff = nan(nfr,size(d,2));
    for iframe = 1:nfr
        if iframe < halfwin
            ind = 1:Fwin_fr;
        elseif iframe > (nfr-halfwin)
            ind = (nfr-Fwin_fr+1):nfr;
        else
            ind = (iframe-halfwin+1):(iframe+halfwin);
        end
        F = d(ind,:);
        
        F_mean = mean(F,1);
        F_gt_2std = 2*std(F,[],1)+F_mean;
        e(iframe,:) = d(iframe,:) > F_gt_2std;    
        
        dff(iframe,:) = (d(iframe,:)-mean(F,1))./mean(F,1);
    end
    events_nostim{idrug} = e;
    dff_nostim{idrug} = dff;
end


resp_lim = [-0.1 0.25];
exCells = [11,34];
figure
% expt time course
for idrug = 1:2
    for icell = 1:2
        subplot(4,2,((icell-1)*2)+idrug)
        tt = (1:nframes(idrug))./(expt(slct_expt).frame_rate./dsFactor);
        y = dff_nostim{idrug}(:,exCells(icell));
        plot(tt,y)
        figXAxis([],'time (s)',[1 max(tt)])
        figYAxis([],'df/f',resp_lim)
        figAxForm([],0)
    end
    
    subplot(4,2,idrug+4)
    y = mean(events_nostim{idrug},2);
    plot(tt,y)
    figXAxis([],'time (s)',[1 max(tt)])
    figYAxis([],'avg events across cells',resp_lim)
    figAxForm([],0)
end
subplot 427
e1 = mean(mean(events_nostim{1}));
e1_err =  ste(mean(events_nostim{1},1),2);
win = round(nframes(2)/3);
e2_bins = {1:win;(win+1):(2*win);((2*win)+1):nframes(2)};
e2 = cellfun(@(x) mean(mean(events_nostim{2}(x,:))),e2_bins);
e2_err = cellfun(@(x) ste(mean(events_nostim{2}(x,:),1),2),e2_bins);
errorbar(1:4,[e1;e2],[e1_err;e2_err],'.','MarkerSize',20)
figXAxis([],'expt',[0 5],1:4,{'pre-drug';'drug t1';'drug t2';'drug t3'})
figYAxis([],'avg events/frame across cells',[])
figAxForm

%% tuning responses
load(fullfile(fnout,'oriTuning'))
nruns = length(expt(slct_expt).dirtuning_runs);
data_mw = cell(1,nruns);
for irun = 1:nruns    
    mw_temp = loadMworksFile(mouse,expDate,expt(slct_expt).dirtuning_time{irun});
    data_mw{irun} = mw_temp;
end
directions = unique(cell2mat_padded(data_mw{1}.tGratingDirectionDeg));
orientations = directions(directions <180);

cellInd = isResp_tagneg{1} & fitReliability_tagneg{1}' < 30;

% change in direction tuning curve

exCell = 1;
conditions = expt(slct_expt).dirtuning_condition;
colors = cat(1,[0,0,0],[0.9922,0.7059,0.3843]);

figure
[nrows,ncols] = optimizeSubplotDim(sum(cellInd));
for icell = 1:length(cellInd)
    if icell == 1
        iplot = 0;
    end
    if cellInd(icell) == 1
        iplot = iplot+1;
        subplot(nrows,ncols,iplot)
        for icond = 1:length(conditions)
            r = avgResponseEaOri_tagneg{icond}(icell,:);
            r_err = semResponseEaOri_tagneg{icond}(icell,:);
            f = vonMisesFitAllCells_tagneg{icond}(:,1,icell);
            
            h = errorbar(orientations,r,r_err,'.','MarkerSize',10);
            h.Color = colors(icond,:);
            hold on
            h = plot(0:180,f,'-');
            h.Color = colors(icond,:);
            figXAxis([],'orientation',[-10 190],0:45:180,0:45:180)
            figYAxis([],'dF/F',[])
            figAxForm
            if icell == exCell
                title(['*' num2str(icell)])
            else
            title(num2str(icell))
            end
        end
    end
end
print(fullfile(fnout,'oriTuning_respCells'),'-dpdf','-fillpage')

% change in tc of response to preferred stim
nstim = length(orientations);
[~,prefstim_ind] = max(avgResponseEaOri_tagneg{1}(exCell,:));

on = data_mw{1}.nScansOn./dsFactor;
off = data_mw{1}.nScansOff./dsFactor;
tt_s = double((-off+1):on)./expt(slct_expt).frame_rate;

figure
[nrows,ncols] = optimizeSubplotDim(length(orientations));
for iori = 1:nstim
    subplot(nrows,ncols,iori)
    hold on
    for icond = 1:length(conditions)
        y = tuningTC_tagneg{icond}(:,exCell,iori);
        h = plot(tt_s,y,'-');
        h.Color = colors(icond,:);
    end
    figXAxis([],'time from stim (s)',[min(tt_s) max(tt_s)])
    figYAxis([],'dF/F',[])
    figAxForm
    if iori == prefstim_ind
        title(sprintf('%s deg *pref',num2str(orientations(iori))))
    else
        title(sprintf('%s deg',num2str(orientations(iori))))
    end
end
print(fullfile(fnout,'respTC_exCell'),'-dpdf','-fillpage')

%% FOV change in response
% load(fullfile(fnout,selection_