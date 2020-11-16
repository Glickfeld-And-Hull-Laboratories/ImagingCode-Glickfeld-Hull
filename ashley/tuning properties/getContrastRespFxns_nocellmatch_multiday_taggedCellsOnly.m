%% which dataset and parameters
clear all
clear global
close all
ds = 'DART_V1_PV_multiday'; %dataset info
dataStructLabels = 'contrastxsize';

rc = behavConstsAV; %directories
eval(ds

%%
size2use = 20;

%% which days to use?
day1_id = 1;
day2_id = 2;
day3_id = 3;

%% mouse and saving info
mouse = expt(day1_id).mouse;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
fnout = fullfile(fn,'multiday');
mkdir(fnout)

%% load time-courses
% day 1
expDate_day1 = expt(day1_id).date;
expTime_day1 = eval(['expt(day1_id).' dataStructLabels '_time']);
load(fullfile(fn,expDate_day1,eval(['expt(day1_id).' dataStructLabels '_runs{1}']),'timecourses_cells'))
day1_tc = tc_subnp_red;

% day 2
expDate_day2 = expt(day2_id).date;
expTime_day2 = eval(['expt(day2_id).' dataStructLabels '_time']);
load(fullfile(fn,expDate_day2,eval(['expt(day2_id).' dataStructLabels '_runs{1}']),'timecourses_cells'))
day2_tc = tc_subnp_red;

% day 3
expDate_day3 = expt(day3_id).date;
expTime_day3 = eval(['expt(day3_id).' dataStructLabels '_time']);
load(fullfile(fn,expDate_day3,eval(['expt(day3_id).' dataStructLabels '_runs{1}']),'timecourses_cells'))
day3_tc = tc_subnp_red;

data = {day1_tc,day2_tc,day3_tc};
%% load mworks data
mw_day1 = loadMworksFile(mouse,expDate_day1,expTime_day1{1});
mw_day2 = loadMworksFile(mouse,expDate_day2,expTime_day2{1});
mw_day3 = loadMworksFile(mouse,expDate_day3,expTime_day3{1});

mw = {mw_day1,mw_day2,mw_day3};
%% get contrast response tuning
[conXsize_resp,contrastXsize_respErr,contrasts,sizes,contrastXsize_respData] = cellfun(@(x,y) getContrastBySizeResp(...
    x,y,expt(day1_id).frame_rate),data,mw,'unif',0);
%% fit contrast response
[fit,R2,c50] = cellfun(@(x,y) fitContrastResp(x,y,expt(day1_id).frame_rate),...
    contrastResp,contrasts,'unif',0);

%% normalized data
fit_norm = cellfun(@(x) x./max(x), fit,'unif',0);
c50_norm = cellfun(@(x) getC50(x),fit_norm,'unif',0);

%% find contrast responsive cells
fitThreshold = 0.6;
[~,size2use_ind] = min(abs(sizes{1}-size2use));
[isResponsive_any, isResponsive_eaStim] = cellfun(@(x) findRespCells_anyContrast(...
    x(size2use_ind,:),expt(day1_id).frame_rate),contrastXsize_respData,'unif',0);
% maxConIsMaxCon = cellfun(@(x) max(x)==x(end,:)| max(x)==x(end-1,:),...
%     contrastResp,'unif',0);
% contrastModulatedCells = cellfun(@(x,y,z) x&y'&z>fitThreshold,...
%     maxConIsMaxCon,isResponsive_any,R2,'unif',0);
isResponsive_maxCon = cellfun(@(x) x(:,end),isResponsive_eaStim,'unif',0);
%% plot tcs of each cell to max contrast, mid size
tt = ((-expt(day1_id).frame_rate+1):expt(day1_id).frame_rate)./expt(day1_id).frame_rate.*1000;
nc = cellfun(@(x) size(x{1,1},2),contrastXsize_respData);

respcell_fig = figure;
for iday = 1:3
    [nrow,ncol] = optimizeSubplotDim(nc(iday));
    dt = eval(['expt(day' num2str(iday),'_id).date']);
    figure
    suptitle(dt)
    for icell = 1:nc(iday)
        tc = mean(contrastXsize_respData{iday}{size2use_ind,end}(:,icell,:),3);
        tc_err = ste(contrastXsize_respData{iday}{size2use_ind,end}(:,icell,:),3);
        subplot(nrow,ncol,icell)
        shadedErrorBar_chooseColor(tt,tc,tc_err,[0 0 0]);
        figXAxis([],'time from stim (ms)',[min(tt) max(tt)])
        figYAxis([],'dF/F',[])
        figAxForm
        if isResponsive_maxCon{iday}(icell) == 1
            title('**')
        end
    end
    print(fullfile(fnout,['eaTaggedCell_' dt]),'-dpdf','-fillpage')
end

%% spontaneous activity
pct_events = cell(1,3);
for iday = 1:3
    ind = find(isResponsive_maxCon{iday}== 1);
    dt = eval(['expt(day' num2str(iday),'_id).date']);
    figure
    suptitle([dt ', resp cells'])
    d = eval(['day' num2str(iday) '_tc']);
    dff = (d-mean(d,1))./mean(d,1);
    tt = (1:size(d,1))./expt(day1_id).frame_rate./60;
    for icell = 1:length(ind)
        tc = dff(:,ind(icell));
        subplot(length(ind),1,icell)
        plot(tt,tc,'k-','LineWidth',1);
        figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
        figYAxis([],'dF/F',[-1 2])
        figAxForm([],0)
    end
    print(fullfile(fnout,['spontaneous_dff_' dt]),'-dpdf','-fillpage')
    dff_3sd = (std(dff) + mean(dff))*3;
    dff_test = dff > dff_3sd;
    pct_events{iday} = sum(dff_test,1)./size(dff,1);
end

figure; hold on
for iday = 1:3
    d = pct_events{iday};
    plot(ones(nc(iday),1)*iday,d,'k.','MarkerSize',10)
    errorbar(iday,mean(d),ste(d,2),'b.','MarkerSize',20)
end
figXAxis([],'',[0 4],1:3,{'pre-DART','2hrs','24hrs'})
figYAxis([],'spontaneous events',[])
figAxForm
print(fullfile(fnout,'spontaneousEvents'),'-dpdf')