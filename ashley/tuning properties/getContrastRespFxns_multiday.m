%% which dataset and parameters
clear all
clear global
close all
ds = 'multiday_test'; %dataset info
dataStructLabels = 'contrast';

rc = behavConstsAV; %directories
eval(ds)
%%
day1_id = 1;
day2_id = 2;
day3_id = nan;

%%
mouse = expt(day1_id).mouse;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
fnout = fullfile(fn,'multiday');
mkdir(fnout)

%% load time-courses
load(fullfile(fnout,'timecourses'))

goodCells = cell2mat({cellTCs.pass}) == 1;
day1_tc = cell2mat({cellTCs.tc_day1});
day2_tc = cell2mat({cellTCs.tc_day2});
data = {day1_tc,day2_tc};
%% load mworks data
expDate_day1 = expt(day1_id).date;
expTime_day1 = eval(['expt(day1_id).' dataStructLabels '_time']);
mw_day1 = loadMworksFile(mouse,expDate_day1,expTime_day1{1});

expDate_day2 = expt(day2_id).date;
expTime_day2 = eval(['expt(day2_id).' dataStructLabels '_time']);
mw_day2 = loadMworksFile(mouse,expDate_day2,expTime_day2{1});

mw = {mw_day1,mw_day2};
%% get contrast response tuning
[contrastResp,contrastRespErr,~,~,contrasts,contrastRespData] = cellfun(@(x,y) getContrastResp(...
    x,y,expt(day1_id).frame_rate),data,mw,'unif',0);
%% fit contrast response
[fit,R2,c50] = cellfun(@(x,y) fitContrastResp(x,y,expt(day1_id).frame_rate),...
    contrastResp,contrasts,'unif',0);

%% normalized data
fit_norm = cellfun(@(x) x./max(x), fit,'unif',0);
c50_norm = cellfun(@(x) getC50(x),fit_norm,'unif',0);

%% find contrast responsive cells
fitThreshold = 0.6;
[isResponsive_any,isResponsive_eaStim] = cellfun(@(x) findRespCells_anyContrast(...
    x,expt(day1_id).frame_rate),contrastRespData,'unif',0);
maxConIsMaxCon = cellfun(@(x) max(x)==x(end,:)| max(x)==x(end-1,:),...
    contrastResp,'unif',0);
contrastModulatedCells = cellfun(@(x,y,z) x&y'&z>fitThreshold,...
    maxConIsMaxCon,isResponsive_any,R2,'unif',0);
%% day to day correlation
day2dayConCorr = nan(1,sum(goodCells));
for icell = 1:sum(goodCells)
    day2dayConCorr(icell) = corr(contrastResp{1}(:,icell),contrastResp{2}(:,icell));
end

%% save analysis
% save(fullfile(fnout,'contrastResponses'),