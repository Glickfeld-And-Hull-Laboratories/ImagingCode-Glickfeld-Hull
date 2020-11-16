%% which dataset and parameters
clear all
clear global
close all
ds = 'multiday_test'; %dataset info
dataStructLabels = {'contrast'};
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doGreenOnly = true;
%% get expt info

mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');

%% load time-courses
runs = expt(slct_expt).contrast_runs;
times = expt(slct_expt).contrast_time;
nruns = length(runs);

if nruns > 1
    conditions = expt(slct_expt).contrast_condition;

data = cell(1,nruns);
data_mw = cell(1,nruns);
for irun = 1:nruns
    load(fullfile(fn,runs{irun},'timecourses_cells.mat'))
    data{irun} = tc_subnp_green;
    data_mw{irun} = loadMworksFile(mouse,expDate,times{irun});
end
%% get contrast response tuning
[contrastResp,contrastRespErr,~,~,contrasts,contrastRespData] = cellfun(@(x,y) getContrastResp(...
    x,y,expt(slct_expt).frame_rate),data,data_mw,'unif',0);
%% fit contrast response
[fit,R2,c50] = cellfun(@(x,y) fitContrastResp(x,y,expt(slct_expt).frame_rate),...
    contrastResp,contrasts,'unif',0);

%% find contrast responsive cells
[isResponsive_any,isResponsive_eaStim] = cellfun(@(x) findRespCells_anyContrast(...
    x,expt(slct_expt).frame_rate),contrastRespData,'unif',0);
%% save analysis
save(fullfile(fnout,'contrastResponses'),