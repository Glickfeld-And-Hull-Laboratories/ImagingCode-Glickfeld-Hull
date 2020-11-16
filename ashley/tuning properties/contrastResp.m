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
data = cell(1,nruns);
data_mw = cell(1,nruns);
for irun = 1:nruns
    load(fullfile(fn,runs{irun},'timecourses_cells.mat'))
    data{irun} = tc_subnp_green;
    data_mw{irun} = loadMworksFile(mouse,expDate,times{irun});
end
%% get contrast response tuning

%% fit contrast response

%% save analysis