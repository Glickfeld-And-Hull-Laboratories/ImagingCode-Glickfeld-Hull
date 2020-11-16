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
mkdir(fnout)

%% load time-courses and ori tuning
load(fullfile(fnout,'timecourses_cells'))

%% load mworks data
nruns = length(expt(slct_expt).dirtuning_runs);
data_mw = cell(1,nruns);
for irun = 1:nruns    
    mw_temp = loadMworksFile(mouse,expDate,expt(slct_expt).dirtuning_time{irun});
    data_mw{irun} = mw_temp;
end
%% get ori tuning

[avgResponseEaOri_tagneg,semResponseEaOri_tagneg,...
    vonMisesFitAllCells_tagneg,fitReliability_tagneg,~,tuningTC_tagneg,isResp_tagneg] =...
    cellfun(@(x,y) getOriTuning(x,y,dsFactor,expt(slct_expt).frame_rate),tc_tagneg_tuning_subnp,...
    data_mw,'unif',0);
if ~isnan(expt(slct_expt).redChannelRun)
[avgResponseEaOri_tagpos,semResponseEaOri_tagpos,...
    vonMisesFitAllCells_tagpos,fitReliability_tagpos,~,tuningTC_tagpos,isResp_tagpos] =...
    cellfun(@(x,y) getOriTuning(x,y,dsFactor),tc_tagpos_tuning_subnp,...
    data_mw,'unif',0);
end
if ~isnan(expt(slct_expt).redChannelRun)
save(fullfile(fnout,'oriTuning'),'avgResponseEaOri_tagneg',...
    'semResponseEaOri_tagneg','vonMisesFitAllCells_tagneg',...
    'fitReliability_tagneg','tuningTC_tagneg','isResp_tagneg',...
    'avgResponseEaOri_tagpos',...
    'semResponseEaOri_tagpos','vonMisesFitAllCells_tagpos',...
    'fitReliability_tagpos','tuningTC_tagpos','isResp_tagpos');
else
save(fullfile(fnout,'oriTuning'),'avgResponseEaOri_tagneg',...
    'semResponseEaOri_tagneg','vonMisesFitAllCells_tagneg',...
    'fitReliability_tagneg','tuningTC_tagneg','isResp_tagneg');
end