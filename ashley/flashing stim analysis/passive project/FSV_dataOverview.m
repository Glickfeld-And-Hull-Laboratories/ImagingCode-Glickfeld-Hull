clear all
close all
ds = 'FSV_V1';
eval(ds)

%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

fnout = fullfile(rc.caOutputDir,ds,[ds '_']);


load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));

%%
% rewwin = 32:39;
% nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
respwin_opt = respwin;
respwin_target_opt = respwin_target;
%% count up number of mice and experiments from each mouse for certain
% categories

% sort PV and SOM expts
pv_expt = [];
som_expt = [];
for iexp = 1:nexp
    if strcmp(expt(iexp).redChannelLabel,'PV') || strcmp(expt(iexp).greenChannelLabel,'PV')
        if isempty(pv_expt)
            pv_expt = expt(iexp);
        else
            pv_expt(length(pv_expt)+1) = expt(iexp);
        end
    end
    if strcmp(expt(iexp).redChannelLabel,'SOM') || strcmp(expt(iexp).greenChannelLabel,'SOM')
        if isempty(som_expt)
            som_expt = expt(iexp);
        else
            som_expt(length(som_expt)+1) = expt(iexp);
        end
    end
end

% counts in PV expts
pv_info = struct;
pv_info.mice = unique({pv_expt.SubNum});
pv_info.nExpts = nan(1,length(pv_info.mice));
for i = 1:length(pv_info.mice)
    pv_info.nExpts(i) = length(pv_expt(strcmp(pv_info.mice{i},{pv_expt.SubNum})));
end
pv_info.isTwoColor = cellfun(@isempty,{pv_expt.redChannelLabel}) == 0;
pv_info.indicator = cellfun(@(x) x{2}(end-6:end),{pv_expt.indicator},'unif',0);
pv_info.isBehav = cell2mat({pv_expt.isBehav}) == 1;
pv_info.isAttnTask = cell2mat_padded({pv_expt.attentionTask}) == 1;
pv_info.isAttn = cell2mat_padded({pv_expt.hasAttention}) == 1;
pv_info.isPassive = cellfun(@isempty,{pv_expt.passExpt}) == 0;

% counts in SOM expts
som_info = struct;
som_info.mice = unique({som_expt.SubNum});
som_info.nExpts = nan(1,length(som_info.mice));
for i = 1:length(som_info.mice)
    som_info.nExpts(i) = length(som_expt(strcmp(som_info.mice{i},{som_expt.SubNum})));
end
som_info.isTwoColor = cellfun(@isempty,{som_expt.redChannelLabel}) == 0;
som_info.indicator = cellfun(@(x) x{2}(end-6:end),{som_expt.indicator},'unif',0);
som_info.isBehav = cell2mat({som_expt.isBehav}) == 1;
som_info.isAttnTask = cell2mat_padded({som_expt.attentionTask}) == 1;
som_info.isAttn = cell2mat_padded({som_expt.hasAttention}) == 1;
som_info.isPassive = cellfun(@isempty,{som_expt.passExpt}) == 0;
%% print out some info about som and pv experiments
ind = som_info.isBehav & som_info.isTwoColor;
fprintf('\n SOM: %s behavior mice have %s two color expts',...
    num2str(length(unique({som_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = som_info.isBehav & som_info.isPassive & som_info.isTwoColor;
fprintf('\n SOM: %s behavior mice have %s passive expts',...
    num2str(length(unique({som_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = som_info.isBehav & ~som_info.isTwoColor;
fprintf('\n SOM: %s behavior mice have %s gcamp-only expts',...
    num2str(length(unique({som_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = ~som_info.isBehav & som_info.isTwoColor;
fprintf('\n SOM: %s naive mice have %s two color expts\n',...
    num2str(length(unique({som_expt(ind).SubNum}))),...
    num2str(sum(ind)))

ind = pv_info.isBehav & pv_info.isTwoColor;
fprintf('\n PV: %s behavior mice have %s two color expts',...
    num2str(length(unique({pv_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = pv_info.isBehav & pv_info.isPassive & pv_info.isTwoColor;
fprintf('\n PV: %s behavior mice have %s passive expts',...
    num2str(length(unique({pv_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = pv_info.isBehav & ~pv_info.isTwoColor;
fprintf('\n PV: %s behavior mice have %s gcamp-only expts',...
    num2str(length(unique({pv_expt(ind).SubNum}))),...
    num2str(sum(ind)))
ind = ~pv_info.isBehav & pv_info.isTwoColor;
fprintf('\n PV: %s naive mice have %s two color expts\n',...
    num2str(length(unique({pv_expt(ind).SubNum}))),...
    num2str(sum(ind)))
%% data overview - engagement


%% SOM analyses
% edit FSV_analysis_adaptation
edit FSV_analysis_visStimAlignedTCs.m
edit FSV_analysis_anticipationOverview.m
edit FSV_analysis_spatial.m
edit FSV_analysis_decoding.m
%% other analyses
edit FSV_analysis_tuning