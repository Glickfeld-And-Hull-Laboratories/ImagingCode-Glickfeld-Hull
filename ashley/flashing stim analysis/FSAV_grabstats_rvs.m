clear all
close all
ds = 'FSAV_attentionV1_noAttn'; % 'FSAV_V1_100ms_naive'  'FSAV_V1_naive_GCaMP6m'  'FSAV_attentionV1'   'FSAV_attentionV1_noAttn'
cellsOrDendrites = 1;
attnAnalysisDate = '191211';

%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms


%% load behavior stats

load(fullfile(rc.ashley, 'Manuscripts','Attention V1','Matlab Figs','bxStats.mat'))

%% for reaction time figure
grpnames = {'modality','attention'};

rt_diff = 

avGroup = cat(1,ones(