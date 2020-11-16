%%
jb_dir = ['Z:\home\ashley\Manuscripts\Attention V1\Matlab Figs\JB Dataset'];
date = '09-Jul-2020';
% mkdir(fullfile(jb_dir,date))
save(fullfile(jb_dir,date,'noAttn_respCellsExpt_decodingDataExpt'),...
    'respCellsExpt','decodeDataExpt')
%%
rc = behavConstsAV;
imgParams_FSAV
%%
load(fullfile(jb_dir,date,'attn_respCellsExpt_decodingDataExpt'),...
    'respCellsExpt','decodeDataExpt')

data_attn = struct;
for iexp = 1:size(respCellsExpt,2)
    
    data_attn(iexp).expt_name = respCellsExpt(iexp).exptName;
    data_attn(iexp).hasAttention = true;
    
    cellInd = respCellsExpt(iexp).lateCycRespCells |...
        respCellsExpt(iexp).targetRespCells;
    
%     r = cat(2,decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
%         decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:))';
    if isfield(decodeDataExpt(iexp).catch,'cycRespEaTrial')
        r = cat(2,cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{visualTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).catch.cycRespEaTrial,'unif',0));     
        
        nvis = size(decodeDataExpt(iexp).cycRespEaTrial{visualTrials},2);
        naud = size(decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},2);
        ncatch = size(decodeDataExpt(iexp).catch.cycRespEaTrial,2);
        outcome = cat(2,decodeDataExpt(iexp).outcomeEaTrial{visualTrials},...
            decodeDataExpt(iexp).outcomeEaTrial{auditoryTrials},...
            decodeDataExpt(iexp).catch.outcomeEaTrial);
        catchoutcome = cat(2,cell(1,nvis+naud),...
            decodeDataExpt(iexp).catch.catchOutcomeEaTrial);

        ori = cat(2,decodeDataExpt(iexp).stimEaTrial{visualTrials},zeros(1,naud),...
            zeros(1,ncatch));
        amp = cat(2,zeros(1,nvis),decodeDataExpt(iexp).stimEaTrial{auditoryTrials},...
            decodeDataExpt(iexp).catch.stimEaTrial);
        catchstimtime = cat(2,nan(1,nvis),nan(1,naud),...
            decodeDataExpt(iexp).catch.catchCycN);
        catchori = cat(2,nan(1,nvis),nan(1,naud),...
            decodeDataExpt(iexp).catch.catchStimEaTrial);
        
        isvis = cat(2,true(1,nvis),false(1,naud),false(1,ncatch));
    else
        r = cat(2,cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{visualTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},'unif',0)); 
        
        nvis = size(decodeDataExpt(iexp).cycRespEaTrial{visualTrials},2);
        naud = size(decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},2);
        outcome = cat(2,decodeDataExpt(iexp).outcomeEaTrial{visualTrials},...
            decodeDataExpt(iexp).outcomeEaTrial{auditoryTrials});
        catchoutcome = cell(1,nvis+naud);

        ori = cat(2,decodeDataExpt(iexp).stimEaTrial{visualTrials},zeros(1,naud));
        amp = cat(2,zeros(1,nvis),decodeDataExpt(iexp).stimEaTrial{auditoryTrials}); 
        catchstimtime = cat(2,nan(1,nvis),nan(1,naud));
        catchori = cat(2,nan(1,nvis),nan(1,naud));   
        
        isvis = cat(2,true(1,nvis),false(1,naud));
    end
    
            
    data_attn(iexp).dff = r; %(trials x neurons)
    data_attn(iexp).visualTrial = isvis;
    data_attn(iexp).stimulus_orientation = ori;
    data_attn(iexp).stimulus_soundVolume = amp;
    data_attn(iexp).mouseResponse = outcome;
    data_attn(iexp).catchStimulus_orientation = catchori;
    data_attn(iexp).catchStimulus_time = catchstimtime;
    data_attn(iexp).mouseResponse_catchStimulus = catchoutcome;
end

load(fullfile(jb_dir,date,'noAttn_respCellsExpt_decodingDataExpt'),...
    'respCellsExpt','decodeDataExpt')

data_noAttn = struct;
for iexp = 1:size(respCellsExpt,2)
    
    data_noAttn(iexp).expt_name = respCellsExpt(iexp).exptName;
    data_noAttn(iexp).hasAttention = false;
    
    cellInd = respCellsExpt(iexp).lateCycRespCells |...
        respCellsExpt(iexp).targetRespCells;
    
%     r = cat(2,decodeDataExpt(iexp).av(visualTrials).resp(cellInd,:),...
%         decodeDataExpt(iexp).av(auditoryTrials).resp(cellInd,:))';
    if isfield(decodeDataExpt(iexp).catch,'cycRespEaTrial')
        r = cat(2,cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{visualTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).catch.cycRespEaTrial,'unif',0));     
        
        nvis = size(decodeDataExpt(iexp).cycRespEaTrial{visualTrials},2);
        naud = size(decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},2);
        ncatch = size(decodeDataExpt(iexp).catch.cycRespEaTrial,2);
        outcome = cat(2,decodeDataExpt(iexp).outcomeEaTrial{visualTrials},...
            decodeDataExpt(iexp).outcomeEaTrial{auditoryTrials},...
            decodeDataExpt(iexp).catch.outcomeEaTrial);
        catchoutcome = cat(2,cell(1,nvis+naud),...
            decodeDataExpt(iexp).catch.catchOutcomeEaTrial);

        ori = cat(2,decodeDataExpt(iexp).stimEaTrial{visualTrials},zeros(1,naud),...
            zeros(1,ncatch));
        amp = cat(2,zeros(1,nvis),decodeDataExpt(iexp).stimEaTrial{auditoryTrials},...
            decodeDataExpt(iexp).catch.stimEaTrial);
        catchstimtime = cat(2,nan(1,nvis),nan(1,naud),...
            decodeDataExpt(iexp).catch.catchCycN);
        catchori = cat(2,nan(1,nvis),nan(1,naud),...
            decodeDataExpt(iexp).catch.catchStimEaTrial);
        
        isvis = cat(2,true(1,nvis),false(1,naud),false(1,ncatch));
    else
        r = cat(2,cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{visualTrials},'unif',0),...
            cellfun(@(x) x(:,cellInd),...
            decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},'unif',0)); 
        
        nvis = size(decodeDataExpt(iexp).cycRespEaTrial{visualTrials},2);
        naud = size(decodeDataExpt(iexp).cycRespEaTrial{auditoryTrials},2);
        outcome = cat(2,decodeDataExpt(iexp).outcomeEaTrial{visualTrials},...
            decodeDataExpt(iexp).outcomeEaTrial{auditoryTrials});
        catchoutcome = cell(1,nvis+naud);

        ori = cat(2,decodeDataExpt(iexp).stimEaTrial{visualTrials},zeros(1,naud));
        amp = cat(2,zeros(1,nvis),decodeDataExpt(iexp).stimEaTrial{auditoryTrials}); 
        catchstimtime = cat(2,nan(1,nvis),nan(1,naud));
        catchori = cat(2,nan(1,nvis),nan(1,naud));   
        
        isvis = cat(2,true(1,nvis),false(1,naud));
    end
    
            
    data_noAttn(iexp).dff = r; %(trials x neurons)
    data_noAttn(iexp).visualTrial = isvis;
    data_noAttn(iexp).stimulus_orientation = ori;
    data_noAttn(iexp).stimulus_soundVolume = amp;
    data_noAttn(iexp).mouseResponse = outcome;
    data_noAttn(iexp).catchStimulus_orientation = catchori;
    data_noAttn(iexp).catchStimulus_time = catchstimtime;
    data_noAttn(iexp).mouseResponse_catchStimulus = catchoutcome;
end

data = cat(2,data_attn,data_noAttn);
save(fullfile(jb_dir,date,'attnData_neuronResp_trialOutcome'),...
    'data')

%%
clear all
jb_dir = ['Z:\home\ashley\Manuscripts\Attention V1\Matlab Figs\JB Dataset'];
load(fullfile(jb_dir,date,'attnData_neuronResp_trialOutcome'),...
    'data')