function getOriTuningFSAV(ds,cellsOrDendrites,varargin)
%cellsOrDendrites: 1 == cells; 2 == dendrites; 3 == axons
rc = behavConstsAV;
if isempty(ds) & strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
    isFSAV = 0;
elseif (strcmp(ds(1:3),'FSA')|strcmp(ds(1:3),'FSV')) & strcmp(rc.name,'ashle')
    dataGroup = ds;
    isFSAV = 1;
elseif strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
    isFSAV = 1;
else
    dataGroup = [];
    isFSAV = 0;
end
eval(dataGroup)

if ~isempty(varargin)
    expt = expt(varargin{1});
end

tuningDownSampFactor = 10;
nexp = size(expt,2);

for iexp = 1:nexp
    disp(iexp)
    dataPath = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,...
        'two-photon imaging', expt(iexp).date, expt(iexp).dirtuning);
    if isFSAV
        if cellsOrDendrites == 1
            load(fullfile(dataPath,'timecourses_tun_cells.mat'))
            try
                data_tc_subnp = data_tun_tc_subnp;
            catch
                try
                    data_tc_subnp = data_g_tc_subnp;
                catch
                    data_tc_subnp = dataTC_npSub;
                end
            end
        elseif cellsOrDendrites == 2
            load(fullfile(dataPath,'timecourses_tun_dendrites.mat'))
            data_tc_subnp = data_tun_den_tc_subnp;
        elseif cellsOrDendrites == 3
            load(fullfile(dataPath,'timecoures.mat'))
        end
    else
        load(fullfile(dataPath,'timecourses.mat'))
    end
    mworks = loadMworksFile(...
        expt(iexp).SubNum,expt(iexp).date,expt(iexp).dirtuning_time);
    
    [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuning(data_tc_subnp,mworks,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));
    if cellsOrDendrites == 1
        save(fullfile(dataPath,'oriTuningAndFits'),...
            'avgResponseEaOri','semResponseEaOri',...
            'vonMisesFitAllCells','fitReliability','R_square')
        save(fullfile(dataPath,'oriTuningTCs'),'tuningTC')
    elseif cellsOrDendrites == 2
        save(fullfile(dataPath,'oriTuningAndFits_den'),...
            'avgResponseEaOri','semResponseEaOri',...
            'vonMisesFitAllCells','fitReliability','R_square')
        save(fullfile(dataPath,'oriTuningTCs_den'),'tuningTC')
    elseif cellsOrDendrites == 3
        save(fullfile(dataPath,'oriTuningAndFits_boutons'),...
            'avgResponseEaOri','semResponseEaOri',...
            'vonMisesFitAllCells','fitReliability','R_square')
        save(fullfile(dataPath,'oriTuningTCs_boutons'),'tuningTC')
    end
end