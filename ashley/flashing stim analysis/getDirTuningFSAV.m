function getDirTuningFSAV(ds,cellsOrDendrites,varargin)
rc = behavConstsAV;
if isempty(ds) & strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
    isFSAV = 0;
elseif strcmp(ds(1:3),'FSA') & strcmp(rc.name,'ashle')
    dataGroup = ds;
    isFSAV = 1;
elseif strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
    isFSAV = 0;
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
            data_tc_subnp = data_tun_tc_subnp;
        elseif cellsOrDendrites == 2
            load(fullfile(dataPath,'timecourses_tun_dendrites.mat'))
            data_tc_subnp = data_tun_den_tc_subnp;
        end
    else
        load(fullfile(dataPath,'timecourses.mat'))
    end
    mworks = loadMworksFile(...
        expt(iexp).SubNum,expt(iexp).date,expt(iexp).dirtuning_time);
    
    [avgResponseEaDir,semResponseEaDir,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getDirTuning(data_tc_subnp,mworks,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));
    if cellsOrDendrites == 1
        save(fullfile(dataPath,'dirTuningAndFits'),...
            'avgResponseEaDir','semResponseEaDir',...
            'vonMisesFitAllCells','fitReliability','R_square')
        save(fullfile(dataPath,'dirTuningTCs'),'tuningTC')
    elseif cellsOrDendrites == 2
        save(fullfile(dataPath,'dirTuningAndFits_den'),...
            'avgResponseEaDir','semResponseEaDir',...
            'vonMisesFitAllCells','fitReliability','R_square')
        save(fullfile(dataPath,'dirTuningTCs_den'),'tuningTC')
    end
end
end