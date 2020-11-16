function getOriTuningFSAV_2color(ds,cellsOrDendrites,varargin)
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
        else
            error('dendrites not available')
        end
    else
        load(fullfile(dataPath,'timecourses.mat'))
    end
    load(fullfile(rc.ashleyData,expt(iexp).mouse,...
        'two-photon imaging', expt(iexp).date, expt(iexp).dirtuning,...
        [expt(iexp).dirtuning '_000_000.mat']))
    if size(data_g_tc,1) == info.config.frames
        data_g_tc_subnp = stackGroupProject_timecourses(data_g_tc_subnp,tuningDownSampFactor);
        data_r_tc_subnp = stackGroupProject_timecourses(data_r_tc_subnp,tuningDownSampFactor);
    end
    mworks = loadMworksFile(...
        expt(iexp).SubNum,expt(iexp).date,expt(iexp).dirtuning_time);
    
    [avgResponseEaOri_g,semResponseEaOri_g,vonMisesFitAllCellsAllBoots_g,...
        fitReliability_g,R_square_g,tuningTC_g] = ...
    getOriTuning(double(data_g_tc_subnp),mworks,tuningDownSampFactor);
    vonMisesFitAllCells_g = squeeze(vonMisesFitAllCellsAllBoots_g(:,1,:));
    [avgResponseEaOri_r,semResponseEaOri_r,vonMisesFitAllCellsAllBoots_r,...
        fitReliability_r,R_square_r,tuningTC_r] = ...
    getOriTuning(double(data_r_tc_subnp),mworks,tuningDownSampFactor);
    vonMisesFitAllCells_r = squeeze(vonMisesFitAllCellsAllBoots_r(:,1,:));
    if cellsOrDendrites == 1
        save(fullfile(dataPath,'oriTuningAndFits_g'),...
            'avgResponseEaOri_g','semResponseEaOri_g',...
            'vonMisesFitAllCells_g','fitReliability_g','R_square_g')
        save(fullfile(dataPath,'oriTuningTCs_g'),'tuningTC_g')
        save(fullfile(dataPath,'oriTuningAndFits_r'),...
            'avgResponseEaOri_r','semResponseEaOri_r',...
            'vonMisesFitAllCells_r','fitReliability_r','R_square_r')
        save(fullfile(dataPath,'oriTuningTCs_r'),'tuningTC_r')
    else
        error('dendrites not available')
    end
end