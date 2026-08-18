% Pool MAD variability outputs across sessions.
% Loads singleday_extraction.mat from each session defined in the session
% list txt file and concatenates across cells, preserving red/green identity.
%
% Session list txt format:
%   Line 1 (uncommented): name of datasheet .m file
%   Remaining uncommented lines: session indices (day_id) to pool

if ~exist('sess_file', 'var')
    sess_file = input('Enter session list file: ', 's');
end

% Parse session list txt file
fid = fopen(sess_file, 'r');
lines = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ischar(line) && ~isempty(line) && ~startsWith(line, '#')
        lines{end+1} = line;
    end
end
fclose(fid);

ds      = lines{1};
day_ids = cellfun(@str2num, lines(2:end));

run(ds);

rc = behavConstsDART;
dataStructLabels = {'contrastxori'};

% Pre-allocate pooled arrays (grow by concatenation)
conBySize_MAD_stat_all       = [];
conBySize_MAD_loc_all        = [];
conBySize_MAD_largePupil_all = [];
conBySize_MAD_smallPupil_all = [];
tc_trial_MAD_stat_all        = [];
tc_trial_MAD_loc_all         = [];
tc_trial_MAD_largePupil_all  = [];
tc_trial_MAD_smallPupil_all  = [];
red_cells_all                = [];
session_id_all               = [];

nOn      = [];
cons     = [];
sizes    = [];
dirs     = [];
resp_win = [];
base_win = [];

for iDay = day_ids
    mouse            = expt(iDay).mouse;
    expDate          = expt(iDay).date;
    ExperimentFolder = expt(iDay).exptType;
    runs             = eval(['expt(iDay).' cell2mat(dataStructLabels) '_runs']);
    runFolder        = runs{1};

    fnout   = fullfile(rc.analysis, ExperimentFolder, mouse, expDate, runFolder);
    matFile = fullfile(fnout, 'singleday_extraction.mat');

    if ~exist(matFile, 'file')
        fprintf('Session %d: singleday_extraction.mat not found, skipping\n', iDay);
        continue
    end

    fprintf('Loading session %d: %s %s\n', iDay, mouse, expDate);
    d = load(matFile, ...
        'conBySize_MAD_stat', 'conBySize_MAD_loc', ...
        'conBySize_MAD_largePupil', 'conBySize_MAD_smallPupil', ...
        'tc_trial_MAD_stat', 'tc_trial_MAD_loc', ...
        'tc_trial_MAD_largePupil', 'tc_trial_MAD_smallPupil', ...
        'red_cells', 'nOn', 'nOff', 'cons', 'sizes', 'dirs', ...
        'resp_win', 'base_win');

    % On first valid session, grab shared params
    if isempty(cons)
        nOn      = d.nOn;
        nOff     = d.nOff;
        cons     = d.cons;   nCon  = length(cons);
        sizes    = d.sizes;  nSize = length(sizes);
        dirs     = d.dirs;
        resp_win = d.resp_win;
        base_win = d.base_win;
    end

    nKeep_day = size(d.conBySize_MAD_stat, 1);

    conBySize_MAD_stat_all       = [conBySize_MAD_stat_all;       d.conBySize_MAD_stat];
    conBySize_MAD_loc_all        = [conBySize_MAD_loc_all;        d.conBySize_MAD_loc];
    conBySize_MAD_largePupil_all = [conBySize_MAD_largePupil_all; d.conBySize_MAD_largePupil];
    conBySize_MAD_smallPupil_all = [conBySize_MAD_smallPupil_all; d.conBySize_MAD_smallPupil];

    tc_trial_MAD_stat_all       = cat(2, tc_trial_MAD_stat_all,       d.tc_trial_MAD_stat);
    tc_trial_MAD_loc_all        = cat(2, tc_trial_MAD_loc_all,        d.tc_trial_MAD_loc);
    tc_trial_MAD_largePupil_all = cat(2, tc_trial_MAD_largePupil_all, d.tc_trial_MAD_largePupil);
    tc_trial_MAD_smallPupil_all = cat(2, tc_trial_MAD_smallPupil_all, d.tc_trial_MAD_smallPupil);

    red_cells_all  = [red_cells_all;  d.red_cells(:)];
    session_id_all = [session_id_all; iDay * ones(nKeep_day, 1)];

    fprintf('  %d cells (cumulative: %d)\n', nKeep_day, length(red_cells_all));
end

nKeep_total = length(red_cells_all);
HTP_cells   = logical(red_cells_all);
fprintf('\nTotal pooled cells: %d (%d HTP+, %d HTP-)\n', ...
    nKeep_total, sum(HTP_cells), sum(~HTP_cells));

%% Diagnostic plot
t_frames = (1:(nOn+nOff)) - nOff/2;

figure;

subplot(1,2,1)
plot(t_frames, mean(tc_trial_MAD_stat_all(:, HTP_cells, 3, 3), 2, 'omitnan'))
title('HTP+  100% contrast 30 deg')
xlim([-30,60])
set(gca, 'TickDir', 'out', 'Box', 'off')

subplot(1,2,2)
plot(t_frames, mean(tc_trial_MAD_stat_all(:, ~HTP_cells, 3, 3), 2, 'omitnan'))
title('HTP- 100% contrast 30 deg'); ylabel('MAD (\DeltaF/F)');
xlim([-30,60])
set(gca, 'TickDir', 'out', 'Box', 'off')



%%
plotNeuralTimecourse_singleDay(tc_trial_MAD_stat_all, HTP_cells, ~HTP_cells, ...
    'Colors', {'r', 'g'}, 'Titles', {'HTP+', 'HTP-'}, ...
    'FrameRate', frame_rate, 'StimStart', stimStart);
figs = findobj('Type', 'figure');
for i = 1:length(figs)
    saveas(figs(i), fullfile(fnout, sprintf('timecourse_stat_size%d.pdf', i)));
end
%% Save
fnpool = fullfile(rc.analysis, ExperimentFolder);
save(fullfile(fnpool, 'pooled_MAD.mat'), ...
    'conBySize_MAD_stat_all', 'conBySize_MAD_loc_all', ...
    'conBySize_MAD_largePupil_all', 'conBySize_MAD_smallPupil_all', ...
    'tc_trial_MAD_stat_all', 'tc_trial_MAD_loc_all', ...
    'tc_trial_MAD_largePupil_all', 'tc_trial_MAD_smallPupil_all', ...
    'red_cells_all', 'HTP_cells', 'session_id_all', ...
    'cons', 'sizes', 'dirs', 'nOn', 'nOff', 'resp_win', 'base_win', ...
    'nKeep_total', '-v7.3');

fprintf('Saved pooled MAD to %s\n', fullfile(fnpool, 'pooled_MAD.mat'));