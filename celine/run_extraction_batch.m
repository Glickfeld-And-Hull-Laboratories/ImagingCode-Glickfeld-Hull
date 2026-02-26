% Batch runner for step2_singleday_extraction.m
% Reads a text file listing a datasheet name and session numbers,
% then runs the extraction script on each session in turn.
% Pupil data is included automatically if pupil.mat is found for that session,
% otherwise it is skipped with a note.
%
% sessions.txt format:
%   Line 1:   datasheet filename (no .m extension)
%   Lines 2+: one session number per line
%   Blank lines and lines starting with # are ignored.

sessions_file = input('Enter sessions filename (no .m extension): ', 's');

% --- Parse sessions file ---
fid = fopen([sessions_file '.txt'], 'r');
if fid == -1
    error('Could not open file: %s.txt', sessions_file);
end

ds       = [];
sessions = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || strncmp(line, '#', 1)
        continue
    end
    if isempty(ds)
        ds = line;
    else
        sessions(end+1) = str2double(line); %#ok<SAGROW>
    end
end
fclose(fid);

if isempty(ds) || isempty(sessions)
    error('Sessions file must contain a datasheet name and at least one session number');
end

fprintf('Datasheet: %s\n', ds);
fprintf('Sessions:  %s\n\n', num2str(sessions));

% Load datasheet and rc paths once so we can check for pupil.mat per session
run(ds);
rc = behavConstsDART;
dataStructLabels = {'contrastxori'};

% --- Loop through sessions ---
failed          = [];
pupil_missing   = [];
pupil_found     = [];

for i = 1:length(sessions)
    day_id = sessions(i);
    fprintf('=== Session %d (%d of %d) ===\n', day_id, i, length(sessions));

    % Construct path to this session's output folder
    mouse            = expt(day_id).mouse;
    expDate          = expt(day_id).date;
    ExperimentFolder = expt(day_id).exptType;
    runs             = eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']);
    runFolder        = runs{1};
    fnout_check      = fullfile(rc.analysis, ExperimentFolder, mouse, expDate, runFolder);

    % Check for pupil data and set includePupil accordingly
    if exist(fullfile(fnout_check, 'pupil.mat'), 'file')
        includePupil = 'y';
        pupil_found(end+1) = day_id; %#ok<SAGROW>
        fprintf('pupil.mat found - including pupil data\n');
    else
        includePupil = 'n';
        pupil_missing(end+1) = day_id; %#ok<SAGROW>
        fprintf('NOTE: pupil.mat not found - skipping pupil analysis\n');
    end

    try
        step2_singleday_extraction;
        fprintf('Completed session %d\n\n', day_id);
    catch ME
        fprintf('FAILED session %d: %s\n\n', day_id, ME.message);
        failed(end+1) = day_id; %#ok<SAGROW>
    end

    % Clear session-specific variables before next iteration
    clearvars -except ds sessions rc expt dataStructLabels failed pupil_missing pupil_found i
end

% --- Summary ---
fprintf('==============================\n');
fprintf('Batch complete: %d/%d succeeded\n', length(sessions)-length(failed), length(sessions));
if ~isempty(pupil_found)
    fprintf('Pupil data included:  sessions %s\n', num2str(pupil_found));
end
if ~isempty(pupil_missing)
    fprintf('Pupil data missing:   sessions %s\n', num2str(pupil_missing));
end
if ~isempty(failed)
    fprintf('Failed sessions:      sessions %s\n', num2str(failed));
end
