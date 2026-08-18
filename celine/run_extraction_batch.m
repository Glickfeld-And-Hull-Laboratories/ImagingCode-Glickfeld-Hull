% Batch runner for step2_singleday_extraction.m
% Reads a text file listing a datasheet name and session numbers,
% then runs the extraction script on each session in turn.
%
% sessions.txt format:
%   Line 1:   datasheet filename (no .m extension)
%   Lines 2+: one session number per line
%   Blank lines and lines starting with # are ignored.
%
% Example sessions.txt:
%   DART_V1_YM90K_Celine
%   1
%   2
%   5
%   6

sessions_file = input('Enter sessions filename (no .m extension): ', 's');
includePupil  = input('Include pupil data for all sessions? (y/n): ', 's');
retino_resp   = input('Run retinotopy alignment for all sessions? (y/n): ', 's');
doRetino      = strcmpi(retino_resp, 'y') || strcmpi(retino_resp, 'yes');
if doRetino
    val_resp         = input('Plot retinotopy validation traces? (y/n): ', 's');
    validation_choice = strcmpi(val_resp, 'y') || strcmpi(val_resp, 'yes');
else
    validation_choice = false;
end

% --- Parse sessions file ---
fid = fopen([sessions_file '.txt'], 'r');
if fid == -1
    error('Could not open file: %s.txt', sessions_file);
end

ds       = [];
sessions = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '#')
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
    error('sessions file must contain a datasheet name and at least one session number');
end

fprintf('Datasheet:  %s\n', ds);
fprintf('Sessions:   %s\n', num2str(sessions));
fprintf('Pupil:      %s\n', includePupil);
fprintf('Retinotopy: %s\n\n', retino_resp);

% --- Loop through sessions ---
failed = [];

for i = 1:length(sessions)
    day_id = sessions(i);
    fprintf('=== Session %d (%d of %d) ===\n', day_id, i, length(sessions));

    try
        step2_singleday_extraction;
        fprintf('Completed session %d\n\n', day_id);
    catch ME
        fprintf('FAILED session %d: %s\n\n', day_id, ME.getReport);
        failed(end+1) = day_id; %#ok<SAGROW>
    end

    % Clear session-specific variables before next iteration
    clearvars -except ds sessions includePupil doRetino validation_choice failed i
end

% --- Summary ---
fprintf('==============================\n');
fprintf('Batch complete: %d/%d succeeded\n', length(sessions)-length(failed), length(sessions));
if ~isempty(failed)
    fprintf('Failed sessions: %s\n', num2str(failed));
end