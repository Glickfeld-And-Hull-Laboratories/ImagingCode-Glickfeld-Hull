function batch_step4b_DART()
% Batch wrapper for step4b_DART_retOnly over multiple sessions.
%
% REQUIRED CHANGES TO step4b_DART_retOnly.m:
%   1. Replace 'clear all' with: clearvars -except instr sessNum
%   2. After run(instr), add:
%        if exist('sessNum','var'), instructions.session = num2str(sessNum); end
%
% IN YOUR INSTR FILE:
%   Use instructions.sessions (plural) with a list, e.g.:
%        instructions.sessions = [2, 4, 6, 7];
%
% Skips sessions where fitGeoTAf is missing from multiday_alignment.mat.

instrFile = 'Z:\home\celine\Repositories\ImagingCode-Glickfeld-Hull\celine\batch_list.m';

rc = behavConstsDART;

% Read sessions list and shared config from instr file
run(instrFile);
run(instructions.ds);
sessionNums = instructions.sessions;

for iSess = 1:length(sessionNums)
    sessNum = sessionNums(iSess);
    fprintf('\n=== Session %d/%d (day_id %d) ===\n', iSess, length(sessionNums), sessNum);

    % Check for retinotopy data
    if expt(sessNum).multiday_timesincedrug_hours > 0
        dart_str = [expt(sessNum).drug '_' num2str(expt(sessNum).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi  = fullfile(rc.analysis, expt(sessNum).exptType, expt(sessNum).mouse, ['multiday_' dart_str]);
    alignFile = fullfile(fn_multi, 'multiday_alignment.mat');

    hasRetino = false;
    if exist(alignFile, 'file')
        alignVars = whos('-file', alignFile);
        if any(strcmp({alignVars.name}, 'fitGeoTAf'))
            tmp       = load(alignFile, 'fitGeoTAf');
            hasRetino = ~isempty(tmp.fitGeoTAf);
        end
    end
    if ~hasRetino
        fprintf('fitGeoTAf not found in %s — skipping\n', alignFile);
        continue
    end

    try
        run_retino_session(instrFile, sessNum);
        fprintf('Session %d/%d complete\n', iSess, length(sessionNums));
    catch ME
        fprintf('ERROR in session %d (day_id %d):\n  %s\n', iSess, sessNum, ME.message);
        if ~isempty(ME.stack)
            fprintf('  In: %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
        end
    end
end

fprintf('\n=== Batch complete: %d sessions ===\n', length(sessionNums));
end


function run_retino_session(instrFile, sessNum) %#ok<INUSD>
% sessNum is intentionally unused here — it survives clearvars in step4b
% because it's named in 'clearvars -except instr sessNum'.
instr = instrFile; %#ok<NASGU>
run('step4b_DART_retOnly');
end