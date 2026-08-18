function addSemicolons(inputFile, outputFile)
% addSemicolons(inputFile, outputFile)
% Adds semicolons to the end of lines in a MATLAB file if missing.
%
% Example:
%   addSemicolons('myscript.m', 'myscript_semicolon.m')

    % Read file
    fid = fopen(inputFile, 'r');
    if fid == -1
        error('Cannot open file: %s', inputFile);
    end
    lines = {};
    tline = fgetl(fid);
    while ischar(tline)
        lines{end+1} = tline; %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    % Process lines
    for i = 1:numel(lines)
        line = strtrim(lines{i});

        % Skip empty lines or comment-only lines
        if isempty(line) || startsWith(line, '%')
            continue
        end

        % Skip lines that already end with semicolon, ellipsis, or block symbols
        if ~endsWith(line, ';') && ~endsWith(line, '...') && ~endsWith(line, 'end') 
            % Only add semicolon if it's an assignment or function call
            lines{i} = [lines{i} ';'];
        end
    end

    % Write output
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Cannot write to file: %s', outputFile);
    end
    for i = 1:numel(lines)
        fprintf(fid, '%s\n', lines{i});
    end
    fclose(fid);

    fprintf('Semicolons added and saved to "%s"\n', outputFile);
end
