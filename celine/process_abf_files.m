function [] = process_abf_files(fileList,outpath)
% PROCESS_ABF_FILES Processes multiple .abf files for EPSC detection
%
% This function takes a list of filenames (without extensions) and processes
% each one through the EPSC_detect_fxn.
%
% Usage:
%   process_abf_files({'24110000', '24110001', '24110002'})
%   - or -
%   process_abf_files('filenames.txt') % Text file with one filename per line

% Input validation
if ischar(fileList) || isstring(fileList)
    % If input is a string, assume it's a text file with filenames
    if exist(fileList, 'file')
        fid = fopen(fileList, 'r');
        fileData = textscan(fid, '%s');
        fclose(fid);
        fileList = fileData{1};
        fprintf('Loaded %d filenames from %s\n', length(fileList), char(fileList(1:min(length(fileList),3))'));
    else
        error('File list not found: %s', fileList);
    end
elseif ~iscell(fileList)
    error('Input must be a cell array of filenames or a text file containing filenames');
end

% Process each file in the list
fprintf('Processing %d files...\n', length(fileList));
for i = 1:length(fileList)
    fileName = fileList{i};
    fprintf('\n==== Processing file %d of %d: %s ====\n', i, length(fileList), fileName);
    
    try
        % Check if the file exists
        if ~exist([fileName, '.abf'], 'file')
            fprintf('ERROR: File not found: %s.abf\n', fileName);
            continue;
        end
        
        % Process the current file with EPSC_detect_fxn
        EPSC_detect_fxn(fileName,outpath);
        fprintf('Successfully processed %s\n', fileName);
        
    catch ME
        % Handle any errors during processing
        fprintf('ERROR processing %s: %s\n', fileName, ME.message);
    end
end

fprintf('\nAll files processed.\n');
close all
end