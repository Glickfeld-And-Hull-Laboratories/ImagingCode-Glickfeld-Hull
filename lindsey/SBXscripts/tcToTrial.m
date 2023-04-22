function [data_tr nOn nOff] = tcToTrial(data,input);
    % transform timecourse into nframes x ntrials
    nOn= input.nScansOn;
    nOff = input.nScansOff;
    nTrial = size(input.tGratingContrast);
    sz = size(data);
    data_tr = reshape(data,[nOn+nOff, nTrial, min(sz)]);
end