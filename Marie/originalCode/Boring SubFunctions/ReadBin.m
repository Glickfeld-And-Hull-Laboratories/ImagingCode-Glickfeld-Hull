function dataArray = ReadBin(ReadMin, ReadMax, meta, binName, path) % read a chunk to smooth and analyze, Batch (in sec) with addtional smoothing edges (in samples)

    nChan = str2double(meta.nSavedChans);

    %nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);

    sizeA = [nChan, (ReadMax-ReadMin)];
    
    fid = fopen(fullfile(path, binName), 'rb');
    fseek(fid, ReadMin * 2 * nChan, 'bof');
     %position1 = ftell(fid) % For troubleshooting
    dataArray = fread(fid, sizeA, 'int16=>double');
     %position2 = ftell(fid)  % For troubleshooting

    fclose(fid);
end % ReadBin