function nfr_run = nFramesSbxDataset(mouse,expDate,runFolder)

if strcmp(expDate,'150508') & strcmp(runFolder,'003')
    nfr_run = 13577;
else

    fn = ['Z:\home\ashley\data\' mouse '\two-photon imaging\' expDate '\' runFolder];
    fName = [runFolder '_000_000'];
    imgMatFile = [fName '.mat'];
    load(fullfile(fn,imgMatFile));
    nfr_run = info.config.frames;
end
end
