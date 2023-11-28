datasetStr = 'FS_V1_passive_LG';

rc = behavConstsAV;
eval(datasetStr)  
fnbase = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ashley','data');

for iexp = 1:size(expt,2)
    disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])
    nrun = size(expt(iexp).runs,1);

    for irun = 1:nrun
        ImgFolder = expt(iexp).runs(irun,:);
        fnout = fullfile(fnbase, expt(iexp).mouse, 'two-photon imaging',expt(iexp).date,ImgFolder);
        if ~exist(fnout)
            fn = ['X:\\All_staff\home\ashley\data\' expt(iexp).mouse '\two-photon imaging\' expt(iexp).date '\' ImgFolder];
            fName = [ImgFolder '_000_000'];
            imgMatFile = [fName '.mat'];
            load(fullfile(fn,imgMatFile));
            mkdir(fnout)
            save(fullfile(fnout,imgMatFile),'info')
        end
    end
end
