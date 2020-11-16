clear all
clear global
close all
ds = 'atro_adaptcon_V1'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doPreviousReg = false;
doGreenOnly = false;
%%
mouse = expt(slct_expt).mouse;
expDate = expt(slct_expt).date;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
fnout = fullfile(fn,'data processing');
%% load and register data
runFolder = expt(slct_expt).texasRedRun;
fName = [runFolder '_000_001'];
data_temp= loadsbx_choosepmt(3,mouse,expDate,runFolder,fName);
data_g_temp = squeeze(data_temp(1,:,:,:));
data_r_temp = squeeze(data_temp(2,:,:,:));
clear data_temp

load(fullfile(fnout,'regImg'))
[outs,data_g_reg] = stackRegister(data_g_temp,regImg);
[~,data_r_reg] = stackRegister_MA(data_r_temp,[],[],double(outs));

writetiff(data_g_reg,fullfile(fnout,'texasRed_greenCh'))
writetiff(data_r_reg,fullfile(fnout,'texasRed_redCh'))

% no stim runs
nruns = length(expt(slct_expt).nostim_runs);
data_nostim = cell(1,nruns);
for irun = 1:nruns
    
    runFolder = expt(slct_expt).nostim_runs{irun};
    mkdir(fullfile(fn,runFolder))
    fName = [runFolder '_000_000'];
        
    if expt(slct_expt).greenredsimultaneous == 1
        data_temp_g = loadsbx_choosepmt(3,mouse,expDate,runFolder,fName);
        data_temp_r = squeeze(data_temp_r(2,:,:,:));
        data_temp_g = squeeze(data_temp_g(1,:,:,:));
    else
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
    end
    % remove negative data by subtraction
    data_sub = data_temp_g-min(min(min(data_temp_g,[],1),[],2),[],3);
    data_temp_g = data_sub;
    if expt(slct_expt).greenredsimultaneous == 1
        data_sub = data_temp_r-min(min(min(data_temp_r,[],1),[],2),[],3);
        data_temp_r = data_sub;
    end
    clear data_sub
    
    % reg images
    if ~doPreviousReg && irun == 1
        regImgStartFrame = compareRegImg_2color(rc,ds,mouse,expDate,data_temp_g);
        regImg = mean(data_temp_g(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
        save(fullfile(fnout,'regImg'),'regImg')
    elseif irun == 1
        load(fullfile(fnout,'regImg'))
    end
    
    if doPreviousReg
        load(fullfile(fn,runFolder,'regOuts&Img'))
        [~,data_g_reg] = stackRegister_MA(data_temp_g,[],[],double(outs));
        figure;imagesc(regImg);colormap gray; title('reg image')
    else
        [outs,data_g_reg] = stackRegister(data_temp_g,regImg);
        figure;imagesc(mean(data_g_reg,3));colormap gray; title(['no vis stim - ' expt(slct_expt).nostim_condition{irun}])
        clear data_temp_g
        save(fullfile(fn,runFolder,'regOuts&Img'),'outs','regImg')
    end
    
    data_nostim{irun} = data_g_reg;
    
    nfr = size(data_g_reg,3);
    nimg = 100;
    frInterval = 1:floor((nfr-100)/nimg):nfr-100;
    F = nan(size(data_g_reg,1),size(data_g_reg,2),nimg);
    for iimg = 1:nimg
        F(:,:,iimg) = mean(data_g_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
    end
    writetiff(F,fullfile(fn,runFolder,'FImages'))
end
clear data_g_reg
