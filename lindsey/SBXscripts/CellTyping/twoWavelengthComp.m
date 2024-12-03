pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
mouse = 'i3317';
date = '240911';

ImgFolder = {'001','002'};
Run = {'000','001','002'};
ImgType = {'flexGCaMP&flpHT','GCaMP&flpHT&flexMCh'};
RunType = {'920','1040','800'};
Channel = {'Green','Red'};
Reg = [1 2 1];
avgImg = zeros(512,796,size(RunType,2),size(ImgFolder,2),2);

for iImg = 1:size(ImgFolder,2)
    CD = fullfile(pn, mouse, date, ImgFolder{iImg});
    cd(CD);
    figure;
    for iRun = 1:size(Run,2)
        imgMatFile = [ImgFolder{iImg} '_000_' Run{iRun} '.mat'];
        load(imgMatFile);
        totframes = info.config.frames;
        data = sbxread([ImgFolder{iImg} '_000_' Run{iRun}],0,totframes);
        [out data_reg] = stackRegister(squeeze(data(Reg(iRun),:,:,:)), squeeze(mean(data(Reg(iRun),:,:,1:100),4)));
        avgImg(:,:,iRun,iImg,Reg(iRun)) = mean(data_reg,3);
        subplot(3,2,Reg(iRun)+(iRun-1)*2)
        imagesc(avgImg(:,:,iRun,iImg,Reg(iRun)))
        title([RunType{iRun} ' nm- ' Channel{Reg(iRun)}])
        run_use = setdiff(1:2,Reg(iRun));
        [out2 data_reg] = stackRegister_MA(squeeze(data(run_use,:,:,:)),[],[],out);
        avgImg(:,:,iRun,iImg,run_use) = mean(data_reg,3);
        subplot(3,2,run_use+(iRun-1)*2)
        imagesc(avgImg(:,:,iRun,iImg,run_use))
        title([RunType{iRun} ' nm- ' Channel{run_use}])
    end
    sgtitle(ImgType{iImg})
end

fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
mkdir(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{1}]));
mkdir(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{2}]));

for iImg = 1:size(ImgFolder,2)
    for iRun = 1:size(Run,2)
        % subplot(3,3,1+(iRun-1)*3)
        % imagesc(avgImg(:,:,iRun,iImg,1))
        writetiff(avgImg(:,:,iRun,iImg,1), fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_' RunType{iRun} 'nm_' Channel{1} '.tif']))
        % title([RunType{iRun} ' nm- ' Channel{Reg(iRun)}])
        % subplot(3,3,2+(iRun-1)*3)
        % imagesc(avgImg(:,:,iRun,iImg,2))
        writetiff(avgImg(:,:,iRun,iImg,2), fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_' RunType{iRun} 'nm_' Channel{2} '.tif']))
        % title([RunType{iRun} ' nm- ' Channel{2}])
        %subplot(3,3,3+(iRun-1)*3)
        %image(temp)
    end
    sgtitle(ImgType{iImg})
end


for iImg = 1:size(ImgFolder,2)
    figure;
    for iRun = 1:size(Run,2)
        temp = zeros(512,796,3);
        temp1 = zeros(512,796,3);
        subplot(3,3,1+(iRun-1)*3)
        max_use =max(max(max(avgImg(:,:,iRun,iImg,:),[],1),[],2),[],5);
        temp(:,:,2) = avgImg(:,:,iRun,iImg,1)./max_use;
        imagesc(temp)
        title([RunType{iRun} ' nm- ' Channel{Reg(iRun)}])
        subplot(3,3,2+(iRun-1)*3)
        temp(:,:,1) = avgImg(:,:,iRun,iImg,2)./max_use;
        temp1(:,:,1) = avgImg(:,:,iRun,iImg,2)./max_use;
        imagesc(temp1)
        title([RunType{iRun} ' nm- ' Channel{2}])
        subplot(3,3,3+(iRun-1)*3)
        imagesc(temp)
    end
    sgtitle([ImgType{iImg} '- norm each wavelength'])
    print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgFolder{iImg} 'FOV_normEachWL.pdf']),'-dpdf','-bestfit')
end

for iImg = 1:size(ImgFolder,2)
    figure;
    for iRun = 1:size(Run,2)
        temp = zeros(512,796,3);
        temp1 = zeros(512,796,3);
        subplot(3,3,1+(iRun-1)*3)
        max_use =max(max(max(avgImg(:,:,iRun,iImg,1),[],1),[],2),[],5);
        temp(:,:,2) = avgImg(:,:,iRun,iImg,1)./max_use;
        imagesc(temp)
        title([RunType{iRun} ' nm- ' Channel{Reg(iRun)}])
        subplot(3,3,2+(iRun-1)*3)
        max_use =max(max(max(avgImg(:,:,iRun,iImg,2),[],1),[],2),[],5);
        temp(:,:,1) = avgImg(:,:,iRun,iImg,2)./max_use;
        temp1(:,:,1) = avgImg(:,:,iRun,iImg,2)./max_use;
        imagesc(temp1)
        title([RunType{iRun} ' nm- ' Channel{2}])
        subplot(3,3,3+(iRun-1)*3)
        imagesc(temp)
    end
    sgtitle([ImgType{iImg} '- norm each channel'])
    print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgFolder{iImg} 'FOV_normEachCh.pdf']),'-dpdf','-bestfit')
end



Red800nm = avgImg(:,:,3,2,2);
Red1040nm = avgImg(:,:,2,2,2);
Green800nm = avgImg(:,:,3,2,1);
Red920nm = avgImg(:,:,1,2,2);
Green920nm = avgImg(:,:,1,2,1);
Green1040nm = avgImg(:,:,2,2,1);
[out Red800nm_reg] = stackRegister(Red800nm,Red1040nm);
Red = cat(3,Red1040nm,Red800nm_reg,zeros(size(Red1040nm)));
[out2 Green800nm_reg] = stackRegister_MA(Green800nm,[],[],out);
[out3 Red920nm_reg] = stackRegister(Red920nm,Red1040nm);
[out4 Green920nm_reg] = stackRegister_MA(Green920nm,[],[],out3);

iImg = 2;
writetiff(Red1040nm, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_1040nm_Red.tif']))
writetiff(Green1040nm, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_1040nm_Green.tif']))
writetiff(Red920nm_reg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_920nm_Red.tif']))
writetiff(Red800nm_reg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_800nm_Red.tif']))
writetiff(Green920nm_reg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_920nm_Green.tif']))
writetiff(Green800nm_reg, fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgType{iImg} '_800nm_Green.tif']))

figure;
iImg = 2;
iRun = 2;
subplot(2,3,1)
max_use =max(Red(:));
temp = zeros(512,796,3);
temp(:,:,1) = Red1040nm./max(Red1040nm(:));
imagesc(2.*temp)
title([RunType{iRun} ' nm- ' Channel{2}])
iRun = 3;
subplot(2,3,2)
temp1 = zeros(512,796,3);
temp1(:,:,1) = Red800nm_reg./max(Red800nm_reg(:));
imagesc(2.*temp1)
title([RunType{iRun} ' nm- ' Channel{2}])
subplot(2,3,3)
imagesc(2.*Red./max_use)
title('1040 (red) vs 800 (green)')
subplot(2,3,4)
iRun = 1;
temp = zeros(512,796,3);
temp(:,:,1) = Red920nm_reg./max(Red920nm_reg(:));
imagesc(temp)
title([RunType{iRun} ' nm- ' Channel{2}])
subplot(2,3,5)
temp = zeros(512,796,3);
temp(:,:,2) = Green920nm./max(Green920nm(:));
imagesc(temp)
title([RunType{iRun} ' nm- ' Channel{1}])
sgtitle(ImgType{iImg})
print(fullfile(fn_out,[date '_' mouse],[date '_' mouse '_' ImgFolder{iImg}],[date '_' mouse '_' ImgFolder{iImg} 'FOV_800vs1040.pdf']),'-dpdf','-bestfit')
