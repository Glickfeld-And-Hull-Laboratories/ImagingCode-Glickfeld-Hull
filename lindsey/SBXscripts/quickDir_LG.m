clear all
clear all global
close all

date = '230222';
mouse = 'i2902';
ImgFolder = '002';
time = ['1531'];
doReg = 1;
nrun = size(ImgFolder,1);
rc = behavConstsAV;
subnum = mouse;

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

fprintf(['Reading run ' num2str(ImgFolder(1,:))])

data = [];
clear temp
for irun = 1:nrun

    CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images\' mouse '\' date '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    nframes = info.config.frames;
    %nframes = 4050;
    nframes = input.counterValues{end}(end);
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

    if size(data_temp,1) > 1
        data_temp = data_temp(1,:,:,:);
    end
    
    
    expt_input = input;
    temp(irun) = expt_input;
    
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    
    
    if size(data_temp,1) == 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
%     if nframes>ntrials*(nOn+nOff)
%         data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
%     elseif nframes<ntrials*(nOn+nOff)
%         temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
%     end
        
    data = cat(3,data,data_temp);
end

clear data_temp
expt_input = concatenateDataBlocks(temp);

    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);

    nOn = expt_input.nScansOn;
    nOff = expt_input.nScansOff;
    ntrials = size(expt_input.tGratingDirectionDeg,2);
    
   %use if pmt 2 was saved
%     data = squeeze(data(1,:,:,:));
    
    if doReg
    data_avg = mean(data(:,:,1000:1500),3);
    [out data_reg] = stackRegister(data,data_avg);
    data = data_reg;
    clear data_reg
    end

    dir_mat = celleqel2mat_padded(expt_input.tGratingDirectionDeg);
    if (nOn+nOff)*ntrials > size(data,3)
        ntrials = floor(size(data,3)/((nOn+nOff)));    
        dir_mat = dir_mat(1:ntrials);
    end
    data = data(:,:,1:(nOn+nOff)*ntrials);
    sz = size(data);
    data_tr = reshape(data,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(permute(data_tr,[1 2 4 3])), squeeze(data_f)); 
    data_dfof = permute(bsxfun(@rdivide,data_df, squeeze(data_f)),[1 2 4 3]); 
    clear data_f data_df data_tr

    dirs = unique(dir_mat);
    nStim = length(dirs);
    Stims = [];
    data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, nStim);
    start = 1;
    for iDir = 1:nStim
        ind = find(dir_mat == dirs(iDir));
        Stims = [Stims; dirs(iDir)];
        data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
        start = start +1;
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));
    %data_dfof_avg_all = squeeze(mean(data_dfof_avg(:,:,nOff:nOff+nOn,:),3));

%         img_avg_resp = zeros(1,nStim);
    figure;
    [n, n2] = subplotn(nStim);
    for i = 1:nStim
        subplot(n, n2,i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        colormap(gray)
        title(num2str(Stims(i,:)))
%             img_avg_resp(i) = mean(mean(mean(data_dfof_avg_all(:,:,i),3),2),1);
        %clim([0 max(data_dfof_avg_all(:))./2])
    end
        mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:)]);
        print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_direction.pdf'], '-dpdf','-bestfit')
    movegui('center')

    figure; imagesc(mean(data,3));
    movegui('center')
    axis off
    title([mouse ' ' date])
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_FOV.pdf'], '-dpdf','-bestfit')
