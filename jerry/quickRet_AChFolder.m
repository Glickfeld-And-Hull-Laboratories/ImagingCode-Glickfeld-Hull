
close all

date = '250627';
mouse = 'i2207';
ImgFolder = '002';
time = '1008';
doReg = 1;
nrun = size(ImgFolder,1);
rc = behavConstsAV;
subnum = mouse;
datemouse = [date '_' mouse]
saveRet = 1;
experimentFolder = 'SST_YM90K';

isilon_home = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\';
analysis_root = [isilon_home 'ACh\Analysis\2p_analysis'];

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

data = [];
clear temp
for irun = 1:nrun
    CD = fullfile('G:\home\ACh\Data\2p_data',mouse,date,ImgFolder(irun,:));
    cd(CD)
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
   
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

    if size(data_temp,1) > 1
        data_temp = data_temp(1,:,:,:);
    end
    
     fName = ['G:\Behavior\Data\data-' subnum '-' date '-' time(irun,:) '.mat'];

    load(fName);
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

% %% plot F std over time
% frame_vector = nan(size(data_temp,3),1);
% 
% for frame = 1:size(data_temp,3)
%     this_frame = data_temp(:,:,frame);
%     frame_pix = reshape(this_frame,[],1);
%     frame_mean = mean(double(frame_pix));
%     frame_vector(frame) = frame_mean;
% end
% 
% std_session = std(frame_vector);
% 
% figure;
% plot(frame_vector);
% sgtitle(['i3328 Retinotopy Mean F; STD = ' num2str(std_session)]);
% ylabel("Average F")
% xlabel('Frame #')
% 
% % clear frame_pix
% %% resume retinotopy






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
    [out data_reg] = stackRegGPU(data,data_avg);
    data = data_reg;
    clear data_reg
    end
    
    sz = size(data);
    ntrials = 45;
    data = data(:,:,1:(nOn+nOff)*ntrials);
    if size(data,3) < 100000
       
        Az = celleqel2mat_padded(expt_input.tGratingAzimuthDeg);
        El = celleqel2mat_padded(expt_input.tGratingElevationDeg);
        if (nOn+nOff)*ntrials > size(data,3)
            ntrials = floor(size(data,3)/((nOn+nOff)));    
            Az = Az(1:ntrials);
            El = El(1:ntrials);
        end
        data = data(:,:,1:(nOn+nOff)*ntrials);
        sz = size(data);
        data_tr = reshape(data,[sz(1), sz(2), nOn+nOff, ntrials]);
        data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
        data_df = bsxfun(@minus, double(permute(data_tr,[1 2 4 3])), squeeze(data_f)); 
        data_dfof = permute(bsxfun(@rdivide,data_df, squeeze(data_f)),[1 2 4 3]); 
        clear data_f data_df data_tr

        Azs = unique(Az);
        Els = unique(El);
        if min(Els,[],2)<0
            Els = fliplr(Els);
        end
        nStim = length(Azs).*length(Els);
        Stims = [];
        data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
        start = 1;
        for iEl = 1:length(Els)
            ind1 = find(El == Els(iEl));
            for iAz = 1:length(Azs)
                Stims = [Stims; Els(iEl) Azs(iAz)];
                ind2 = find(Az == Azs(iAz));
                ind = intersect(ind1,ind2);
                data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
                start = start +1;
            end
        end
        clear data_dfof
        myfilter = fspecial('gaussian',[20 20], 0.5);
        %data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:75,:),myfilter),3));
         data_dfof_avg_all = imgaussfilt(squeeze(mean(data_dfof_avg(:,:,nOff:nOff+nOn,:),3)),1.5);
        
%         img_avg_resp = zeros(1,nStim);
        figure; 
        for i = 1:nStim
            subplot(length(Els),length(Azs),i); 
            imagesc(data_dfof_avg_all(:,:,i)); 
            colormap(gray)
            title(num2str(Stims(i,:)))
%             img_avg_resp(i) = mean(mean(mean(data_dfof_avg_all(:,:,i),3),2),1);
            %clim([0 max(data_dfof_avg_all(:))./2])
        end
        fnout = fullfile(analysis_root,experimentFolder,mouse,date,ImgFolder);
        mkdir(fnout);
        if saveRet == 1
            print(fullfile(fnout,'Retinotopy.pdf'),'-dpdf','-bestfit');
        end
        
%         if strcmp(rc.name,'linds')
%             mkdir(['\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:)]);
%             print(['\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_retinotopy.pdf'], '-dpdf','-bestfit')
%         elseif strcmp(rc.name,'robin')
%             mkdir(['R:\home\robin\Imaging\Analysis\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:)]);
%             print(['R:\home\robin\Imaging\Analysis\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_retinotopy.pdf'], '-dpdf','-bestfit')    
%         end
% 
%         if strcmp(rc.name,'ashle')
%             if ~exist(fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',date,ImgFolder))
%                 mkdir(fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',date,ImgFolder))
%             end
%             print(fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',date,ImgFolder,'quickRet.pdf'),'-dpdf','-fillpage')
%         end
%         
        
        pixThreshold = 0.2*max(data_dfof_avg_all(:));
        img_avg_resp = zeros(1,nStim);
        for i = 1:nStim 
            img = data_dfof_avg_all(:,:,i);
            img_thresh = zeros(size(img));
            img_thresh(img>pixThreshold) = img(img>pixThreshold);
            img_avg_resp(i) = mean(img(img(:)>pixThreshold));
        end
        figure
        heatmap = imagesc(fliplr(rot90(reshape(img_avg_resp,length(Els),length(Azs)),3)));
        heatmap.Parent.YTick = 1:length(Els);
        heatmap.Parent.YTickLabel = strread(num2str(Els),'%s');    
        heatmap.Parent.XTick = 1:length(Azs);
        heatmap.Parent.XTickLabel = strread(num2str(Azs),'%s');
        xlabel('Azimuth');
        ylabel('Elevation');
        colorbar
%         caxis([0 pixThreshold/.4])
%         caxis([-0.1 0.1])
    
    else
        data_tc = squeeze(mean(mean(data,1),2));
        data_tr = reshape(data_tc,[nOn+nOff, ntrials]);
        data_f = mean(data_tr(nOff/2:nOff,:),1);
        data_df = bsxfun(@minus, data_tr, data_f); 
        data_dfof = bsxfun(@rdivide,data_df(nOff:nOff+nOn,:), squeeze(data_f)); 

        Az = celleqel2mat_padded(expt_input.tGratingAzimuthDeg);
        El = celleqel2mat_padded(expt_input.tGratingElevationDeg);
        Azs = unique(Az);
        Els = unique(El);
        if min(Els,[],2)<0
            Els = fliplr(Els);
        end
        for i=1:length(Els)
            for j = 1:length(Azs)
                ind = (El==Els(i)).*(Az==Azs(j));
                respMap(i, j) = mean(data_dfof(find(ind)));
            end
        end
        figure;
        ax=gca;
        imagesc(respMap);
        ax.XTickLabel = Azs;
        ax.YTickLabel = Els;
        title('Total image average dF/F response map')
        xlabel('Azimuth (deg)')
        ylabel('Elevation (deg)')

    end
    
    figure; imagesc(mean(data,3));
    axis off
    title([mouse ' ' date])

    if isfield(info,'frame')
        fprintf('Photodiode data collected')
    else
        fprintf('WARNING- NOT PHOTODIODE DATA!!!!!')
    end
%     if strcmp(rc.name,'linds')
%     print(['\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_FOV.pdf'], '-dpdf','-bestfit')
%     elseif strcmp(rc.name,'robin')
%            print(['R:\home\robin\Imaging\Analysis\' date '_' mouse '\' date '_' mouse '_' ImgFolder(irun,:) '\' date '_' mouse '_' ImgFolder(irun,:) '_FOV.pdf'], '-dpdf','-bestfit')
%     end