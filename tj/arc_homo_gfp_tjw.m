clear all;
clear global;
clc;

%% get path names
ref_date = '230924';
date = '230924';
time = strvcat('1023');
alignToRef = 1;
ImgFolder = strvcat('003');
mouse = 'i2566';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-003';
run_str = catRunName(ImgFolder, nrun);
tj_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %MAKE SURE TO SET THIS
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';

%% load and register - same as d1 code
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse  '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            nframes = ntrials*(nOn+nOff);
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

 %% Choose register interval - averaging 500 frames, also skipping 2000
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(data(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end


%% Register data - identify clearest stack and align frames to that
%try a smaller interval here
% data_avg = mean(data(:,:,34001:34500),3); 
data_avg = mean(data(:,:,34001:34100),3); 


if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg,3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%data_avg = selected stack to register to
%data_reg = all frames registered
%data_reg_avg = mean of all registered frames

%% test stability - see how the image looks averaged across all frames now
figure; imagesq(data_reg_avg); truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')

%%
t = 2500;
nep = floor(size(data_reg,3)./t);
[n n2] = subplotn(nep);
a = figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data_reg(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); 
    colormap("gray")
end
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_every500.pdf']),'-dpdf','-bestfit')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_every500.mat']), 'a')

%%
t = 14400;
nep = floor(size(data_reg,3)./t);
[n n2] = subplotn(nep);
a = figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data_reg(:,:,1+((i-1)*t):14400+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(14400+((i-1)*t))]); 
    colormap("gray")
end
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_every14400.pdf']),'-dpdf','-bestfit')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_every14400.mat']), 'a')

%%

x = sqrt(out(:,3).^2 + out(:,4).^2);
figure; plot(x)

clean_data_reg = data_reg(:,:,x <= 10);
fprintf(['Removed ' num2str(size(data_reg,3)-size(clean_data_reg,3)) ' frames \r\n'])

t = 2500;
nep = floor(size(clean_data_reg,3)./t);
[n n2] = subplotn(nep);
a = figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(clean_data_reg(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); 
    colormap("gray")
end


%%
figure; 
if nep<5
    [n n2] = subplotn(nep);
else
    [n, n2] = subplotn(4);
end
start = 1;
x = 0;
for i = 1:nep
    if start > 4 %?***
        %suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        % print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stacksof500' num2str(start-4) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,i-(x.*4))
    imagesc(mean(clean_data_reg(:,:,1+((i-1)*t):500+((i-1)*t)),3));
    title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(clean_data_reg(:,:,1+((i-1)*t):100+((i-1)*t)),3));
    % title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    colormap("gray")
    start = start+1;
end
% suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
% print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stacksof500' num2str(start-4) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')


