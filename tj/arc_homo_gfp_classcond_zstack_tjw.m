clear all;
clear global;
close all;
clc;

%% get path names
ref_date = '240118';
date = '240118';
time = strvcat('1037');
alignToRef = 1;
ImgFolder = strvcat('003');
mouse = 'i2570';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
ref_str = 'runs-002';
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
    nframes = max(info.frame);
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

%%

data_avg = squeeze(mean(reshape(data,[529 796 100 11]),3));
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_avg(:,:,i))
end

%%
data_by_z = reshape(data,[529 796 100 11]);

% figure;
% for i = 1:100
%     subplot(10,10,i)
%     imagesc(data_by_z(:,:,i,1))
% end



%%
data_reg = [];
outs = [];

for i = 1:11;
    [out, reg] = stackRegister(data_by_z(:,:,:,i), data_avg(:,:,i));
    data_reg = cat(3,data_reg,reg);
end

%%
data_reg_avg = squeeze(mean(reshape(data_reg,[529 796 100 11]),3));
figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_reg_avg(:,:,i))
end

%%

data_reg_reg = [];

for i = 1:11;
    if i == 1
        data_reg_reg = data_reg_avg(:,:,1);
    else
        [out, reg] = stackRegister(data_reg_avg(:,:,i), data_reg_avg(:,:,i-1));
        data_reg_reg = cat(3,data_reg_reg,reg);
    end
end

figure;
for i = 1:9
    subplot(3,3,i)
    imagesc(data_reg_reg(:,:,i))
end

%%
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
CD = fullfile(fnout, [date '_' mouse '\' date '_' mouse '_' 'runs-' ImgFolder(irun,:)]);
cd(CD);
writetiff(data_reg_reg, 'block2.tiff', 'double');


%%
x = sqrt(out(:,3).^2 + out(:,4).^2);
figure; plot(x)
%%

 %% Choose register interval - averaging 500 frames, also skipping 2000
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end


%% Register data - identify clearest stack and align frames to that
%try a smaller interval here
data_avg = mean(data(:,:,32001:32250),3); 

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
    %clear out outs
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
figure; 
imagesq(data_reg_avg); 
colormap('gray')
truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')

%%
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
data_size = round(size(data_reg,3)/4);
data_1 = 1:data_size;
data_2 = data_size+1:data_size*2;
data_3 = data_size*2+1:data_size*3;
data_4 = data_size*3+1:data_size*4;
%%
% big images
figure;
imagesc(mean(data_reg(:,:,1:100),3))
title('1st min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_1st_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,data_2:100+data_2),3))
title('15th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_15th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,data_3:100+data_3),3))
title('30th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_30th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,data_4:100+data_4),3))
title('45th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_45th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,57500:57600),3))
title('60th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_60th_min.pdf']),'-dpdf','-bestfit')

%%

x = sqrt(out(:,3).^2 + out(:,4).^2);
figure; plot(x)

%%

clean_data_reg = data_reg(:,:,x <= 30);
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



%%

CD = fullfile(fnout, [date '_' mouse '\' date '_' mouse '_' 'runs-' ImgFolder(irun,:)]);
cd(CD);

clean_size = round(size(clean_data_reg,3)/4)

%%
data_1 = clean_data_reg(:,:,1:clean_size);
data_2 = clean_data_reg(:,:,clean_size+1:clean_size*2);
data_3 = clean_data_reg(:,:,clean_size*2+1:clean_size*3);
data_4 = clean_data_reg(:,:,clean_size*3+1:clean_size*4);


writetiff(data_1, 'cleandata1.tiff', 'double');
writetiff(data_2, 'cleandata2.tiff', 'double');
writetiff(data_3, 'cleandata3.tiff', 'double');
writetiff(data_4, 'cleandata4.tiff', 'double');


%%
figure;
imagesc(img(:,:,2))
colormap("gray")


%%

irun = 1;
WL = '920';
ImgFolder = strvcat('002');
run = catRunName(ImgFolder, nrun);
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);

data_920_green = squeeze(data_temp(1,:,:,:)); %PMT 0 (green)
% data_920_red = squeeze(data_temp(2,:,:,:)); %PMT 1 (red)

[out_920_green_regtoself data_920_green_regtoself] = stackRegister(data_920_green,mean(data_920_green,3)); 
% [out_920_red_regtoself data_920_red_regtoself] = stackRegister(data_920_red,mean(data_920_red,3));

%%

irun = 1;
WL = '1040';
ImgFolder = strvcat('003');
run = catRunName(ImgFolder, nrun);
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(tj_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);

data_1040_green = squeeze(data_temp(1,:,:,:)); %PMT 0 (green)
% data_1040_red = squeeze(data_temp(2,:,:,:)); %PMT 1 (red)

[out_1040_green_regtoself data_1040_green_regtoself] = stackRegister(data_1040_green,mean(data_1040_green,3)); 
% [out_1040_red_regtoself data_1040_red_regtoself] = stackRegister(data_1040_red,mean(data_1040_red,3));

%%
figure; 
subplot(2,2,1)
imagesc(mean(data_920_green_regtoself,3))
subplot(2,2,2)
imagesc(mean(data_920_red_regtoself,3))
subplot(2,2,3)
imagesc(mean(data_1040_green_regtoself,3))
subplot(2,2,4)
imagesc(mean(data_1040_green_regtoself,3))

%%
figure; 
subplot(2,2,1)
imagesc(mean(data_920_green_regtoself(:,:,1:250),3))
subplot(2,2,2)
imagesc(mean(data_920_green_regtoself(:,:,251:500),3))
subplot(2,2,3)
imagesc(mean(data_920_green_regtoself(:,:,501:750),3))
subplot(2,2,4)
imagesc(mean(data_920_green_regtoself(:,:,751:1000),3))
colormap('gray')

%%
figure;
imagesc(mean(data_920_green_regtoself,3))
colormap('gray')

figure;
imagesc(mean(data_920_red_regtoself,3))
colormap('gray')
%%
CD = fullfile(fnout, [date '_' mouse '\' date '_' mouse '_runs-001']);
cd(CD);
writetiff(data_920_green_regtoself, 'green_snap_reg.tiff', 'double');
