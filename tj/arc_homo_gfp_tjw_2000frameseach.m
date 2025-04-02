clear all;
clear global;
close all;
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


%%

tc = squeeze(mean(mean(data,1),2));

min_tc = min(tc);
max_tc = max(tc);
mean_tc = mean(tc);

new_tc = zeros(1, nframes);
new_tc(find(tc>mean_tc))=1;

data_laser = data(:,:,(new_tc==1));

ind = find(new_tc==1);
ind_diff = find(diff(ind)>2);

 %% Choose register interval - averaging 500 frames, also skipping 2000
t = 500;
nep = floor(size(data_laser,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(data_laser(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end



%% Register data - identify clearest stack and align frames to that
%try a smaller interval here
% data_avg = mean(data(:,:,34001:34500),3); 
data_avg = mean(data_laser(:,:,5501:5600),3); 


if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data_laser(:,:,:),[],[],out);
    %clear out outs
else
    [out, data_reg] = stackRegister(data_laser,data_avg);
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
% figure('Renderer', 'painters', 'Position', [200 200 1200 800])
% subplot(2,2,1)
% imagesc(mean(data_reg(:,:,1:ind_diff(1)),3))
% title('16dir')
% hold on
% subplot(2,2,2)
% imagesc(mean(data_reg(:,:,1+ind_diff(1):ind_diff(2)),3))
% title('Blank')
% subplot(2,2,3)
% imagesc(mean(data_reg(:,:,1+ind_diff(2):ind_diff(3)),3))
% title('1dir')
% subplot(2,2,4)
% imagesc(mean(data_reg(:,:,1+ind_diff(3):length(ind)),3))
% title('Blank')
% colormap("gray")

%5 'blocks'
figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(mean(data_reg(:,:,1:ind_diff(1)),3))
title('1st min')
hold on
subplot(2,3,2)
imagesc(mean(data_reg(:,:,1+ind_diff(1):ind_diff(2)),3))
title('15th min')
subplot(2,3,3)
imagesc(mean(data_reg(:,:,1+ind_diff(2):ind_diff(3)),3))
title('30th min')
subplot(2,3,4)
imagesc(mean(data_reg(:,:,1+ind_diff(3):ind_diff(4)),3))
title('45th min')
subplot(2,3,5)
imagesc(mean(data_reg(:,:,1+ind_diff(4):length(ind)),3))
title('58th min')
colormap("gray")
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_blocks_avg.pdf']),'-dpdf','-bestfit')


%only 100 frames
figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(mean(data_reg(:,:,1:100),3))
title('1st min')
hold on
subplot(2,3,2)
imagesc(mean(data_reg(:,:,1+ind_diff(1):101+ind_diff(1)),3))
title('15th min')
subplot(2,3,3)
imagesc(mean(data_reg(:,:,1+ind_diff(2):101+ind_diff(2)),3))
title('30th min')
subplot(2,3,4)
imagesc(mean(data_reg(:,:,1+ind_diff(3):101+ind_diff(3)),3))
title('45th min')
subplot(2,3,5)
imagesc(mean(data_reg(:,:,1+ind_diff(4):101+ind_diff(4)),3))
title('58th min')
colormap("gray")



%%
% big images
figure;
imagesc(mean(data_reg(:,:,1:100),3))
title('1st min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_1st_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,1+ind_diff(1):101+ind_diff(1)),3))
title('15th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_15th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,1+ind_diff(2):101+ind_diff(2)),3))
title('30th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_30th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,1+ind_diff(3):101+ind_diff(3)),3))
title('45th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_45th_min.pdf']),'-dpdf','-bestfit')

figure;
imagesc(mean(data_reg(:,:,1+ind_diff(4):101+ind_diff(4)),3))
title('60th min')
colormap gray
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_60th_min.pdf']),'-dpdf','-bestfit')

%%

% x = sqrt(out(:,3).^2 + out(:,4).^2);
% figure; plot(x)
% hold on;
% xline([ind_diff(1), ind_diff(2), ind_diff(3)])

x = sqrt(out(:,3).^2 + out(:,4).^2);
figure; plot(x)
hold on;
xline([ind_diff(1), ind_diff(2), ind_diff(3), ind_diff(4)])



%%
t = 250;
nep = floor(size(data_reg(1:ind_diff(1)),2)./t);
[n n2] = subplotn(nep);
a = figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data_reg(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]); 
    colormap("gray")
end

%%
t = 2000;
nep = floor(size(data_reg,3)./t);
[n n2] = subplotn(nep);
a = figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    imagesc(mean(data_reg(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]); 
    colormap("gray")
end
%%

%%
% block1_data = data_reg(:,:,(min(ind):ind_diff(1)));
% block2_data = data_reg(:,:,ind_diff(1)+1:ind_diff(2));
% block3_data = data_reg(:,:,ind_diff(2)+1:ind_diff(3));
% block4_data = data_reg(:,:,ind_diff(3)+1:size(data_reg,3));
% 
% clean_block1_reg = block1_data(:,:,x(min(ind):ind_diff(1))'<=40);
% clean_block2_reg = block2_data(:,:,x(ind_diff(1)+1:ind_diff(2))'<=40);
% clean_block3_reg = block3_data(:,:,x(ind_diff(2)+1:ind_diff(3))'<=20);
% clean_block4_reg = block4_data(:,:,x(ind_diff(3)+1:size(data_reg,3))'<=20);

block1_data = data_reg(:,:,(min(ind):ind_diff(1)));
block2_data = data_reg(:,:,ind_diff(1)+1:ind_diff(2));
block3_data = data_reg(:,:,ind_diff(2)+1:ind_diff(3));
block4_data = data_reg(:,:,ind_diff(3)+1:ind_diff(4));
block5_data = data_reg(:,:,ind_diff(4)+1:size(data_reg,3));

clean_block1_reg = block1_data(:,:,x(min(ind):ind_diff(1))'<=10);
clean_block2_reg = block2_data(:,:,x(ind_diff(1)+1:ind_diff(2))'<=10);
clean_block3_reg = block3_data(:,:,x(ind_diff(2)+1:ind_diff(3))'<=20);
clean_block4_reg = block4_data(:,:,x(ind_diff(3)+1:ind_diff(4))'<=5);
clean_block5_reg = block5_data(:,:,x(ind_diff(4)+1:size(data_reg,3))'<=10);

%%
% figure('Renderer', 'painters', 'Position', [200 200 1200 800])
% subplot(2,2,1)
% imagesc(mean(clean_block1_reg,3))
% title('16dir')
% hold on
% subplot(2,2,2)
% imagesc(mean(clean_block2_reg,3))
% title('Blank')
% subplot(2,2,3)
% imagesc(mean(clean_block3_reg,3))
% title('1dir')
% subplot(2,2,4)
% imagesc(mean(clean_block4_reg,3))
% title('Blank')
% colormap("gray")

figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(mean(clean_block1_reg,3))
title('1st min')
hold on
subplot(2,3,2)
imagesc(mean(clean_block2_reg,3))
title('15th min')
subplot(2,3,3)
imagesc(mean(clean_block3_reg,3))
title('30th min')
subplot(2,3,4)
imagesc(mean(clean_block4_reg,3))
title('45th min')
subplot(2,3,5)
imagesc(mean(clean_block5_reg,3))
title('58th min')
colormap("gray")
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_blocks_clean_avg.pdf']),'-dpdf','-bestfit')


%%
%smaller stack
% figure('Renderer', 'painters', 'Position', [200 200 1200 800])
% subplot(2,2,1)
% imagesc(mean(clean_block1_reg(:,:,1:100),3))
% title('16dir')
% hold on
% subplot(2,2,2)
% imagesc(mean(clean_block2_reg(:,:,1:100),3))
% title('Blank')
% subplot(2,2,3)
% imagesc(mean(clean_block3_reg(:,:,1:100),3))
% title('1dir')
% subplot(2,2,4)
% imagesc(mean(clean_block4_reg(:,:,1:100),3))
% title('Blank')
% colormap("gray")

figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(mean(clean_block1_reg(:,:,1:100),3))
title('1st min')
hold on
subplot(2,3,2)
imagesc(mean(clean_block2_reg(:,:,1:100),3))
title('15th min')
subplot(2,3,3)
imagesc(mean(clean_block3_reg(:,:,1:100),3))
title('30th min')
subplot(2,3,4)
imagesc(mean(clean_block4_reg(:,:,1:100),3))
title('45th min')
subplot(2,3,5)
imagesc(mean(clean_block5_reg(:,:,1:100),3))
title('58th min')
colormap("gray")

%%

cd \\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\231006_i2566\231006_i2566_runs-001

%WRITETIFF Writes array as TIFF file
%   b=writetiff(array,filename, typestr)
%
%   Uses ImageJ code to do write
%   typestr is a matlab type name, supplied to the cast function
%
%   array should be size [nRows,nCols,nFrames,nPlanes] where nPlanes is 1
%   (indexed color) or 3 (RGB color)

% tiff_block1 = writetiff(clean_block1_reg, 'tiff_block1.tiff',)
% 
% imwrite(clean_block1_reg,"my_graphics_file.tif","tif")
% 
% mytiff1 = Tiff('myfile.tif','w');  
% zzz = write(mytiff1,clean_block1_reg);

writetiff(clean_block1_reg, 'block1tiff.tiff', 'double');
writetiff(clean_block2_reg, 'block2tiff.tiff', 'double');
writetiff(clean_block3_reg, 'block3tiff.tiff', 'double');
writetiff(clean_block4_reg, 'block4tiff.tiff', 'double');
writetiff(clean_block5_reg, 'block5tiff.tiff', 'double');


%%

%%
%try to register each block to itself for more clarity***


%%
block1_raw = data_laser(:,:,(min(ind):ind_diff(1)));
block2_raw = data_laser(:,:,ind_diff(1)+1:ind_diff(2));
block3_raw = data_laser(:,:,ind_diff(2)+1:ind_diff(3));
block4_raw = data_laser(:,:,ind_diff(3)+1:size(data_laser,3));

% block1_raw = data_laser(:,:,(min(ind):ind_diff(1)));
% block2_raw = data_laser(:,:,ind_diff(1)+1:ind_diff(2));
% block3_raw = data_laser(:,:,ind_diff(2)+1:ind_diff(3));
% block4_raw = data_laser(:,:,ind_diff(3)+1:ind_diff(4));
% block5_raw = data_laser(:,:,ind_diff(4)+1:size(data_laser,3));

%%
%BLOCK 1
t = 500;
nep = floor(size(block1_raw,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(block1_raw(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end

data_avg_block1_only = mean(block1_raw(:,:,1001:1100),3); 
[out_block1_only, data_reg_block1_only] = stackRegister(block1_raw,data_avg_block1_only);
data_reg_avg_block1_only = mean(data_reg_block1_only,3);

%%
%BLOCK 2
t = 500;
nep = floor(size(block2_raw,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(block2_raw(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end

data_avg_block2_only = mean(block2_raw(:,:,501:600),3); 
[out_block2_only, data_reg_block2_only] = stackRegister(block2_raw,data_avg_block2_only);
data_reg_avg_block2_only = mean(data_reg_block2_only,3);


%%
%BLOCK 3
t = 500;
nep = floor(size(block3_raw,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(block3_raw(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end

data_avg_block3_only = mean(block3_raw(:,:,1001:1100),3); 
[out_block3_only, data_reg_block3_only] = stackRegister(block3_raw,data_avg_block3_only);
data_reg_avg_block3_only = mean(data_reg_block3_only,3);

%%
%BLOCK 4
t = 500;
nep = floor(size(block4_raw,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(block4_raw(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end

data_avg_block4_only = mean(block4_raw(:,:,501:600),3); 
[out_block4_only, data_reg_block4_only] = stackRegister(block4_raw,data_avg_block4_only);
data_reg_avg_block4_only = mean(data_reg_block4_only,3);

%%
%BLOCK 5
t = 500;
nep = floor(size(block5_raw,3)./t);
[n n2] = subplotn(nep);
figure; 
for i = 1:nep; 
    subplot(n,n2,i); 
    % imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):250+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(250+((i-1)*t))]);
    imagesc(mean(block5_raw(:,:,1+((i-1)*t):100+((i-1)*t)),3)); 
    title([num2str(1+((i-1)*t)) '-' num2str(100+((i-1)*t))]);
    % imagesc(mean(data(:,:,1+((i-1)*t):1000+((i-1)*t)),3)); 
    % title([num2str(1+((i-1)*t)) '-' num2str(1000+((i-1)*t))]);
end

data_avg_block5_only = mean(block5_raw(:,:,501:600),3); 
[out_block5_only, data_reg_block5_only] = stackRegister(block5_raw,data_avg_block5_only);
data_reg_avg_block5_only = mean(data_reg_block5_only,3);



%%

% figure('Renderer', 'painters', 'Position', [200 200 1200 800])
% subplot(2,2,1)
% imagesc(data_reg_avg_block1_only)
% title('16dir')
% hold on
% subplot(2,2,2)
% imagesc(data_reg_avg_block2_only)
% title('Blank')
% subplot(2,2,3)
% imagesc(data_reg_avg_block3_only)
% title('1dir')
% subplot(2,2,4)
% imagesc(data_reg_avg_block4_only)
% title('Blank')
% colormap("gray")

figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(data_reg_avg_block1_only)
title('1st min')
hold on
subplot(2,3,2)
imagesc(data_reg_avg_block2_only)
title('15th min')
subplot(2,3,3)
imagesc(data_reg_avg_block3_only)
title('30th min')
subplot(2,3,4)
imagesc(data_reg_avg_block4_only)
title('45th min')
subplot(2,3,5)
imagesc(data_reg_avg_block5_only)
title('58th min')
colormap("gray")
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_blocks_regtoself_avg.pdf']),'-dpdf','-bestfit')

%%
x_block1_only = sqrt(out_block1_only(:,3).^2 + out_block1_only(:,4).^2);
x_block2_only = sqrt(out_block2_only(:,3).^2 + out_block2_only(:,4).^2);
x_block3_only = sqrt(out_block3_only(:,3).^2 + out_block3_only(:,4).^2);
x_block4_only = sqrt(out_block4_only(:,3).^2 + out_block4_only(:,4).^2);
x_block5_only = sqrt(out_block5_only(:,3).^2 + out_block5_only(:,4).^2);


figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
plot(x_block1_only)
title('1st min')
hold on
subplot(2,3,2)
plot(x_block2_only)
title('15th min')
subplot(2,3,3)
plot(x_block3_only)
title('30th min')
subplot(2,3,4)
plot(x_block4_only)
title('45th min')
subplot(2,3,5)
plot(x_block5_only)
title('58th min')
colormap("gray")

%%
clean_data_reg_block1_only = data_reg_block1_only(:,:,x_block1_only'<=10);
clean_data_reg_block2_only = data_reg_block2_only(:,:,x_block2_only'<=10);
clean_data_reg_block3_only = data_reg_block3_only(:,:,x_block3_only'<=10);
clean_data_reg_block4_only = data_reg_block4_only(:,:,x_block4_only'<=10);
clean_data_reg_block5_only = data_reg_block5_only(:,:,x_block5_only'<=10);


%%
figure('Renderer', 'painters', 'Position', [200 200 1200 800])
subplot(2,3,1)
imagesc(mean(clean_data_reg_block1_only,3))
title('1st min')
hold on
subplot(2,3,2)
imagesc(mean(clean_data_reg_block2_only,3))
title('15th min')
subplot(2,3,3)
imagesc(mean(clean_data_reg_block3_only,3))
title('30th min')
subplot(2,3,4)
imagesc(mean(clean_data_reg_block4_only,3))
title('45th min')
subplot(2,3,5)
imagesc(mean(clean_data_reg_block5_only,3))
title('58th min')
colormap("gray")
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_blocks_regtoself_clean_avg.pdf']),'-dpdf','-bestfit')

%%


%%

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


