%% Parameters
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 150./down;
nOFF = 150./down;
pre_win = [floor(0.75*nOFF) nOFF];
post_win = [nOFF+1 nOFF+nON];
nCond = 3;
% TF = [1 2 4 8 16];
% SF = [0.32 0.16 0.08 0.04 0.02];
Az = [0 15 30];
El = [0];

date = '140209';
mouse = 'G008';

%% load data
data = readrawfile;

%if two-channel, choose just green channel
data = data(:,:,:,1);

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down
%% register
data_avg = mean(data_sub(:,:,200:210),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% create resp matrix
%if you need to resort by trial type for randomization, do it first.
%code assumes same number of reps per stimulus
siz = size(data_reg);
nRep = size(data_reg,3)./((nON+nOFF)*nCond);
roi_stim = zeros(siz(1), siz(2), nOFF+nON,nRep);
start = 1;
rep = 1;
for iCond = 1:nCond; 
    for iRep = 1:nRep;
        roi_stim(:,:,:,rep) = data_reg(:,:,start:start-1+nON+nOFF);
        start = start+nON+nOFF;
        rep = rep+1;
    end
end
resp_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
resp_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));

%% pixel based ttest

alpha = .05./(nCond);
Info_ttest_mat = zeros(siz(1),siz(2),nCond);
b= 5;

f1 = fspecial('average');
siz = size(resp_on);
resp_on_long = reshape(resp_on, [siz(1) siz(2)*siz(3)]);
resp_off_long = reshape(resp_off, [siz(1) siz(2)*siz(3)]);

resp_on_sm = reshape(filter2(f1,resp_on_long),[siz(1) siz(2) siz(3)]);
resp_off_sm = reshape(filter2(f1,resp_off_long),[siz(1) siz(2) siz(3)]);

clear resp_on
clear resp_off
clear resp_on_long
clear resp_off_long

for iy = b+1:siz(1)-b
    fprintf([num2str(iy) ' '])
    for ix = b+1:siz(2)-b
        start = 1;
        p_ttestB = zeros(1,1,nCond);
        for iCond = 1:nCond
            [h_ttestB1,p_ttestB1] = ttest(resp_off_sm(iy,ix,start:start-1+nRep),resp_on_sm(iy,ix,start:start-1+nRep),alpha,'left');
            p_ttestB(1,1,iCond) = p_ttestB1;
            start = start+nRep;
        end
    Info_ttest_mat(iy,ix,:) = p_ttestB;
    end
end

clear resp_on_sm
clear resp_off_sm

siz = size(Info_ttest_mat);
Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]);
ttest_smooth = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) siz(3)]);
ttest_mask = min(ttest_smooth,[],3) < alpha;

ttest_mask(1:5,1:end) = 0;
ttest_mask(1:end, 1:5) = 0;
ttest_mask(1:end, 251:end) = 0;
ttest_mask(235:end,1:end) = 0;

figure; imagesq(ttest_mask);
        
%% create dF/F movie to find active boutons
nRep = size(data_reg,3)./((nON+nOFF)*nCond);

%find off and on frames
base_ind = zeros(1,(length(pre_win(1):pre_win(2))*nCond*nRep));
start = 1;
for iCond = 1:(nRep*nCond)
    base_ind(1, start:start+length(pre_win(1):pre_win(2))-1) = pre_win(1)+((iCond-1)*(nON+nOFF)):pre_win(2) + ((iCond-1)*(nOFF+nON));
    start = start+length(pre_win(1):pre_win(2));
end

base_avg = mean(data_reg(:,:,base_ind),3);

%dF/F
dF_data = bsxfun(@minus,data_reg, base_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, base_avg);
max_dF = max(dFoverF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS- local maxima

max_interp = interp2(max_dF);
f1 = fspecial('average');   
max_interp_sm = filter2(f1, max_interp);
siz2 = size(max_interp);
Xi = 1:2:siz2(1);
Yi = 1:2:siz2(2);
stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);

%searches 7x7 for brightest - could  be 5x5 or 3x3
siz = size(max_dF);
local_max = zeros(siz(1), siz(2));
border = 5;
for iy = border:(siz(1)-border);
    for ix = border:(siz(2)-border);            
        sub = stack_max_interp_sm(iy-3:iy+3,ix-3:ix+3);
        sub_long = reshape(sub, [1 49]);
        [sub_long_order ind_order] = sort(sub_long);
        if ind_order(end)==25
            local_max(iy,ix) = 1;
        end
    end
end
local_max_long = reshape(local_max, [siz(1)*siz(2) 1]);
ind_local_max = find(local_max_long==1);
figure; imagesq(local_max);

siz = size(Info_ttest_mat);
Info_ttest_mat_long = interp2(reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]));
Info_ttest_smooth = filter2(f1,Info_ttest_mat_long);
siz_interp = size(Info_ttest_smooth);
Xi = 1:2:siz_interp(1);
Yi = 1:2:siz_interp(2);
ttest_smooth_siz = interp2(Info_ttest_smooth, Yi', Xi);
ttest_smooth = min(reshape(ttest_smooth_siz, [siz(1) siz(2) siz(3)]),[],3);
ttest_long = reshape(ttest_smooth, [siz(1)*siz(2) 1]);
alpha = 0.05/nCond;
ind_highP = find(ttest_long(ind_local_max,:)>=(alpha));
local_max_long(ind_local_max(ind_highP,:),:) = 0;

local_max = reshape(local_max_long, [siz(1) siz(2)]);
n_pix = sum(sum(local_max));
[i, j] = find(local_max ==1);

outDir = '\\crash.dhe.duke.edu\ashley';
fn_local = fullfile(outDir,'analysis',[date '_' mouse '_local_max.mat']);
        save(fn_local, 'local_max', 'n_pix', 'i', 'j'); 
        
%expand local maxima (goes to 5x5 -  could be 3x3)
FOV = zeros(size(local_max));
for ipix = 1:n_pix
    FOV(i(ipix)-2:i(ipix)+2,j(ipix)-2:j(ipix)+2) = 1;
end

figure; imagesq(FOV)

%timecourses
data_TC = zeros(size(dFoverF_data,3),n_pix);
for ipix = 1:n_pix
    data_TC(:,ipix) = mean(mean(dFoverF_data(i(ipix)-2:i(ipix)+2,j(ipix)-2:j(ipix)+2,:),1),2);
end
figure; tcOffsetPlot(data_TC)

%% plot responses
stim_mat = zeros(nCond,nRep,length(pre_win(1):post_win(2)));
for iRep = 1:nRep
    for iCond = 1:nCond 
        stim_mat(iCond,iRep,:) = pre_win(1)+((iCond-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nCond)): nOFF+nON +((iCond-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nCond));
    end
end

%for retinotopy
for iCell = 7:15;
    figure;
    for iCond = 1:nCond
        subplot(1,nCond,iCond)
        rep_mat = zeros(length(pre_win(1):post_win(2)),nRep);
        for iRep = 1:nRep
            plot(pre_win(1)-post_win(1):post_win(2)-post_win(1),data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell);
        end
        plot(pre_win(1)-post_win(1):post_win(2)-post_win(1), mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_Az = rem(iCond,size(Az,2));
        if stim_Az == 0
            stim_Az = size(Az,2);
        end
        stim_El= ceil(iCond./size(Az,2));
        title(['Az = ' num2str(Az(1,stim_Az)) ' El = ' num2str(El(1,stim_El))]);
    end
end

%for speed
for iCell = 7:15;
    figure;
    for iCond = 1:nCond
        subplot(1,nCond,iCond)
        rep_mat = zeros(length(pre_win(1):post_win(2)),nRep);
        for iRep = 1:nRep
            plot(pre_win(1)-post_win(1):post_win(2)-post_win(1),data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iCond,iRep,:))',iCell);
        end
        plot(pre_win(1)-post_win(1):post_win(2)-post_win(1), mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_TF = TF(iCond);
        stim_SF= SF(iCond);
        stim_speed = stim_TF./stim_SF;
        title(['Speed = ' num2str(stim_speed)]);
    end
end
