clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirType2_ExptList';
eval(ds)
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\CrossOri\RandDirType2';
nexp = length(expt);
dirs = deg2rad(0:22.5:360-22.5);
angle = 45;
area = 'PM';

ind_angle  = find([expt.angle]==angle);
ind_area  = find(strcmp([expt.img_loc],area));
ind = intersect(ind_angle,ind_area);
avg_resp_dir_all = [];
offset = 0;
resp_ind_all = [];
for iexp = ind
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    saveLoc = expt(iexp).saveLoc;
    
    LG_base = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' saveLoc];
    
    fprintf([mouse ' ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    avg_resp_dir_all = cat(1, avg_resp_dir_all, squeeze(avg_resp_dir));
    length(resp_ind)
    resp_ind_all = [resp_ind_all; resp_ind+offset];
    offset = offset+size(avg_resp_dir,1);
end

nDir = size(avg_resp_dir_all,2);
dirs = deg2rad(0:22.5:360-22.5);
nCells = size(avg_resp_dir_all,1);
for iC = 1:nCells
    [b_g(iC),k1_g(iC),R1_g(iC),R2_g(iC),u1_g(iC),u2_g(iC),sse_g(iC),R_square_g(iC)] = miaovonmisesfit_dir(dirs,avg_resp_dir_all(iC,:,1,1));
    [b_p(iC),k1_p(iC),R1_p(iC),R2_p(iC),u1_p(iC),u2_p(iC),sse_p(iC),R_square_p(iC)] = miaovonmisesfit_dir(dirs,avg_resp_dir_all(iC,:,2,1));
end

DSI_g = (R1_g-R2_g)./(R1_g+R2_g);
DSI_p = (R1_p-R2_p)./(R1_p+R2_p);

indUse = intersect(intersect(find(DSI_g>0.3),find(DSI_p>0.3)),resp_ind_all);

prefDiff = u1_g-u1_p;
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),4,deg2rad(angle),1,1);

save(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_' area '_Summary.mat']), 'avg_resp_dir_all','resp_ind_all','prefDiff','DSI_g','DSI_p','indUse','ioc_ang','av_ang','b_g','k1_g','R1_g','R2_g','u1_g','u2_g','sse_g','R_square_g','b_p','k1_p','R1_p','R2_p','u1_p','u2_p','sse_p','R_square_p')

%%
areaList = {'V1','PM'};
angle = 45;
figure
for iarea = 1:length(areaList)
    load(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_' areaList{iarea} '_Summary.mat']))
    subplot(3,2,1+(iarea-1))
    scatter(rad2deg(u1_g(indUse)),rad2deg(u1_p(indUse)))
    refline(1)
    refline(1,-angle)
    refline(1,-rad2deg(ioc_ang))
    refline(1,-rad2deg(av_ang))
    xlabel('Pref dir- Grating')
    ylabel('Pref dir- Plaid')
    xlim([0 360])
    ylim([0 360])
    axis square
    title([areaList{iarea} '- ' num2str(length(indUse)) ' cells'])
    %legend('data','fast comp','slow comp','IOC','VA','Location','northeastoutside')
    
    subplot(3,2,3+(iarea-1))
    scale = length(indUse)/5;
    polarhistogram(prefDiff(indUse),32);
    hold on
    polarplot([0;0]*pi/180,[0;1]*scale)
    polarplot([angle;angle]*pi/180,[0;1]*scale)
    polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
    polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)

    diff_thresh = 7.5;
    prefDiff_deg = rad2deg(prefDiff);
    ind_fast = intersect(indUse, intersect(find(prefDiff_deg>-diff_thresh), find(prefDiff_deg<diff_thresh)));
    ind_slow = intersect(indUse,  intersect(find(prefDiff_deg>(angle-diff_thresh)), find(prefDiff_deg<(angle+diff_thresh))));
    ind_ioc = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(ioc_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(ioc_ang)+diff_thresh))));
    ind_av = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(av_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(av_ang)+diff_thresh))));
    ind_fast_av = intersect(ind_fast,ind_av);
    ind_fast = setdiff(ind_fast, ind_fast_av);
    ind_av = setdiff(ind_av, ind_fast_av);

    subplot(3,2,5+(iarea-1))
    piechart([length(ind_fast) length(ind_slow) length(ind_ioc) length(ind_av) length(ind_fast_av)],{'Fast','Slow','IOC','VA','Fast+AV'})
    title([num2str(angle) ' deg'])
    
    sgtitle(['Plaid angle: ' num2str(angle)])
end
print(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) '_deg_' cell2mat(areaList) '_Summary.pdf']), '-fillpage', '-dpdf')

%% stimuli
[s2d_test, p2d] = sin2d2(-15:0.1:15,10,0,90);
[s2d_mask, p2d] = sin2d2(-15:0.1:15,10,deg2rad(0),135);
plaid = s2d_test+s2d_mask;
subplot(2,2,1)
imagesc(s2d_test)
clim([-2 2])
colormap gray
axis square
axis off
subplot(2,2,2)
imagesc(s2d_mask)
clim([-2 2])
colormap gray
axis square
axis off
subplot(2,2,3)
imagesc(plaid)
clim([-2 2])
colormap gray
axis square
axis off
print(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) '_stimuli.pdf']), '-bestfit', '-dpdf')


%%
[prefResp prefDir] = max(avg_resp_dir_all(:,:,1,1),[],2);
nCells = size(avg_resp_dir_all,1);
avg_resp_dir_all_shift = zeros(size(avg_resp_dir_all));
for i = 1:nCells
    avg_resp_dir_all_shift(i,:,:,:) = circshift(avg_resp_dir_all(i,:,:,:),1-prefDir(i),2);
end
avg_resp_dir_all_shift_circ = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:));
figure(2);
subplot(2,2,1)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,2,1),1))
title(['Fast comp- n =' num2str(length(ind_fast))])
rlim([0 0.2])
subplot(2,2,2)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,2,1),1))
title(['Slow comp- n =' num2str(length(ind_slow))])
rlim([0 0.2])
subplot(2,2,3)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,2,1),1))
title(['IOC- n =' num2str(length(ind_ioc))])
rlim([0 0.2])
subplot(2,2,4)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,2,1),1))
title(['VA- n =' num2str(length(ind_av))])
rlim([0 0.2])
sgtitle([num2str(angle) ' deg'])
print(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_' area '_TuningCurves.pdf']), '-fillpage', '-dpdf')
%%

angle = 135;

if ~exist(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_Summary.mat']))

    ind  = find([expt.angle]==angle);
    
    avg_resp_dir_all = [];
    offset = 0;
    resp_ind_all = [];
    for iexp = ind
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc;
        ImgFolder = expt(iexp).coFolder;
        time = expt(iexp).coTime;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        
        %LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';
        
        fsaveasf([mouse ' ' date '\n'])
    
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        avg_resp_dir_all = cat(1, avg_resp_dir_all, squeeze(avg_resp_dir));
        resp_ind_all = [resp_ind_all; resp_ind+offset];
        offset = offset+size(avg_resp_dir,1);
    end
    
    nDir = size(avg_resp_dir_all,2);
    dirs = deg2rad(0:22.5:360-22.5);
    nCells = size(avg_resp_dir_all,1);
    for iC = 1:nCells
        [b_g(iC),k1_g(iC),R1_g(iC),R2_g(iC),u1_g(iC),u2_g(iC),sse_g(iC),R_square_g(iC)] = miaovonmisesfit_dir(dirs,avg_resp_dir_all(iC,:,1,1));
        [b_p(iC),k1_p(iC),R1_p(iC),R2_p(iC),u1_p(iC),u2_p(iC),sse_p(iC),R_square_p(iC)] = miaovonmisesfit_dir(dirs,avg_resp_dir_all(iC,:,2,1));
    end
    
    DSI_g = (R1_g-R2_g)./(R1_g+R2_g);
    DSI_p = (R1_p-R2_p)./(R1_p+R2_p);
    
    indUse = intersect(intersect(find(DSI_g>0.5),find(DSI_p>0.5)),resp_ind_all);
    
    prefDiff = u1_g-u1_p;
    [ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),4,deg2rad(angle),1);
    
    save(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_Summary.mat']), 'avg_resp_dir_all','resp_ind_all','prefDiff','DSI_g','DSI_p','indUse','ioc_ang','av_ang','b_g','k1_g','R1_g','R2_g','u1_g','u2_g','sse_g','R_square_g','b_p','k1_p','R1_p','R2_p','u1_p','u2_p','sse_p','R_square_p')
else
    load(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_Summary.mat']))
end

figure(1);
subplot(2,2,3)
scatter(rad2deg(u1_g(indUse)),rad2deg(u1_p(indUse)))
refline(1)
refline(1,-angle)
refline(1,-rad2deg(ioc_ang))
refline(1,-rad2deg(av_ang))
xlabel('Pref dir- Grating')
ylabel('Pref dir- Plaid')
xlim([0 360])
ylim([0 360])
axis square
title(['Plaid angle: ' num2str(angle)])
legend('data','fast comp','slow comp','IOC','VA','Location','northeastoutside')


subplot(2,2,4)
polarhistogram(prefDiff(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*60)
polarplot([angle;angle]*pi/180,[0;1]*60)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*60)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*60)

print(fullfile(fnout,'crossOriRandDir_Type2_Summary.pdf'), '-fillpage', '-dpdf')

diff_thresh = 7.5;
prefDiff_deg = rad2deg(prefDiff);
ind_fast = intersect(indUse, intersect(find(prefDiff_deg>-diff_thresh), find(prefDiff_deg<diff_thresh)));
ind_slow = intersect(indUse,  intersect(find(prefDiff_deg>(angle-diff_thresh)), find(prefDiff_deg<(angle+diff_thresh))));
ind_ioc = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(ioc_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(ioc_ang)+diff_thresh))));
ind_av = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(av_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(av_ang)+diff_thresh))));
ind_fast_av = intersect(ind_fast,ind_av);
ind_fast = setdiff(ind_fast, ind_fast_av);
ind_av = setdiff(ind_av, ind_fast_av);

[prefResp prefDir] = max(avg_resp_dir_all(:,:,1,1),[],2);
nCells = size(avg_resp_dir_all,1);
avg_resp_dir_all_shift = zeros(size(avg_resp_dir_all));
for i = 1:nCells
    avg_resp_dir_all_shift(i,:,:,:) = circshift(avg_resp_dir_all(i,:,:,:),1-prefDir(i),2);
end
avg_resp_dir_all_shift_circ = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:));

figure(4);
subplot(2,2,1)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_fast,:,2,1),1))
title(['Fast comp- n =' num2str(length(ind_fast))])
rlim([0 0.2])
subplot(2,2,2)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_slow,:,2,1),1))
title(['Slow comp- n =' num2str(length(ind_slow))])
rlim([0 0.2])
subplot(2,2,3)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_ioc,:,2,1),1))
title(['IOC- n =' num2str(length(ind_ioc))])
rlim([0 0.2])
subplot(2,2,4)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_av,:,2,1),1))
title(['VA- n =' num2str(length(ind_av))])
rlim([0 0.2])
sgtitle([num2str(angle) ' deg'])
print(fullfile(fnout,['crossOriRandDir_Type2_' num2str(angle) 'deg_TuningCurves.pdf']), '-fillpage', '-dpdf')

figure(3)
subplot(1,2,2)
piechart([length(ind_fast) length(ind_slow) length(ind_ioc) length(ind_av) length(ind_fast_av)],{'Fast','Slow','IOC','VA','Fast+AV'})
title([num2str(angle) ' deg'])
print(fullfile(fnout,['crossOriRandDir_Type2_FractDir.pdf']), '-fillpage', '-dpdf')

