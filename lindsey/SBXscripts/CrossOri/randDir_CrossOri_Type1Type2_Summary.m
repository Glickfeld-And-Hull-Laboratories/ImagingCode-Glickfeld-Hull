clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirType2_ExptList';
eval(ds)
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\CrossOri\RandDirType2';
nexp = length(expt);
dirs = deg2rad(0:22.5:360-22.5);
angle = 120;
area = 'V1';

DSI_thresh = 0.3;

ind_ang  = find([expt.angle]==angle);
ind_area = find(strcmp([expt.img_loc],area));
ind = intersect(ind_ang,ind_area);

avg_resp_dir_all = [];
offset = 0;
resp_ind_all = cell(1,2);
u1_dir_all = [];
DSI_all = [];
TF_pref_all = [];
Zp_all = [];
Zc_all = [];

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
    resp_ind_all{1} = [resp_ind_all{1}; find(sum(h_resp(:,:,1),2)) + offset];
    resp_ind_all{2} = [resp_ind_all{2}; find(sum(h_resp(:,:,2),2)) + offset];
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']))
    u1_dir_all = cat(1,u1_dir_all, u1_dir);
    DSI_all = cat(1,DSI_all,DSI);
    TF_pref_all = [TF_pref_all TF_pref];
    Zp_all = [Zp_all; Zp];
    Zc_all = [Zc_all; Zc];
    offset = offset+size(avg_resp_dir,1);
end
save(fullfile(fnout,['crossOriRandDir_Type1Type2_' num2str(angle) 'deg_' area '_Summary.mat']), 'avg_resp_dir_all','resp_ind_all','DSI_all','u1_dir_all','TF_pref_all','Zp_all','Zc_all')
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))


figure(1)
for iTF = 1:nTF
    subplot(2,2,iTF)
    scatter(Zc_all(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),resp_ind_all{iTF}),iTF), Zp_all(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),resp_ind_all{iTF}),iTF))
    hold on
    scatter(Zc_all(intersect(find(TF_pref_all==iTF),intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),resp_ind_all{iTF})),iTF), Zp_all(intersect(find(TF_pref_all==iTF),intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),resp_ind_all{iTF})),iTF))
    title(['TF = ' num2str(TFs(iTF))])
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    plotZcZpBorders
end


prefDiff = cell(1,nTF);
ind_grat = cell(1,nTF);
ind_test = cell(1,nTF);
ind_av = cell(1,nTF);
ind_ioc  = cell(1,nTF);
ind_grat_av = cell(1,nTF);
ind_test_av = cell(1,nTF);
ind_rest = cell(1,nTF);
for iTF = 1:nTF
    if iTF == 1
        TF2 = 2;
    else
        TF2 = 1;
    end
    indUse = intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF});
    prefDiff{iTF} = u1_dir_all(:,iTF,1,1)-u1_dir_all(:,iTF,TF2,2);
    [ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),TFs(iTF),deg2rad(angle),TFs(TF2));
    figure(2) 
    subplot(3,2,iTF)
    scatter(rad2deg(u1_dir_all(indUse,iTF,1,1)),rad2deg(u1_dir_all(indUse,iTF,TF2,2)))
    refline(1)
    refline(1,-angle)
    refline(1,-rad2deg(ioc_ang))
    refline(1,-rad2deg(av_ang))
    xlabel('Pref dir- Grating')
    ylabel('Pref dir- Plaid')
    xlim([0 360])
    ylim([0 360])
    axis square
    title(['Grating- ' num2str(TFs(iTF)) ' Hz'])
    %legend('data',[num2str(TFs(iTF)) ' Hz comp'],[num2str(TFs(TF2)) ' Hz comp'],'IOC','VA','Location','northeastoutside')

    subplot(3,2,iTF+2)
    scale = 20;
    polarhistogram(prefDiff{iTF}(indUse),32);
    hold on
    polarplot([0;0]*pi/180,[0;1]*scale)
    polarplot([angle;angle]*pi/180,[0;1]*scale)
    polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
    polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)

    subplot(3,2,iTF+4)
    diff_thresh = 7.5;
    indUse = intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF});
    prefDiff_deg = rad2deg(prefDiff{iTF});
    ind_grat{iTF} = intersect(indUse, intersect(find(prefDiff_deg>-diff_thresh), find(prefDiff_deg<diff_thresh)));
    ind_test{iTF} = intersect(indUse,  intersect(find(prefDiff_deg>(angle-diff_thresh)), find(prefDiff_deg<(angle+diff_thresh))));
    ind_ioc{iTF} = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(ioc_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(ioc_ang)+diff_thresh))));
    ind_av{iTF} = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(av_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(av_ang)+diff_thresh))));
    if iTF == 2
        ind_grat_av{iTF} = intersect(ind_grat{iTF},ind_av{iTF});
        ind_grat{iTF} = setdiff(ind_grat{iTF}, ind_grat_av{iTF});
        ind_av{iTF} = setdiff(ind_av{iTF}, ind_grat_av{iTF});
        ind_rest{iTF} = setdiff(indUse,[ind_av{iTF}; ind_grat{iTF}; ind_grat_av{iTF}; ind_test{iTF}; ind_ioc{iTF}]);
        piechart([length(ind_grat{iTF}) length(ind_test{iTF}) length(ind_ioc{iTF}) length(ind_av{iTF}) length(ind_grat_av{iTF})],{[num2str(TFs(iTF)) ' Hz'],[num2str(TFs(TF2)) ' Hz'],'IOC','VA',[num2str(TFs(iTF)) ' Hz + AV']})
    else
        ind_test_av{iTF} = intersect(ind_test{iTF},ind_av{iTF});
        ind_test{iTF} = setdiff(ind_test{iTF}, ind_test_av{iTF});
        ind_av{iTF} = setdiff(ind_av{iTF}, ind_test_av{iTF});
        ind_rest{iTF} = setdiff(indUse,[ind_av{iTF}; ind_grat{iTF}; ind_test_av{iTF}; ind_test{iTF}; ind_ioc{iTF}]);
        piechart([length(ind_grat{iTF}) length(ind_test{iTF}) length(ind_ioc{iTF}) length(ind_av{iTF}) length(ind_test_av{iTF})],{[num2str(TFs(iTF)) ' Hz'],[num2str(TFs(TF2)) ' Hz'],'IOC','VA',[num2str(TFs(TF2)) ' Hz + AV']})
    end
end
print(fullfile(fnout,['crossOriRandDir_Type1Type2_' num2str(angle) 'deg_' area '_Summary.pdf']), '-bestfit','-dpdf')


figure; 
for iTF = 1:nTF
    subplot(2,1,iTF)
    scatter(Zc_all(ind_av{iTF},iTF), Zp_all(ind_av{iTF},iTF))
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    plotZcZpBorders
    title(['Grating- ' num2str(TFs(iTF)) ' Hz'])
end

scale = 15;
for iTF = 1:nTF
    if iTF == 1
        TF2 = 2;
    else
        TF2 = 1;
    end
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),TFs(iTF),deg2rad(angle),TFs(TF2));
figure
indUse = intersect(intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF}),find(Zp_all(:,iTF)-Zc_all(:,iTF)>-1.25));
subplot(2,2,1)
scatter(Zc_all(indUse,iTF), Zp_all(indUse,iTF))
title('Zp-Zc>0')
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
subplot(2,2,3)
p{iTF,1} = polarhistogram(prefDiff{iTF}(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*scale)
polarplot([angle;angle]*pi/180,[0;1]*scale)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)
indUse = intersect(intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF}),find(Zp_all(:,iTF)-Zc_all(:,iTF)<-1.25));
subplot(2,2,2)
scatter(Zc_all(indUse,iTF), Zp_all(indUse,iTF))
title('Zp-Zc<0')
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
hold on
plotZcZpBorders
subplot(2,2,4)
p{iTF,2} = polarhistogram(prefDiff{iTF}(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*scale)
polarplot([angle;angle]*pi/180,[0;1]*scale)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)
suptitle(['Grating TF ' num2str(TFs(iTF)) ' Hz'])
% subplot(3,1,3)
% plot(rad2deg(p{iTF,1}.BinEdges(1:32)),p{iTF,1}.Values)
% hold on
% plot(rad2deg(p{iTF,2}.BinEdges(1:32)),p{iTF,2}.Values)
% xlabel('Direction (deg)')
% legend({'Zp-Zc > -1.25','Zp-Zc < -1.25'})
print(fullfile(fnout,['crossOriRandDir_Type1Type2_' num2str(angle) 'deg_' area '_Summary_byZpZc_' num2str(TFs(iTF)) 'Hz.pdf']), '-bestfit','-dpdf')
end


figure;
scale = 15;
for iTF = 1:nTF
    if iTF == 1
        TF2 = 2;
    else
        TF2 = 1;
    end
[ioc_ang ioc_sp av_ang av_sp] = iocVectorCalc(deg2rad(0),TFs(iTF),deg2rad(angle),TFs(TF2));
indUse = intersect(intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF}),find(TF_pref_all==iTF));
subplot(2,2,1+((iTF-1)*2))
polarhistogram(prefDiff{iTF}(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*scale)
polarplot([angle;angle]*pi/180,[0;1]*scale)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)
title({['Pref TF- ' num2str(TFs(iTF)) ' Hz'], ['Grating TF- ' num2str(TFs(iTF)) ' Hz']})
indUse = intersect(intersect(intersect(find(DSI_all(:,iTF,1,1)>DSI_thresh),find(DSI_all(:,iTF,TF2,2)>DSI_thresh)),resp_ind_all{iTF}),find(TF_pref_all==TF2));
subplot(2,2,2+((iTF-1)*2))
polarhistogram(prefDiff{iTF}(indUse),32);
hold on
polarplot([0;0]*pi/180,[0;1]*scale)
polarplot([angle;angle]*pi/180,[0;1]*scale)
polarplot([rad2deg(ioc_ang);rad2deg(ioc_ang)]*pi/180,[0;1]*scale)
polarplot([rad2deg(av_ang);rad2deg(av_ang)]*pi/180,[0;1]*scale)
title({['Pref TF- ' num2str(TFs(TF2)) ' Hz'], ['Grating TF- ' num2str(TFs(iTF)) ' Hz']})
end
print(fullfile(fnout,['crossOriRandDir_Type1Type2_' num2str(angle) 'deg_' area '_Summary_prefTF.pdf']), '-bestfit','-dpdf')




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
ind_grat = intersect(indUse, intersect(find(prefDiff_deg>-diff_thresh), find(prefDiff_deg<diff_thresh)));
ind_test = intersect(indUse,  intersect(find(prefDiff_deg>(angle-diff_thresh)), find(prefDiff_deg<(angle+diff_thresh))));
ind_ioc = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(ioc_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(ioc_ang)+diff_thresh))));
ind_av = intersect(indUse,  intersect(find(prefDiff_deg>(rad2deg(av_ang)-diff_thresh)), find(prefDiff_deg<(rad2deg(av_ang)+diff_thresh))));
ind_grat_av = intersect(ind_grat,ind_av);
ind_grat = setdiff(ind_grat, ind_grat_av);
ind_av = setdiff(ind_av, ind_grat_av);

[prefResp prefDir] = max(avg_resp_dir_all(:,:,1,1),[],2);
nCells = size(avg_resp_dir_all,1);
avg_resp_dir_all_shift = zeros(size(avg_resp_dir_all));
for i = 1:nCells
    avg_resp_dir_all_shift(i,:,:,:) = circshift(avg_resp_dir_all(i,:,:,:),1-prefDir(i),2);
end
avg_resp_dir_all_shift_circ = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:));

figure(4);
subplot(2,2,1)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_grat,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_grat,:,2,1),1))
title(['Fast comp- n =' num2str(length(ind_grat))])
rlim([0 0.2])
subplot(2,2,2)
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_test,:,1,1),1))
hold on
polarplot([dirs dirs(1)], mean(avg_resp_dir_all_shift_circ(ind_test,:,2,1),1))
title(['Slow comp- n =' num2str(length(ind_test))])
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
piechart([length(ind_grat) length(ind_test) length(ind_ioc) length(ind_av) length(ind_grat_av)],{'Fast','Slow','IOC','VA','Fast+AV'})
title([num2str(angle) ' deg'])
print(fullfile(fnout,['crossOriRandDir_Type2_FractDir.pdf']), '-fillpage', '-dpdf')

