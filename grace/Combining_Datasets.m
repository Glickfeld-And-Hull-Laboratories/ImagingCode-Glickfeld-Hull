clear all;
clear global;
%% 
doRedChannel = 0;
% ds = 'CrossOriRandDir_ExptList_GL';
% eval(ds)
% rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
GL_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\grace';
summaryDir = fullfile(GL_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';

grna_list = strvcat('arc','lacz');
ngrna = size(grna_list,1);
% SF_list = strvcat('no','yes');
% nSF = size(SF_list,1);

% for igrna = 1:ngrna;
%     grna = grna_list(igrna,:);
%     fprintf([grna '\n'])
%     totCells = 0;
    avgRespD1all = [];
    avgRespD2all = [];
    avgRespD3all = [];
    maxRespD1all = [];
    maxRespD2all = [];
    maxRespD3all = [];
    signal_half1_all = [];
    signal_half2_all = [];
    signal_half3_all = [];
    signal_day_half2_all = [];
    signal_day_half3_all = [];
    prefOriD1all = [];
    prefOriD2all = [];
    prefOriD3all = [];
    k1_all = [];
    k2_all = [];
    k3_all = [];
    s1_all = [];
    s2_all = [];
    s3_all = [];
    h1_all = [];
    h2_all = [];
    h3_all = [];
    OSI1_all = [];
    OSI2_all = [];
    OSI3_all = [];
    DSI1_all = [];
    DSI2_all = [];
    DSI3_all = [];
    NR1_all = [];
    NR2_all = [];
    NR3_all = [];
    mouse_list = [];
    for iexp = 9:13
%         if sum(strcmp(expt(iexp).driver,grna))
            mouse = expt(iexp).mouse;
            mouse_list = strvcat(mouse_list, mouse);
            date = expt(iexp).date;
            ImgFolder = expt(iexp).coFolder;
%             time = expt(iexp).coTime;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);

            fprintf([mouse ' ' date '\n'])

            % load data
            corr = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_signalCorr.mat']));
            signal_half1 = corr.signal_corr_half1;
            signal_day_half2 = corr.signal_corr_day_half2;
            signal_day_half3 = corr.signal_corr_day_half3;
            OSI = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_OSI.mat']));
            osi1 = OSI.OSI1;
            osi2 = OSI.OSI2;
            osi3 = OSI.OSI3;
            DSI = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_DSI.mat']));
            dsi1 = DSI.DSI1;
            dsi2 = DSI.DSI2;
            dsi3 = DSI.DSI3;
            vonmises = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_vonmises.mat']));
            k1 = vonmises.k1_hat;
            k2 = vonmises.k1_hat2;
            k3 = vonmises.k1_hat3;
            thresh = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_threshold.mat']));
            s = thresh.S;
            s2 = thresh.S2;
            s3 = thresh.S3;
            prefOri = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_prefOri.mat']));
            prefOriD1 = prefOri.prefOri_D1;
            prefOriD2 = prefOri.prefOri_D2;
            prefOriD3 = prefOri.prefOri_D3;
            avg1 = prefOri.avgResp_D1;
            avgRespD1 = mean(avg1,2);
            maxRespD1 = max(avg1,[],2);
            avg2 = prefOri.avgResp_D2;
            avgRespD2 = mean(avg2,2);
            maxRespD2 = max(avg2,[],2);
            avg3 = prefOri.avgResp_D3;
            avgRespD3 = mean(avg3,2);
            maxRespD3 = max(avg3,[],2);
            ttest = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ttest.mat']));
            h1 = ttest.h1;
            h2 = ttest.h2;
            h3 = ttest.h3;
            newResp = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_newAvg.mat']));
            NR1 = newResp.newAvgD1;
            NR2 = newResp.newAvgD2;
            NR3 = newResp.newAvgD3;
            
            prefOriD1all = [prefOriD1all prefOriD1];
            prefOriD2all = [prefOriD2all prefOriD2];
            prefOriD3all = [prefOriD3all prefOriD3];
            avgRespD1all = [avgRespD1all avgRespD1'];
            avgRespD2all = [avgRespD2all avgRespD2'];
            avgRespD3all = [avgRespD3all avgRespD3'];
            maxRespD1all = [maxRespD1all maxRespD1'];
            maxRespD2all = [maxRespD2all maxRespD2'];
            maxRespD3all = [maxRespD3all maxRespD3'];
            NR1_all = vertcat(NR1_all,NR1);
            NR2_all = vertcat(NR2_all,NR2);
            NR3_all = vertcat(NR3_all,NR3);
            k1_all = [k1_all k1];
            k2_all = [k2_all k2];
            k3_all = [k3_all k3];
            s1_all = [s1_all s'];
            s2_all = [s2_all s2'];
            s3_all = [s3_all s3'];
            h1_all = [h1_all h1];
            h2_all = [h2_all h2];
            h3_all = [h3_all h3];
            signal_half1_all = [signal_half1_all signal_half1'];
            signal_day_half2_all = [signal_day_half2_all signal_day_half2'];
            signal_day_half3_all = [signal_day_half3_all signal_day_half3'];
            OSI1_all = [OSI1_all osi1'];
            OSI2_all = [OSI2_all osi2'];
            OSI3_all = [OSI3_all osi3'];
            DSI1_all = [DSI1_all dsi1'];
            DSI2_all = [DSI2_all dsi2'];
            DSI3_all = [DSI3_all dsi3'];
%             end
%         end
    end

%% 
deltaori = abs(prefOriD1all-prefOriD2all);
for ori = 1:length(deltaori)
if deltaori(ori)>90
    deltaori(ori) = 180-deltaori(ori);
else
    deltaori(ori) = deltaori(ori);
end
end

% deltaori2 = abs(prefOriD1all-prefOriD3all);
% for ori = 1:length(deltaori2)
% if deltaori2(ori)>90
%     deltaori2(ori) = 180-deltaori2(ori);
% else
%     deltaori2(ori) = deltaori2(ori);
% end
% end

deltamax = abs(maxRespD1all-maxRespD2all);
for ori = 1:length(deltamax)
if deltamax(ori)>90
    deltamax(ori) = 180-deltamax(ori);
else
    deltamax(ori) = deltamax(ori);
end
end

% deltamax2 = abs(maxRespD1all-maxRespD3all);
% for ori = 1:length(deltamax2)
% if deltamax2(ori)>90
%     deltamax2(ori) = 180-deltamax2(ori);
% else
%     deltamax2(ori) = deltamax2(ori);
% end
% end

%% converting to separate arc and lacz
k1_L = k1_all;
k2_L = k2_all;
h1_L = h1_all;
h2_L = h2_all;
s1_L = s1_all;
s2_L = s2_all;
signal_half1_L = signal_half1_all;
signal_day_half2_L = signal_day_half2_all;
signal_day_half3_L = signal_day_half3_all;
deltamaxL = deltamax;
deltaoriL = deltaori;
prefOriD1L = prefOriD1all;
prefOriD2L = prefOriD2all;
prefOriD3L = prefOriD3all;
OSI1L = OSI1_all;
OSI2L = OSI2_all;
OSI3L = OSI3_all;
DSI1L = DSI1_all;
DSI2L = DSI2_all;
DSI3L = DSI3_all;
NR1_L = NR1_all;
NR2_L = NR2_all;
NR3_L = NR3_all;
%% simple figures

% average tuning curve
figure;
subplot(1,2,1); 
fast_errbar(1:8,NR1_L,1,'color',[0.3 0 0.5]);hold on;fast_errbar(1:8,NR2_L,1,'color',[0.5430 0 0.5430]);hold on;fast_errbar(1:8,NR3_L,1,'color',[0.8633 0.6250 0.8633]);
fix_axes(gcf,16,'Orientation','dF/F'); axis square
subplot(1,2,2);
fast_errbar(1:8,NR1_all,1);hold on;fast_errbar(1:8,NR2_all,1,'color',[0.1328 0.5430 0.1328]);hold on;fast_errbar(1:8,NR3_all,1,'color',[0.5625 0.9297 0.5625]);
fix_axes(gcf,16,'Orientation','dF/F'); axis square

% DSI and OSI
sz = 100;
figure;scatter(OSI1L,OSI2L,sz,[0.5430 0 0.5430],'filled')
hold on;scatter(OSI1L,OSI3L,sz,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('OSI D1')
ylabel('OSI DX')
refline(1,0)
legend({'day 2' 'day 3'})
title('LacZ')

sz = 100;
figure;scatter(DSI1L,DSI2L,sz,[0.5430 0 0.5430],'filled')
hold on;scatter(DSI1L,DSI3L,sz,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('DSI D1')
ylabel('DSI DX')
refline(1,0)
legend({'day 2' 'day 3'})
title('LacZ')

sz = 100;
figure;scatter(OSI1_all(1:29),OSI2_all(1:29),sz,[0.1328 0.5430 0.1328],'filled')
hold on;scatter(OSI1_all(1:29),OSI3_all,sz,[0.5625 0.9297 0.5625],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('OSI D1')
ylabel('OSI DX')
refline(1,0)
legend({'day 2' 'day 3'})
title('Arc')

sz = 100;
figure;scatter(DSI1_all(1:29),DSI2_all(1:29),sz,[0.1328 0.5430 0.1328],'filled')
hold on;scatter(DSI1_all(1:29),DSI3_all,sz,[0.5625 0.9297 0.5625],'filled')
ax = gca;xlim([0,1]);ylim([0,1])
ax.FontSize = 20;
xlabel('DSI D1')
ylabel('DSI DX')
refline(1,0)
legend({'day 2' 'day 3'})
title('Arc')

sz = 100;
figure;scatter(signal_half1_all,signal_day_half2_all,sz,[0.5430 0 0.5430],'filled')
hold on;scatter(signal_half1_all,signal_day_half3_all,sz,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('signal correlation within D1')
ylabel('signal correlation D1 vs DX')
refline(1,0)
legend({'day 2' 'day 3'})

figure;scatter(maxRespD1all,maxRespD2all,sz,[0.5430 0 0.5430],'filled');hold on;scatter(maxRespD1all,maxRespD3all,sz,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('max resp D1')
ylabel('max resp DX')
refline(1,0)
legend({'day 2' 'day 3'})

figure;scatter(prefOriD1all,prefOriD2all,sz,[0.5430 0 0.5430],'filled');hold on;scatter(prefOriD1all,prefOriD3all,sz,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('pref ori D1')
ylabel('pref ori DX')
refline(1,0)

figure;
scatter(avgRespD1all,avgRespD2all,100,[0.2500 0.8750 0.8125],'filled');hold on;scatter(avgRespD1all,avgRespD3all,100,[0.9297 0.5078 0.9297],'filled')
xlabel('avg resp D1','FontSize',20);ylabel('avg resp DX','FontSize',20)
legend({'Day 2' 'Day 3'});refline(1,0)

% delta stuff
figure;
scatter(deltamax,deltaori,100,[0.5430 0 0.5430],'filled')
hold on
scatter(deltamax2,deltaori2,100,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('delta maxResp');ylabel('delta prefOri')
legend({'day 2' 'day 3'})

figure;
scatter(avgRespD1all,deltaori,100,'filled')
hold on
scatter(avgRespD1all,deltaori2,100,'filled')
xlabel('avgResp D1','FontSize',20);ylabel('delta prefOri','FontSize',20)
legend({'day 2' 'day 3'})

figure;
scatter(maxRespD1all,deltaori,100,'filled')
hold on
scatter(maxRespD1all,deltaori2,100,'filled')
xlabel('maxResp D1','FontSize',20);ylabel('delta prefOri','FontSize',20)
legend({'day 2' 'day 3'})

% ecdfs
figure;
cdfplot(k1_all);hold on;cdfplot(k2_all);hold on;cdfplot(k3_all)
ax = gca;
ax.FontSize = 20;
xlabel('k value');ylabel('fraction of cells')
legend({'day1','day2','day3'})

figure;
cdfplot(s_all);hold on;cdfplot(s2_all);hold on;cdfplot(s3_all)
ax = gca;
ax.FontSize = 20;
xlabel('selectivity');ylabel('fraction of cells')
legend({'day1','day2','day3'})

figure;
cdfplot(sum(h1_all,1));hold on;cdfplot(sum(h2_all,1));hold on;cdfplot(sum(h3_all,1))
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 20;
xlabel('ndirs')
ylabel('fraction of cells')
title('h value')

figure;
cdfplot(OSI1_all);hold on;cdfplot(OSI2_all);hold on;cdfplot(OSI3_all)
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 20;
xlabel('OSI')
ylabel('fraction of cells')
title('Arc')

figure;
cdfplot(OSI1L);hold on;cdfplot(OSI2L);hold on;cdfplot(OSI3L)
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 20;
xlabel('OSI')
ylabel('fraction of cells')
title('LacZ')   

figure;
cdfplot(DSI1_all);hold on;cdfplot(DSI2_all);hold on;cdfplot(DSI3_all)
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 20;
xlabel('DSI')
ylabel('fraction of cells')
title('Arc')

figure;
cdfplot(DSI1L);hold on;cdfplot(DSI2L);hold on;cdfplot(DSI3L)
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 20;
xlabel('DSI')
ylabel('fraction of cells')
title('LacZ')
%% arc vs lacz

figure;
scatter(OSI1_all,deltaori,100,[0.1328 0.5430 0.1328],'filled')
hold on
scatter(OSI1L,deltaoriL,100,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('delta maxResp');ylabel('delta prefOri')
legend({'arc' 'lacz'})

figure;
scatter(deltamax,deltaori,100,[0.1328 0.5430 0.1328],'filled')
hold on
scatter(deltamaxL,deltaoriL,100,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('delta maxResp');ylabel('delta prefOri')
legend({'arc' 'lacz'})

figure;
scatter(OSI1_all,OSI2_all,100,[0.1328 0.5430 0.1328],'filled')
hold on
scatter(OSI1L,OSI2L,100,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('OSI D1');ylabel('OSI D2')
legend({'arc' 'lacz'})

figure;
scatter(DSI1_all,DSI2_all,100,[0.1328 0.5430 0.1328],'filled')
hold on
scatter(DSI1L,DSI2L,100,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('DSI D1');ylabel('DSI D2')
legend({'arc' 'lacz'})

figure;
scatter(prefOriD1all,prefOriD2all,100,[0.1328 0.5430 0.1328],'filled')
hold on
scatter(prefOriD1L,prefOriD2L,100,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('pref ori D1');ylabel('pref ori D2')
legend({'arc' 'lacz'})

figure;
scatter(prefOriD1all,prefOriD3all,100,[0.5938 0.9833 0.5938],'filled')
hold on
scatter(prefOriD1L,prefOriD3L,100,[0.8633 0.6250 0.8633],'filled')
ax = gca;
ax.FontSize = 20;
refline(1,0)
xlabel('pref ori D1');ylabel('pref ori D3')
legend({'arc' 'lacz'})

sz = 60;
figure;scatter(signal_half1_all,signal_day_half2_all,sz,[0.1328 0.5430 0.1328],'LineWidth',1.3)
hold on;scatter(signal_half1_L,signal_day_half2_L,sz,[0.5430 0 0.5430],'LineWidth',1.3)
ax = gca;
ax.FontSize = 20;
ylim([-.5,1])
xlabel('signal correlation within D1')
ylabel('signal correlation D1 vs D2')
refline(1,0)
legend({'arc' 'lacz'})

sz = 60;
figure;scatter(signal_half1_all,signal_day_half3_all,sz,[0.5938 0.9833 0.5938],'LineWidth',1.3)
hold on;scatter(signal_half1_L,signal_day_half3_L,sz,[0.8633 0.6250 0.8633],'LineWidth',1.3)
ax = gca;
ax.FontSize = 20;
xlabel('signal correlation within D1')
ylabel('signal correlation D1 vs D3')
refline(1,0)
legend({'arc' 'lacz'})

sz = 100;
figure;scatter(k1_all,k2_all,sz,[0.1328 0.5430 0.1328],'filled');hold on;scatter(k1_L,k2_L,sz,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('k value D1')
ylabel('k value D2')
refline(1,0)
legend({'arc','lacz'})

sz = 100;
figure;scatter(sum(h1_all,1),sum(h2_all,1),sz,[0.1328 0.5430 0.1328],'filled');hold on;scatter(sum(h1_L,1),sum(h2_L,1),40,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('h value D1')
ylabel('h value D2')
refline(1,0)
legend({'arc','lacz'})

sz = 100;
figure;scatter(s1_all,s2_all,sz,[0.1328 0.5430 0.1328],'filled');hold on;scatter(s1_L,s2_L,sz,[0.5430 0 0.5430],'filled')
ax = gca;
ax.FontSize = 20;
xlabel('selectivity D1')
ylabel('selectivity D2')
refline(1,0)
legend({'arc','lacz'})

