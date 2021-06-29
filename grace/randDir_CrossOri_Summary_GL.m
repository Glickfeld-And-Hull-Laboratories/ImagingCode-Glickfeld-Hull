clear all;
clear global;
%% 
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList_GL';
eval(ds)
rc = behavConstsAV;
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
    totCells = 0;
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
    s_all = [];
    s2_all = [];
    s3_all = [];
    h_all = [];
    h2_all = [];
    h3_all = [];
    mouse_list = [];
    for iexp = 1:7
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
%             vonmises = load(fullfile(GL_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_vonmises.mat']));
%             k1 = vonmises.k1_hat;
%             k2 = vonmises.k1_hat2;
%             k3 = vonmises.k1_hat3;
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
            
            prefOriD1all = [prefOriD1all prefOriD1];
            prefOriD2all = [prefOriD2all prefOriD2];
            prefOriD3all = [prefOriD3all prefOriD3];
            avgRespD1all = [avgRespD1all avgRespD1'];
            avgRespD2all = [avgRespD2all avgRespD2'];
            avgRespD3all = [avgRespD3all avgRespD3'];
            maxRespD1all = [maxRespD1all maxRespD1'];
            maxRespD2all = [maxRespD2all maxRespD2'];
            maxRespD3all = [maxRespD3all maxRespD3'];
%             k1_all = [k1_all k1];
%             k2_all = [k2_all k2];
%             k3_all = [k3_all k3];
            s_all = [s_all s'];
            s2_all = [s2_all s2'];
            s3_all = [s3_all s3'];
            h1_all = [h_all h1'];
            h2_all = [h2_all h2];
            h3_all = [h3_all h3];
            signal_half1_all = [signal_half1_all signal_half1'];
            signal_day_half2_all = [signal_day_half2_all signal_day_half2'];
            signal_day_half3_all = [signal_day_half3_all signal_day_half3'];
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

deltaori2 = abs(prefOriD1all-prefOriD3all);
for ori = 1:length(deltaori2)
if deltaori2(ori)>90
    deltaori2(ori) = 180-deltaori2(ori);
else
    deltaori2(ori) = deltaori2(ori);
end
end

deltamax = abs(maxRespD1all-maxRespD2all);
for ori = 1:length(deltamax)
if deltamax(ori)>90
    deltamax(ori) = 180-deltamax(ori);
else
    deltamax(ori) = deltamax(ori);
end
end

deltamax2 = abs(maxRespD1all-maxRespD3all);
for ori = 1:length(deltamax2)
if deltamax2(ori)>90
    deltamax2(ori) = 180-deltamax2(ori);
else
    deltamax2(ori) = deltamax2(ori);
end
end

%% simple stuff
sz = 100;
figure;scatter(signal_half1_all,signal_day_half2_all,sz,'filled')
hold on;scatter(signal_half1_all,signal_day_half3_all,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('signal correlation within D1')
ylabel('signal correlation D1 vs DX')
refline(1,0)
legend({'within D1 session' 'between sessions'})

figure;scatter(maxRespD1all,maxRespD2all,sz,'filled');hold on;scatter(maxRespD1all,maxRespD3all,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('max resp D1')
ylabel('max resp DX')
refline(1,0)
legend({'day 2' 'day 3'})

figure;scatter(prefOriD1all,prefOriD2all,sz,'filled');hold on;scatter(prefOriD1all,prefOriD3all,sz,'filled')
ax = gca;
ax.FontSize = 20;
xlabel('pref ori D1')
ylabel('pref ori DX')
legend({'Day 2' 'Day 3'})
refline(1,0)

figure;
scatter(avgRespD1all,avgRespD2all,60,'filled');hold on;scatter(avgRespD1all,avgRespD3all,60,'filled')
xlabel('avg resp D1','FontSize',20);ylabel('avg resp DX','FontSize',20)
legend({'Day 2' 'Day 3'});refline(1,0)

% delta stuff
figure;
scatter(deltamax,deltaori,60,'filled')
hold on
scatter(deltamax2,deltaori2,60,'filled')
xlabel('delta maxResp','FontSize',20);ylabel('delta prefOri','FontSize',20)
legend({'Day 2' 'Day 3'})

figure;
scatter(avgRespD1all,deltaori,60,'filled')
hold on
scatter(avgRespD1all,deltaori2,60,'filled')
xlabel('avgResp D1','FontSize',20);ylabel('delta prefOri','FontSize',20)
legend({'Day 2' 'Day 3'})

figure;
scatter(maxRespD1all,deltaori,60,'filled')
hold on
scatter(maxRespD1all,deltaori2,60,'filled')
xlabel('maxResp D1','FontSize',20);ylabel('delta prefOri','FontSize',20)
legend({'Day 2' 'Day 3'})

% ecdfs
figure;
ecdf(k_all,'Bounds','on');hold on;ecdf(k2_all,'Bounds','on');hold on;ecdf(k3_all,'Bounds','on');
xlabel('k value');ylabel('fraction of cells')
legend({'goodCells'})

figure;
cdfplot(s_all);hold on;cdfplot(s2_all);hold on;cdfplot(s3_all)
ax = gca;
ax.FontSize = 16;
xlabel('selectivity');ylabel('fraction of cells')
legend({'day1','day2','day3'})

figure;
cdfplot(sum(h1_all,1));hold on;cdfplot(sum(h2_all,1));hold on;cdfplot(sum(h3_all,1))
legend({'day1','day2','day3'})
ax = gca;
ax.FontSize = 16;
xlabel('ndirs')
ylabel('fraction of cells')
title('h value')
        
sz = 100;
figure;scatter(k1_A(cell_list),k2_A,sz,'filled');hold on;scatter(k1_L(cell_list),k2_L,sz,'filled')
ax = gca;
ax.FontSize = 16;
xlabel('k value D1','FontSize',20)
ylabel('k value D2','FontSize',20)
refline(1,0)
legend({'arc','lacz'})

sz = 100;
figure;scatter(s1_A,s2_A,sz,'filled');hold on;scatter(s1_L,s2_L,sz,'filled')
ax = gca;
ax.FontSize = 16;
xlabel('selectivity D1','FontSize',20)
ylabel('selectivity D2','FontSize',20)
refline(1,0)
legend({'arc','lacz'})

