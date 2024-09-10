close all; clear all; clc;
doRedChannel = 1;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'PopulationDecoding');

doPlot = 1;
ds = ['CrossOriRandDirFourPhase_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SLC';
img_area = {'V1';'L2/3'}; %LM
inj_area = 'V1';
img_layer = 'L2/3';

max_dist = 5;

rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

mouse_list = [];
totCells = zeros(nexp,1);
ind_all=[];
max_avg_grat_all=[];
max_avg_plaid_all=[];
pref_dirs_all=[];

pattind_all = [];

%V1 L2/3 - 5 24 46 47 79
%V1 L4 - 62 66
%LM L2/3 - 38 49 50 77
% AL L2/3 - 76 82 83 

start=1;
for iexp = [5 24 46 47 79]
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    if isfield(expt,'copFolder') 
        ImgFolder = expt(iexp).copFolder;
    else
        ImgFolder = expt(iexp).coFolder;
    end
    nrun = length(ImgFolder);
    run_str = 'runs-002';
        
    fprintf([mouse ' ' date '\n'])

    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))    
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_bootstrapFits.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']))
    

    nCells = size(avg_resp_dir,1);
    totCells(iexp,:) = nCells;
    ind_all = [ind_all; ind+sum(totCells(1:iexp-1,:),1)];

    max_avg_grat = [];
    max_avg_plaid = [];
    pref_dirs=[];
    for ic = 1:nCells
        [max_avg_grat(ic), pref_dirs(ic)] = max(avg_resp_dir(ic,:,1,1,1));
        [max_avg_plaid(ic), pref_dirs_p(ic)] = max(avg_resp_dir(ic,:,1,2,1));
    end
    
    max_avg_grat_all = [max_avg_grat_all, max_avg_grat];
    max_avg_plaid_all = [max_avg_plaid_all, max_avg_plaid];
    pref_dirs_all = [pref_dirs_all, pref_dirs];

    for ic = 1:nCells
        avg_resp_dir_all((ic+sum(totCells(1:iexp-1,:),1)),:,:) = avg_resp_dir(ic,:,1,:,1);
    end

    stimDirs=[0    30    60    90   120   150   180   210   240   270   300   330];
    nStimDir=length(stimDirs);
    nMaskPhas=1;
    int = unique(diff(stimDirs));
    component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
    pattern = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);
    
    %Calculate pattern and component prediction
    comp_corr = zeros(nMaskPhas,nCells);
    patt_corr = zeros(nMaskPhas,nCells);
    comp_patt_corr = zeros(nMaskPhas,nCells);
    plaid_corr = zeros(1,nCells);
    plaid_corr_rand = zeros(1,nCells);
    
    for iCell = 1:nCells
        for ip = 1:nMaskPhas
            comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
            patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
            comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
        end
    end
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));
    
    ZcZp_diff = Zc-Zp;
    pattind = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));

    pattind_all = [pattind_all, pattind+sum(totCells(1:iexp-1,:),1)];

end

if doPlot == 0
    stop 
end


%% Plot pref directions

    figure;
    subplot(4,2,1)
        histogram(pref_dirs_all)
        hold on
        ylabel('# of cells')
        xlabel('direction')
        subtitle(['All neurons n=' num2str(length(pref_dirs_all))])
    subplot(4,2,2)
        histogram(pref_dirs_all(ind_all))
        hold on
        ylabel('# of cells')
        xlabel('direction')
        subtitle(['vis resp/DS neurons only n=' num2str(length(ind_all))])    
   
    movegui('center')
    sgtitle('Preferred directions')
    print(fullfile(outDir, ['OnePhaseData_' inj_area '_' driver '_PrefDirHistograms.pdf']),'-dpdf', '-fillpage')



%% Plot population tuning for grating direction 

 %All neurons
    for id=1:12
        for ipd = 1:12
            if 
            prefneurons = find(pref_dirs_all==ipd);
            vector(id,ipd) = sum(avg_resp_dir_all(prefneurons,id,1)'./max_avg_grat_all(prefneurons))/length(prefneurons);
        end
    end

    figure;
        x=[-150:30:180];
        x_rad = deg2rad(x);
        for id = 1:12
           subplot(3,4,id)
                polarplot([x_rad x_rad(1)], [vector(id,:) vector(id,1)])
                hold on
        end

    movegui('center')
    sgtitle(['Grating direction decoding n=' num2str(length(pref_dirs_all))])
    print(fullfile(outDir, ['OnePhaseData_' inj_area '_' driver '_DecodeGratingDirection_RateCode_AllNeurons.pdf']),'-dpdf', '-fillpage')


 %Inclusion criteria
    for id=1:12
        for ipd = 1:12
            prefneurons = intersect(find(pref_dirs_all==ipd), ind_all);
            vector(id,ipd) = sum(avg_resp_dir_all(prefneurons,id,1)'./max_avg_grat_all(prefneurons))/length(prefneurons);
        end
    end

    figure;
        x=[-150:30:180];
        x_rad = deg2rad(x);
        for id = 1:12
           subplot(3,4,id)
                polarplot([x_rad x_rad(1)], [vector(id,:) vector(id,1)])
                hold on
        end

    movegui('center')
    sgtitle(['Grating direction decoding n=' num2str(length(ind_all))])
    print(fullfile(outDir, ['OnePhaseData_' inj_area '_' driver '_DecodeGratingDirection_RateCode_DSonly.pdf']),'-dpdf', '-fillpage')



%% Plot population tuning for plaid direction 


 %All neurons
    for id=1:12
        for ipd = 1:12
            if ipd<11 
                prefneurons = find(pref_dirs_all==(ipd+2));
            else
                prefneurons = find(pref_dirs_all==(ipd-10));
            end
            vector(id,ipd) = sum(avg_resp_dir_all(prefneurons,id,2)'./max_avg_plaid_all(prefneurons))/length(prefneurons);
        end
    end

    figure;
        x=[-150:30:180];
        x_rad = deg2rad(x);
        for id = 1:12
           subplot(3,4,id)
                polarplot([x_rad x_rad(1)], [vector(id,:) vector(id,1)])
                hold on
        end

    movegui('center')
    sgtitle(['Plaid direction decoding n=' num2str(length(pref_dirs_all))])
    % print(fullfile(outDir, ['OnePhaseData_' inj_area '_' driver '_DecodeGratingDirection_RateCode_AllNeurons.pdf']),'-dpdf', '-fillpage')


 %Inclusion criteria
    for id=1:12
        for ipd = 1:12
            if ipd<11 
                prefneurons = intersect(find(pref_dirs_all==(ipd+2)), ind_all);
            else
                prefneurons = intersect(find(pref_dirs_all==(ipd-10)), ind_all);
            end
            vector(id,ipd) = sum(avg_resp_dir_all(prefneurons,id,2)'./max_avg_plaid_all(prefneurons))/length(prefneurons);
        end
    end

    figure;
        x=[-150:30:180];
        x_rad = deg2rad(x);
        for id = 1:12
           subplot(3,4,id)
                polarplot([x_rad x_rad(1)], [vector(id,:) vector(id,1)])
                hold on
        end

    movegui('center')
    sgtitle(['Plaid direction decoding n=' num2str(length(ind_all))])
    %print(fullfile(outDir, ['OnePhaseData_' inj_area '_' driver '_DecodeGratingDirection_RateCode_DSonly.pdf']),'-dpdf', '-fillpage')
 
    
 %% inclusion criteria related to being PDS/CDS/unclassified

pattind_all = pattind_all';

    for id=1:12
        for ipd = 1:12
            if ipd<11 
                prefneurons = intersect(find(pref_dirs_all==(ipd+2)), pattind_all);
            else
                prefneurons = intersect(find(pref_dirs_all==(ipd-10)), pattind_all);
            end
            vector(id,ipd) = sum(avg_resp_dir_all(prefneurons,id,2)'./max_avg_plaid_all(prefneurons))/length(prefneurons);
        end
    end

    figure;
        x=[-150:30:180];
        x_rad = deg2rad(x);
        for id = 1:12
           subplot(3,4,id)
                polarplot([x_rad x_rad(1)], [vector(id,:) vector(id,1)])
                hold on
        end

    movegui('center')
    sgtitle(['Plaid direction decoding n=' num2str(length(pattind_all))])


