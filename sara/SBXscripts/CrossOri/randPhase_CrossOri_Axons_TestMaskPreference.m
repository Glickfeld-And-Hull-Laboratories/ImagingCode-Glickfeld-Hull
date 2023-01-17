%Validating V1 boutons compared to V1 somas

clear all; close all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'axons');
svName = 'randPhase';
dateOfAnalysis = '220614';
driver_list = strvcat('SLC', 'SLC', 'syn'); 
area = 'V1';
inj_list = strvcat('V1', 'LM', 'V1');
narea = length(inj_list);

figure;
    for iA = 1:narea
        
        fprintf([inj_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_' inj_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[inj_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(testPI_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('Test/Mask Preference Index')     
    end
 movegui('center')


%obtain lists of responsive somas and boutons in V1
for iA = 1:narea
    load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_' inj_list(iA,:) '_' driver_list(iA,:) '.mat'])))
    
    ind = resp_ind_all;

    if iA == 1
        somaPI = testPI_all(ind,:); %PIs of responsive somas
    else 
        boutonPI = testPI_all(ind,:); %PIs of responsive boutons
    end
end

%find the indices of the closest bouton PIs to every soma PI
index = zeros(280,1);
sampInd = [];
for i = 1:length(somaPI)
    [c, index] = min(abs(boutonPI-somaPI(i)));
    boutonPI(index) = NaN;
    sampInd = [sampInd; index];
end 

%Plot PI cdf and phase modulation amplitude cdf of the subsampled group compared to V1 and LM somas
for iA = 1:narea
    load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_' inj_list(iA,:) '_' driver_list(iA,:) '.mat'])))
    ind = resp_ind_all;
    nCells = length(ind);
    leg_str{iA}=[inj_list(iA,:) ' n=' num2str(length(ind))];
        
    subplot(2,2,2)
        if iA == 3
            leg_str{iA}=[inj_list(iA,:) 'boutons n=' num2str(length(sampInd))];
            amp_all_all = amp_all_all(ind,:);
            cdfplot(amp_all_all(sampInd,:));
        else
            cdfplot(amp_all_all(ind,:));            
        end
            hold on
            legend(leg_str, 'location', 'southeast')
            ylabel('Fraction of cells')
            xlabel('Phase modulation amplitude')
            title('V1 projections subsampled to match V1 soma PIs')
            
    subplot(2,2,3)
            if iA == 3
                boutonPI = testPI_all(ind,:);
                leg_str{iA}=[inj_list(iA,:) 'boutons n=' num2str(length(sampInd))];
                cdfplot(boutonPI(sampInd,:));
            else    
                cdfplot(testPI_all(ind,:));
            end
        hold on
        legend(leg_str, 'location', 'southeast')
        ylabel('Fraction of cells')
        xlabel('Test/Mask Preference Index')
        title('V1 projections subsampled to match V1 soma PIs')            
end
 movegui('center')

    print(fullfile(summaryDir, [svName '_' area '_CompareTestMaskPI.pdf']),'-dpdf', '-fillpage') 

