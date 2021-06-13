close all;clear all;clc;

%ds = 'CrossOriSingleStimAdapt_ExptList';
ds = 'CrossOriRandPhase_lowSF_ExptList';

rc = behavConstsAV;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

eval(ds)
nexp = length(expt);

totCells = 0;

resp_cell_all= cell(3,3);
h_resp_all = [];
resp_ind_all = cell(1,5);

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    fprintf([mouse ' ' date '\n']);
    run_str = ['runs-' cell2mat(expt(iexp).coFolder)];
    
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    
    nMask = length(maskCons);
    nTest = length(stimCons);

    for im = 1:nMask
        for it = 1:nTest
            if size(resp_cell{im,it,1},1)>1
                resp_cell_all{im,it} = [resp_cell_all{im,it}; nanmean(resp_cell{im,it,1},2)];
            else
                resp_cell_all{im,it} = [resp_cell_all{im,it}; resp_cell{im,it,1}'];
            end
        end
    end
    nCells = size(resp_cell{nTest,nTest,1},1);
    p = zeros(nCells,nTest,2);
    h_resp = zeros(nCells,nTest,2);
    for i = 2:nTest
        if size(resp_cell{1,1,1},1)>1
            [h_resp(:,i,1) p(:,i,1)] = ttest2(resp_cell{1,1,1}',resp_cell{i,1,1}');
            [h_resp(:,i,2) p(:,i,2)]= ttest2(resp_cell{1,1,1}',resp_cell{1,i,1}');
        else
            [h_resp(:,i,1) p(:,i,1)] = ttest(resp_cell{i,1,1}');
            [h_resp(:,i,2) p(:,i,2)]= ttest(resp_cell{1,i,1}');
        end
        resp_ind_all{i} = [resp_ind_all{i}; find(sum(h_resp(:,i,:),3))+totCells];
    end
    h_resp_all = cat(1,h_resp_all,h_resp);
    
    totCells = totCells+size(resp_cell{nTest,nTest,1},1); 
end


%% cross ori measures
fnout = fullfile(LG_base, 'Analysis\2P\CrossOri\CrossOri_Figures');
if ~exist(fnout)
    mkdir(fnout);
end
resp_cell_all_sub = cell(5,5);
for it = 1:nTest
    for im = 1:nMask
        resp_cell_all_sub{im,it} = resp_cell_all{im,it}-resp_cell_all{1,1};
        resp_cell_all_sub{im,it}(find(resp_cell_all_sub{im,it}<0)) = 0;
    end
end

SI_all = zeros(totCells, nTest-1);
MI_all = zeros(totCells, nTest-1);
for it = 2:nTest
    SI_all(:,it) = abs((resp_cell_all_sub{1,it}-resp_cell_all_sub{it,1})./(resp_cell_all_sub{1,it}+resp_cell_all_sub{it,1}));
    MI_all(:,it) = (resp_cell_all_sub{it,it}-(resp_cell_all_sub{1,it}+resp_cell_all_sub{it,1}))./...
        (resp_cell_all_sub{it,it}+(resp_cell_all_sub{1,it}+resp_cell_all_sub{it,1}));
end

figure;
for it = 2:nTest
    subplot(2,nTest-1,it-1)
    histogram(SI_all(resp_ind_all{it},it),[0:0.2:1])
    xlabel('Selectivity index')
    title(num2str(stimCons(it)))
    subplot(2,nTest-1,it-1+nTest-1)
    histogram(MI_all(resp_ind_all{it},it),[-1:0.2:1])
    xlabel('Masking index')
    title(num2str(chop(nanmean(MI_all(resp_ind_all{it},it)),2)))
end
suptitle('All responsive cells')

figure;
resp_all = zeros(totCells,nTest-1);
for it = 2:nTest
    subplot(2,2,1)
    cdfplot(SI_all(resp_ind_all{it},it))
    title('All responsive cells')
    xlabel('Selectivity index')
    hold on
    subplot(2,2,2)
    cdfplot(MI_all(resp_ind_all{it},it))
    title('All responsive cells')
    xlabel('Masking index')
    hold on
    subplot(2,2,3)
    resp_all(:,it) = max([resp_cell_all_sub{1,it} resp_cell_all_sub{it,1}],[],2);
    cdfplot(resp_all(:,it))
    xlabel('dF/F')
    title('Response amplitude- All cells')
    hold on
    subplot(2,2,4)
    cdfplot(resp_all(resp_ind_all{it},it))
    xlabel('dF/F')
    title('Response amplitude- Resp cells')
    hold on
end
subplot(2,2,1)
legend(num2str(stimCons(2:end)'),'location','northwest')

[h p_MI] = ttest2(MI_all(resp_ind_all{2},2),MI_all(resp_ind_all{3},3));
[h p_resp] = ttest2(resp_all(resp_ind_all{2},2),resp_all(resp_ind_all{3},3));

print(fullfile(fnout,'CrossOriRandPhase_ContrastCompare_SIMI.pdf'),'-dpdf')