close all;clear all;clc;

ds = 'CrossOriSingleStimAdapt_ExptList';
rc = behavConstsAV;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

eval(ds)
nexp = length(expt);

totCells = 0;

noadapt_resp_cell_all= cell(5,5);
h_noadapt_resp_all = [];
dir_resp_all = [];
dirresp_ind_all = [];
noadapt_resp_ind_all = cell(1,5);

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    fprintf([mouse ' ' date '\n']);
    run_str = ['runs-' cell2mat(expt(iexp).coFolder)];
    dir_run_str = ['runs-' cell2mat(expt(iexp).dirFolder)];
    
    fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
    dir_fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str]);
    load(fullfile(fn,[date '_' mouse '_' run_str '_allCellResp.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_dataStim.mat']));
    nMask = length(maskCons);
    nTest = length(testCons);

    for im = 1:nMask
        for it = 1:nTest
            noadapt_resp_cell_all{im,it} = [noadapt_resp_cell_all{im,it}; nanmean(noadapt_resp_cell{im,it},2)];
        end
    end
    nCells = size(noadapt_resp_cell{1,1},1);
    p = zeros(nCells,nTest,2);
    h_noadapt_resp = zeros(nCells,nTest,2);
    for i = 2:nTest
        [h_noadapt_resp(:,i,1) p(:,i,1)] = ttest2(noadapt_resp_cell{1,1}',noadapt_resp_cell{i,1}');
        [h_noadapt_resp(:,i,2) p(:,i,2)]= ttest2(noadapt_resp_cell{1,1}',noadapt_resp_cell{1,i}');
%         h_noadapt_resp(:,i,1) = p(:,i,1)<0.05./7;
%         h_noadapt_resp(:,i,2) = p(:,i,2)<0.05./7;
        noadapt_resp_ind_all{i} = [noadapt_resp_ind_all{i}; find(sum(h_noadapt_resp(:,i,:),3))+totCells];
    end
    h_noadapt_resp_all = cat(1,h_noadapt_resp_all,h_noadapt_resp);
    load(fullfile(dir_fn,[date '_' mouse '_' dir_run_str '_stimData.mat']));
    load(fullfile(dir_fn,[date '_' mouse '_' dir_run_str '_dfofData.mat']));
    dir_resp_all = cat(2, dir_resp_all, data_dfof_dir);
    dirresp_ind_all = [dirresp_ind_all  dirresp_ind+totCells];
    
    totCells = totCells+size(noadapt_resp_cell{1,1},1); 
end


%% cross ori measures
fnout = fullfile(LG_base, 'Analysis\2P\CrossOri\CrossOri_Figures');
if ~exist(fnout)
    mkdir(fnout);
end
noadapt_resp_cell_all_sub = cell(5,5);
for it = 1:nTest
    for im = 1:nMask
        noadapt_resp_cell_all_sub{im,it} = noadapt_resp_cell_all{im,it}-noadapt_resp_cell_all{1,1};
        noadapt_resp_cell_all_sub{im,it}(find(noadapt_resp_cell_all_sub{im,it}<0)) = 0;
    end
end

SI_all = zeros(totCells, nTest-1);
MI_all = zeros(totCells, nTest-1);
for it = 2:nTest
    SI_all(:,it) = abs((noadapt_resp_cell_all_sub{1,it}-noadapt_resp_cell_all_sub{it,1})./(noadapt_resp_cell_all_sub{1,it}+noadapt_resp_cell_all_sub{it,1}));
    MI_all(:,it) = (noadapt_resp_cell_all_sub{it,it}-(noadapt_resp_cell_all_sub{1,it}+noadapt_resp_cell_all_sub{it,1}))./...
        (noadapt_resp_cell_all_sub{it,it}+(noadapt_resp_cell_all_sub{1,it}+noadapt_resp_cell_all_sub{it,1}));
end

figure;
for it = 2:nTest
    subplot(2,nTest-1,it-1)
    histogram(SI_all(noadapt_resp_ind_all{it},it),[0:0.2:1])
    xlabel('Selectivity index')
    title(num2str(testCons(it)))
    subplot(2,nTest-1,it-1+nTest-1)
    histogram(MI_all(noadapt_resp_ind_all{it},it),[-1:0.2:1])
    xlabel('Masking index')
    title(num2str(chop(mean(MI_all(noadapt_resp_ind_all{it},it)),2)))
end
suptitle('All responsive cells')

[max_val max_dir] = max(squeeze(mean(dir_resp_all(61:90,:,:),1)),[],2);
ind0 = find(max_dir==1);
ind90 = find(max_dir==5);
ind_use = intersect(dirresp_ind_all,unique([ind90; ind0]));
figure;
for it = 2:nTest
    subplot(2,nTest-1,it-1)
    histogram(SI_all(intersect(ind_use,noadapt_resp_ind_all{it}),it),[0:0.2:1])
    xlabel('Selectivity index')
    title(num2str(testCons(it)))
    subplot(2,nTest-1,it-1+nTest-1)
    histogram(MI_all(intersect(ind_use,noadapt_resp_ind_all{it}),it),[-1:0.2:1])
    xlabel('Masking index')
    title(num2str(chop(mean(MI_all(intersect(ind_use,noadapt_resp_ind_all{it}),it)),2)))
end
suptitle('Test/Mask preferred direction')

figure;
for it = 2:nTest
    subplot(2,2,1)
    cdfplot(SI_all(noadapt_resp_ind_all{it},it))
    title('All responsive cells')
    xlabel('Selectivity index')
    hold on
    subplot(2,2,2)
    cdfplot(SI_all(intersect(ind_use,noadapt_resp_ind_all{it}),it))
    title('Test/Mask preferred direction')
    xlabel('Selectivity index')
    hold on
    subplot(2,2,3)
    cdfplot(MI_all(noadapt_resp_ind_all{it},it))
    title('All responsive cells')
    xlabel('Masking index')
    hold on
    subplot(2,2,4)
    cdfplot(MI_all(intersect(ind_use,noadapt_resp_ind_all{it}),it))
    title('Test/Mask preferred direction')
    xlabel('Masking index')
    hold on
end
subplot(2,2,1)
legend(num2str(testCons(2:end)'),'location','northwest')
print(fullfile(fnout,'CrossOri_ContrastCompare_SIMI.pdf'),'-dpdf')