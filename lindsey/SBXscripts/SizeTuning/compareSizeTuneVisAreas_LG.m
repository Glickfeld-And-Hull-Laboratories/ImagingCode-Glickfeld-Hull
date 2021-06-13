%% compare size-tuning across visual areas
% script that loads given experiments and sorts by visual area to compare
% size tuning parameters across V1 + HVAs

% given csv .txt file for experiment list
% includes: date, mouse, run_str(size-tuning exp), HVA, indicator

% data to load:
% ?ret results - lbub_fits (w/ RF location), goodfit cells (maybe?)
% ?size tuning run input - stim location (maybe?)
% _sizeTuneData.mat - sizeTune, sizeMean, sizeSEM, cellDists
% _sizeFitResults_SP.mat - sizeFits struct (all cons)
% _Fit_struct.mat - Fit_struct (highest con, with shuffles)
% _lbub_fits - lbub_fits, goodfit_ind_size

%% load csv to define exp

clear all; clc;
close all;
ds = 'retAndSzTuning_ExptList';
cellType = 'SLC';
rc = behavConstsAV;
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
fndata = fullfile(lg_fn, 'Analysis\2P');
fnout = fullfile(fndata,'SizeTuning',['SzTuning_' cellType]);
if ~exist(fnout)
    mkdir(fnout)
end
eval(ds)
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);
nExp = size(expt,2);
fprintf(['Size-tuning visual-area comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading sizeTuneData (raw size-tuning data)\n')
sizeTune_all = cell(0);
sizeMean_all = [];
sizeSEM_all = [];
cellDists_all = [];
nCellsExp = zeros(1,nExp);
sizeFits_all = []; 
fitout_all = []; 
fitC1_all = []; 
fitC2_all = []; 
data_all = [];
lbub_fits_all = [];
goodfit_ind_size_all = [];
szs0 = [];

t = 0;
for i=1:nExp
    if strcmp(cellType,expt(i).driver{1})
        fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
        date = expt(i).date;
        mouse = expt(i).mouse;
        run_str = cell2mat(expt(i).szFolder);
        filename = fullfile(fndata, [date '_' mouse], [date '_' mouse '_runs-' run_str], [date '_' mouse '_runs-' run_str '_Tuning.mat']);
        if ~exist(filename, 'file')
            fprintf([[date '_' mouse '_' run_str '_Tuning.mat'] ' not found! Please remove from list\n'])
        end
        load(filename, 'sizeTune', 'tuning_mat', 'cellDists')
        filename = fullfile(fndata, [date '_' mouse], [date '_' mouse '_runs-' run_str], [date '_' mouse '_runs-' run_str '_input.mat']);
        load(filename);
        filename = fullfile(fndata, [date '_' mouse], [date '_' mouse '_runs-' run_str], [date '_' mouse '_runs-' run_str '_Fit_struct.mat']);
        if ~exist(filename, 'file')
            fprintf([[date '_' mouse '_' run_str '_Fit_struct.mat'] ' not found! Please remove from list\n'])
        end
        load(filename, 'Fit_struct')
        filename = fullfile(fndata, [date '_' mouse], [date '_' mouse '_runs-' run_str], [date '_' mouse '_runs-' run_str '_lbub_fits.mat']);
        if ~exist(filename, 'file')
            fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
        end
        load(filename, 'lbub_fits', 'goodfit_ind_size')


        sizeTune_all = cat(3,sizeTune_all,sizeTune);
        sizeMean = squeeze(tuning_mat(:,:,1,:));
        sizeSEM = squeeze(tuning_mat(:,:,2,:));
        sizeMean_all = cat(3,sizeMean_all,sizeMean);
        sizeSEM_all = cat(3,sizeSEM_all,sizeSEM);
        cellDists_all = [cellDists_all;cellDists];
        nCellsExp(i) = length(cellDists);

        cons = unique(celleqel2mat_padded(input.tGratingContrast));
        nCon = length(cons);
        fit_true_vec = NaN(nCellsExp(i),11,nCon);
        fitout = NaN(nCellsExp(i),nCon,2,100);
        fitC1 = NaN(nCellsExp(i),nCon,3);
        fitC2 = NaN(nCellsExp(i),nCon,6);
        data = cell(nCellsExp(i),nCon,2);
        for iCon = 1:nCon
            for iCell = 1:nCellsExp(i)
                if ~isempty(Fit_struct(iCell,iCon).True)
                    eval('tmp = Fit_struct(iCell,iCon).True.s_.prefSize;');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.prefSize1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.prefSize2];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.suppInd2];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.Fscore];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.Ftest];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.maxResp1];');
                    eval('tmp = [tmp Fit_struct(iCell,iCon).True.s_.maxResp2];');
                    eval('tmp = [tmp exist(Fit_struct(iCell,iCon).True.s_.Fstr1)];');
                    eval('tmp2 = Fit_struct(iCell,iCon).True.s_.fitout1;');
                    eval('tmp2 = [tmp2; Fit_struct(iCell,iCon).True.s_.fitout2];');
        % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest maxR1 maxR2 Fstr
                    fit_true_vec(iCell,:,iCon) = tmp;
                    fitout(iCell,iCon,:,:) = tmp2;
                    fitC1(iCell,iCon,:) = Fit_struct(iCell,iCon).True.s_.fit1.c1;
                    fitC2(iCell,iCon,:) = Fit_struct(iCell,iCon).True.s_.fit2.c2;
                    data{iCell,iCon,1} = Fit_struct(iCell,iCon).True.s_.data;
                    data{iCell,iCon,2} = Fit_struct(iCell,iCon).True.s_.szs0;
                end
            end
        end
        sizeFits_all = cat(1,sizeFits_all,fit_true_vec);
        fitout_all = cat(1,fitout_all,fitout);
        fitC1_all = cat(1,fitC1_all,fitC1);
        fitC2_all = cat(1,fitC2_all,fitC2);
        data_all = cat(1,data_all,data);
        lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
        tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind_size; % offset by # cells in previous exps
        goodfit_ind_size_all = [goodfit_ind_size_all tempinds];

        fprintf('done\n')
        
    end
end

m1 = Fit_struct(1,1).True.s_.m1;
m2 = Fit_struct(1,1).True.s_.m2;
nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit_size)\n',nCellsTot,length(goodfit_ind_size_all))
expInd = [];
expArea = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
    expArea = [expArea repmat(expt(i).img_loc, 1,nCellsExp(i))];
end

szs = unique(celleqel2mat_padded(input.tGratingDiameterDeg)); 
nSize = length(szs);
szRng = linspace(0,max(szs));

filename = fullfile(fnout, ['sizeDataAll'  cellType '.mat']);
save(filename,'sizeFits_all','fitout_all', 'fitC1_all', 'fitC2_all', 'data_all', 'lbub_fits_all', 'goodfit_ind_size_all', 'nCellsTot', 'm1', 'm2', 'expInd', 'expArea', 'szs', 'nSize', 'szRng');
%% con resp for 4 con data
fprintf('Extracting contrast response of each cell at prefSize\n')

s_all = zeros(1,nCon);
s = zeros(1);
conStruct_all = struct('resp',s_all,'fit',s_all,'C50r',s,'Rsq',s,'x0',s_all);
conStruct_all(nCellsTot) = conStruct_all;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
usePrefSize = 1;
useRFSize = 0;

for iCell=1:nCellsTot
    if ~sum(iCell==goodfit_ind_size_all)
        if sum(iCell==[1 nCellsTot]) % check first and last to reset zeros to blank
            conStruct_all(iCell).resp = [];conStruct_all(iCell).fit = [];conStruct_all(iCell).C50r = [];conStruct_all(iCell).Rsq = [];conStruct_all(iCell).x0 = [];
        end
        continue % do not fit unless goodfit_size
    end
    
    if usePrefSize
        pS = sizeFits_all(iCell,1,nCon);
        pSind = find(szRng==pS);
    elseif useRFSize
        switch cell2mat(expArea(iCell))
            case 'V1'
                RFsize = 20; 
            case {'LM','AL'}
                RFsize = 20;
            case 'PM'
                RFsize = 20;
        end
        pS = RFsize;
        [val pSind] = min(abs(szRng-pS));
    end
    
    for iCon = 1:nCon
        if sizeFits_all(iCell,8,iCon)
            conStruct_all(iCell).resp(iCon) = fitout_all(iCell,iCon,2,pSind);
        else
            conStruct_all(iCell).resp(iCon) = fitout_all(iCell,iCon,1,pSind);
        end
    end
    
    cRi = conStruct_all(iCell).resp;
    lb = [0 0 0.1 1];
    ub = [Inf Inf 0.8 Inf];
    SStot = sum((cRi-mean(cRi)).^2);
    R2best = -Inf;
    for i=1%1:4
        x0 = [cRi(1) max(cRi) 0.1+0.1*i 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,cons,cRi,lb,ub,opts);
        R2 = 1-res/SStot;
        if R2>R2best
            R2best = R2;
            cFbest = cF;
            x0best = x0;
        end
    end
    cF = cFbest;
    R2 = R2best;
    
    conStruct_all(iCell).fit = cF;
    conStruct_all(iCell).Rsq = R2;
    conStruct_all(iCell).x0 = x0best;
    
    fitoutz = conModelH(cF,conRng);
    R50 = fitoutz(1)+(fitoutz(end)-fitoutz(1))/2;
    fitout50rect = abs(fitoutz - R50);
    i50 = find(fitout50rect == min(fitout50rect),1);
    C50 = conRng(i50);
    conStruct_all(iCell).C50r = C50;
    
    fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCellsTot,cF(1),cF(2),cF(3),cF(4),R2,C50)
   
end
fprintf('Done, saving...\n')
filename = fullfile(fnout, ['conStruct_'  cellType '.mat']);
save(filename,['conStruct_all']);

%% load contrast response instead of compute
filename = fullfile(fnout, ['conStruct_' cellType '.mat']);
load(filename);
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));

filename = fullfile(fnout, ['sizeDataAll'  cellType '.mat']);
load(filename);
%% present each exp to examine cells+fits and choose example cells
% show each cell rawdata w/ model1+2 overlaid

close all
for iExp = unique(expInd)
    % select only cells in exp i with goodfits and RF<10
    
    switch expt(iExp).img_loc{1}
        case 'V1'
            RFcutoff = find(cellDists_all<15);
        case {'LM', 'AL'}
            RFcutoff = find(cellDists_all<15);
        case 'PM'
            RFcutoff = find(cellDists_all<15);
    end
    
    inds = intersect(intersect(find(ismember(expInd,iExp)),goodfit_ind_size_all),RFcutoff);
    n = length(inds);
    
    figure;
    start = 1;
    ifig = 1;
    for i=1:n
        iCell = inds(i);
        if start ==37
            suptitle(['Exp: ' num2str(iExp) ', area: ' char(expt(iExp).img_loc{1}) ', mouse: ' char(expt(iExp).mouse) ', date: ' char(expt(iExp).date)])
            set(gcf, 'Position', [0 0 800 1000]);
            fn_out = fullfile(fnout, ['exp' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);
            print(fn_out,'-dpdf')
            figure;
            ifig = 1+ifig;
            start = 1;
        end
        h = subplot(6,6,start);
        errorbar([0 szs],[0 sizeMean_all(:,nCon,iCell)'],[0 sizeSEM_all(:,nCon,iCell)'])
        hold on
        plot(data_all{iCell,iCon,2},data_all{iCell,iCon,1},'.k')
        plot(szRng,squeeze(fitout_all(iCell,iCon,1,:)),'-b')
        plot(szRng,squeeze(fitout_all(iCell,iCon,2,:)),'-r')
        hold off
        ylim([min([-0.5*sizeFits_all(iCell,9,iCon) min(data_all{iCell,iCon})]) 1.2*max([sizeFits_all(iCell,10,iCon) max(data_all{iCell,iCon})])])
        if sizeFits_all(iCell,11,iCon)
            title(['#' num2str(iCell), ' m1 ' num2str(cellDists_all(iCell),3)]);
        else
            title(['#' num2str(iCell), ' m2 ' num2str(cellDists_all(iCell),3)]);
        end
        
        start = start+1;
    end
    suptitle(['Exp: ' num2str(iExp) ', area: ' char(expt(iExp).img_loc{1}) ', mouse: ' char(expt(iExp).mouse) ', date: ' char(expt(iExp).date)])
    set(gcf, 'Position', [0 0 800 1000]);
    fn_out = fullfile(fnout, ['exp' num2str(iExp) '_SizeTuneFits' num2str(ifig) '.pdf']);    
    print(fn_out,'-dpdf')
end

%% compare areas for con data
% making figs 1-7
fprintf('Examine cells from each area:\n')
areas = unique(expArea);
nArea = length(areas);

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

close all
choosefig = [1:5];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
legStrs = strings(1,length(areas)); legStrs2=legStrs;
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = intersect(unique(expInd),find(cellfun(@(x) strcmp(x,areas(i)), [expt.img_loc], 'UniformOutput', 1)));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas{i}
        case 'V1'
            cutoff = 15; %v1 cutoff at 10
            excells = [631 2128];
            % case {'LM','AL'}
            %     cutoff = 15; %alm cutoff at 15
        case 'LM'
            cutoff = 15; %lm cutoff at 15
            excells = [1861 1863];
        case 'AL'
            cutoff = 15; %al cutoff at 15
            excells = [2395 1777];
        case 'PM'
            cutoff = 15; %pm cutoff at 20
            excells = [1952 2292];
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_area(i) = nExpi;
    nCells_area(i) = nCellsi;
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = ~(squeeze(sizeFits(:,8,:)));
    ism2 = squeeze(sizeFits(:,8,:));
    
    cons_c = categorical(cellstr(num2str(cons'))');
    szs_c = cellstr(num2str(chop(szs,2)'))';
    conStruct = conStruct_all(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas{i},nCellsi,nExpi);
    
    if sum(choosefig==1)
        figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
        subplot(2,nArea,i)
        modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
        bar(cons_c,modelcounts,'stacked')
        title({sprintf('Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Contrast')
        ylabel('Frac. cells')
        if i==nArea;legend('m1','m2','location','best');end
        subplot(2,nArea,nArea+i)
        histogram(cellDists_all(ind),[0:4:20])
        xlabel('Cell RF distance')
        ylabel('Number of cells')
    end
    
    if sum(choosefig==2)
        figure(2);if i==1;clf;end %figure 2 = average size tuning curves (normalized)
        % change now to normalize all cells, then plot 3 subplots of m1/m2/all
        sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
        for iCell = 1:nCellsi
            dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
            %dum = sizeMean(:,nCon,iCell); % only at highest con
            norm = max(dum(:)); % take max of all dF/F's including all cons
            sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
            sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
        end
        sizeMean_normall = mean(sizeMean_norm,3);
        norm = max(sizeMean_normall(:,nCon));
        sizeMean_normall = sizeMean_normall/norm;
        sizeSEM_normall = geomean(sizeSEM_norm,3)/norm;
        % split by model
        %subplot(4,3,3*(i-1)+1)
        %subplot(2,4,2*(i-1)+1)
        % for iCon = 1:nCon
        %     errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism1(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism1(:,iCon))),3))
        %     hold on
        % end
        % title({sprintf('Model1: Area:%s',areas{i});['(n=' num2str(mean(sum(ism1))) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Size (deg)')
        % ylabel('dF/F (norm)')
        % ylim([0 1.2])
        % subplot(2,4,2*(i-1)+2)
        % for iCon = 1:nCon
        %     errorbar(szs,mean(sizeMean_norm(:,iCon,find(ism2(:,iCon))),3),geomean(sizeSEM_norm(:,iCon,find(ism2(:,iCon))),3))
        %     hold on
        % end
        % title({sprintf('Model2: Area:%s',areas{i});['(n=' num2str(mean(sum(ism2))) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Size (deg)')
        % ylabel('dF/F (norm)')
        % ylim([0 1.2])
        % collapse models
        subplot(1,nArea,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon))
            hold on
        end
        title({sprintf('Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F (norm)')
        ylim([0 2])
        if i==nArea;legend(num2str(cons'),'location','best');end
    end
    
    if sum(choosefig==3)
        figure(3);if i==1;clf;end %figure 3 = prefSize vs con
        prefSize = squeeze(sizeFits(:,1,:));
        % subplot(2,2,i)
        % prefMean1=zeros(1,nCon);prefSEM1=prefMean1;prefMean2=prefMean1;prefSEM2=prefMean1;prefMeanAll=prefMean1;prefSEMAll=prefMean1;
        % for iCon=1:nCon
        %     prefMean1(iCon) = mean(prefSize(find(ism1(:,iCon)),iCon));
        %     prefSEM1(iCon) = std(prefSize(find(ism1(:,iCon)),iCon))./sqrt(sum(ism1(:,iCon)));
        %     prefMean2(iCon) = mean(prefSize(find(ism2(:,iCon)),iCon));
        %     prefSEM2(iCon) = std(prefSize(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        %     prefMeanAll(iCon) = mean(prefSize(:,iCon));
        %     prefSEMAll(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
        % end
        % errorbar(cons,prefMean1,prefSEM1,'s-');
        % hold on
        % errorbar(cons,prefMean2,prefSEM2,'^-');
        % errorbar(cons,prefMeanAll,prefSEMAll,'kx-');
        % hold off
        % title({sprintf('Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Contrast')
        % ylabel('PrefSize')
        % xlim([0 1])
        % ylim([0 60])
        % if i==4;legend('m1','m2','all','location','best');end
        prefMean=zeros(1,nCon);prefSEM=prefMean;
        for iCon=1:nCon
            prefMean(iCon) = mean(prefSize(:,iCon));
            prefSEM(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
        end
        errorbar(cons,prefMean,prefSEM);
        hold on
        title('Mean Preferred Size by Area')
        xlabel('Contrast')
        ylabel('PrefSize')
        xlim([0 1])
        ylim([0 100])
        if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==4)
        figure(4);if i==1;clf;end %figure 4 = suppInd vs con
        suppInd = squeeze(sizeFits(:,4,:));
        suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
        % subplot(2,2,i)
        % suppMean2=zeros(1,nCon);suppSEM2=suppMean2;suppMeanAll=suppMean2;suppSEMAll=suppMean2;
        % for iCon=1:nCon
        %     suppMean2(iCon) = mean(suppInd(find(ism2(:,iCon)),iCon));
        %     suppSEM2(iCon) = std(suppInd(find(ism2(:,iCon)),iCon))./sqrt(sum(ism2(:,iCon)));
        %     suppMeanAll(iCon) = mean(suppInd(:,iCon));
        %     suppSEMAll(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
        % end
        % errorbar(cons,suppMean2,suppSEM2,'^-');
        % hold on
        % errorbar(cons,suppMeanAll,suppSEMAll,'kx-');
        % hold off
        % title({sprintf('Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        % xlabel('Contrast')
        % ylabel('Supp Ind')
        % xlim([0 1])
        % ylim([0 1])
        % if i==4;legend('m2 only','all','location','best');end
        
        suppMean=zeros(1,nCon);suppSEM=suppMean;
        for iCon=1:nCon
            suppMean(iCon) = mean(suppInd(:,iCon));
            suppSEM(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
        end
        errorbar(cons,suppMean,suppSEM);
        hold on
        title('Mean Suppression Index by area')
        xlabel('Contrast')
        ylabel('SI')
        xlim([0 1])
        ylim([0 1])
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas{i},nCellsi,nExpi);
        if i==nArea;legend(legStrs,'location','eastoutside');end %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==5) %figure 5: average contrast response in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',areas{i},length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(5);if i==1;clf;end
        ax = gca;
        ax.ColorOrderIndex = i;
        %subplot(2,2,i)
        %for iCell = 1:nCellsi
        %    p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %    p1.Color(4) = 0.1;
        %    hold on
        %end
        hold on
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response by area')
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==nArea;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
        % fit
        conResp_norm = conResp_norm';
        cRi = conResp_norm(:);
        cons_exp = repmat(cons,1,length(cut));
        lb = [0 0 0.1 1];
        ub = [Inf Inf 0.8 Inf];
        SStot = sum((cRi-mean(cRi)).^2);
        R2best = -Inf;
        x0 = [cRi(1) mean(cRi) 0.2 3]; %BL Rmax C50 n
        [cF, res] = lsqcurvefit(conModelH,x0,cons_exp',cRi,lb,ub,opts);
        R2 = 1-res/SStot;
        
        fitout = conModelH(cF,conRng);
        R50 = fitout(1)+(fitout(end)-fitout(1))/2;
        fitout50rect = abs(fitout - R50);
        i50 = find(fitout50rect == min(fitout50rect),1);
        C50 = conRng(i50);
        
        ax = gca;
        ax.ColorOrderIndex = i;
        plot(conRng,fitout,':','HandleVisibility','off')
        ax = gca;
        ax.ColorOrderIndex = i;
        plot(C50,R50,'x','HandleVisibility','off')
        ax = gca;
        ax.ColorOrderIndex = i;
        plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    end

    if sum(choosefig==6) %figure 6: example cells from each area, with fits
        figure(6);if i==1;clf;end
        subplot(2,4,2*(i-1)+1)
        dum = sizeMean_all(:,:,excells(1)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(1))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(1))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Non-suppressed cell in %s',areas{i}))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        if sum(i==[1 3]);ylabel('dF/F (norm)');end
        ylim([0 1.2])
        subplot(2,4,2*(i-1)+2)
        dum = sizeMean_all(:,:,excells(2)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(2))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(2))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Suppressed cell in %s',areas{i}))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        %ylabel('dF/F (norm)')
        ylim([0 1.2])
        if i==4;legend(num2str(cons'));end
    end
    
    if sum(choosefig==7) % median model parameters curve
        par1 = zeros(nCellsi,nCon,3);
        par2 = zeros(nCellsi,nCon,6);
        medpar1 = zeros(nCon,3);
        medpar2 = zeros(nCon,6);
        for iCell = 1:nCellsi
            for iCon = 1:nCon
                par1(iCell,iCon,:) = sizeFits(iCell,iCon).True.s_.fit1.c1;
                par2(iCell,iCon,:) = sizeFits(iCell,iCon).True.s_.fit2.c2;
            end
        end
        for iCon = 1:nCon
            ism1_i = find(~sizeFits_all(:,8,iCon));
            ism2_i = find(sizeFits_all(:,8,iCon));
            medpar1(iCon,:) = squeeze(median(par1(ism1_i,iCon,:),1));
            medpar2(iCon,:) = squeeze(median(par2(ism2_i,iCon,:),1));
            meanpar1(iCon,:) = squeeze(mean(par1(ism1_i,iCon,:),1));
            meanpar2(iCon,:) = squeeze(mean(par2(ism2_i,iCon,:),1));
        end
        medpar1(:,1) = medpar1(:,1)./max(medpar1(nCon,1));
        medpar2(:,1) = medpar2(:,1)./max(medpar2(nCon,1));
        meanpar1(:,1) = meanpar1(:,1)./max(meanpar1(nCon,1));
        meanpar2(:,1) = meanpar2(:,1)./max(meanpar2(nCon,1));
        medpar1 = (medpar1+meanpar1)/2;
        medpar2 = (medpar2+meanpar2)/2;
        m1 = str2func(sizeFits(1,1).True.s_.m1);
        m2 = str2func(sizeFits(1,1).True.s_.m2);
        
        figure(7);if i==1;clf;end
        subplot(2,2,i)
        for iCon=1:nCon
            plot(szRng,m1(medpar1(iCon,:),szRng))
            hold on
        end
        ax = gca;
        ax.ColorOrderIndex = 1;
        for iCon=1:nCon
            plot(szRng,m2(medpar2(iCon,:),szRng))
            hold on
        end
        legend(num2str(cons'))
        title({sprintf('Area:%s',areas{i});['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        xlabel('Size (deg)')
        ylabel('dF/F')
        xlim([0 max(szs)+1])
    end
    
    if sum(choosefig==8) % contrast-resp fits
        % look at different data:
        % histogram of iC50 (best guess) vs C50fit and C50 rec
        % cross plots of these, C50rec-C50fit
        % Rsq histogram, choose cutoff, look at ratio by area
        % C50rec across areas (boxplot + mean)
        C50f = 0*ind;
        for iCell=1:length(ind)
            C50f(iCell) = conStruct(iCell).fit(3);
        end
        C50r = [conStruct.C50r];
        Rsq = [conStruct.Rsq];
        cut = find(Rsq>0.9);
        figure(8);if i==1;clf;end
        subplot(2,2,i)
        plot(C50f,C50r,'.')
        xlabel('C50f')
        ylabel('C50r')
        title({sprintf('Area:%s',areas{i});['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})
        
    end
    
end
figure(1)
fn_out = fullfile(fnout, ['allAreaModelID_' cellType '.pdf']);    
print(fn_out,'-dpdf')
figure(2)
fn_out = fullfile(fnout, ['allAreaSizeTuningByCon_' cellType '.pdf']);    
print(fn_out,'-dpdf')
figure(3)
fn_out = fullfile(fnout, ['allAreaPrefSizeByCon_' cellType '.pdf']);    
print(fn_out,'-dpdf') 
figure(4)
fn_out = fullfile(fnout, ['allAreaSuppIndexByCon_' cellType '.pdf']);    
print(fn_out,'-dpdf') 
figure(5)
fn_out = fullfile(fnout, ['allAreaContrastResp_' cellType '.pdf']);    
print(fn_out,'-dpdf') 

%% examine contrast fits - C50 and BLr
Rsqcutoff = 0.9;
% extract C50i and
C50i = 0*goodfit_ind_size_all;
C50f = 0*goodfit_ind_size_all;
for i=1:length(goodfit_ind_size_all)
    C50i(i) = conStruct_all(goodfit_ind_size_all(i)).x0(3);
    C50f(i) = conStruct_all(goodfit_ind_size_all(i)).fit(3);
end
C50r = [conStruct_all.C50r];

Rsq = [conStruct_all.Rsq];
figure(9);clf;
subplot(1,2,1)
histogram(Rsq,-0.05:0.1:1.05)
title('R^2 all cells')
cut = find(Rsq>Rsqcutoff);
subplot(1,2,2)
histogram(Rsq(cut),(Rsqcutoff-0.01):0.01:1.01)
title(['R^2 cutoff at >' num2str(Rsqcutoff,2)])

nCut = 0*nCells_area;
rCut=nCut;
x =[]; yC50=x; yBLr=x;
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas{i}) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas{i}), [expt.img_loc], 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    switch areas{i}
        case 'V1'
            cutoff = 15; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 15; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeFits = sizeFits_all(ind,:,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = ~sizeFits(:,8,:);
    ism2 = sizeFits(:,8,:);
    
    conStruct = conStruct_all(ind);
    C50f = 0*ind;
    baseline = C50f;
    Rmax = C50f;
    for iCell = 1:nCellsi
        C50f(iCell) = conStruct(iCell).fit(3);
        baseline(iCell) = conStruct(iCell).fit(1);
        Rmax(iCell) = conStruct(iCell).fit(2);
    end
    C50r = [conStruct.C50r]';
    
    Rsq = [conStruct.Rsq];
    cut = find(Rsq>Rsqcutoff);
    nCut(i) = length(cut);
    rCut(i) = nCut(i)/nCellsi;
    
    x = [x; i*ones(size(C50r(cut)))];
    yC50 = [yC50; C50r(cut)];
    C50mean(i) = mean(C50r(cut));
    C50SEM(i) = std(C50r(cut))./sqrt(length(cut));
    BLr = baseline(cut)./(baseline(cut)+Rmax(cut));
    yBLr = [yBLr; BLr];
    BLrmean(i) = mean(BLr);
    BLrSEM(i) = std(BLr)./sqrt(length(cut));
    
    figure(11);if i==1;clf;end
    subplot(2,2,i)
    histogram(C50r,0.15:0.05:0.85)
    hold on
    histogram(C50r(cut),0.15:0.05:0.85)
    title(['C50r, area:' char(areas{i})])
    xlabel('C50_r')
    if i==4;legend('all',['R^2>' num2str(Rsqcutoff,2)]);end
    
    %figure(12);if i==1;clf;end
    %subplot(2,2,i)
    %histogram(C50r(cut)./C50f(cut),0.15:0.1:1.05)
    %title(['Ratio rC_{50}/fC_{50}, area:' char(areas{i})])
    %xlabel('rC50/fC50')
    %ylabel('#cells')
    
    %figure(13);if i==1;clf;end
    %subplot(2,2,i)
    %plot(Rsq(cut),C50r(cut)./C50f(cut),'.')
    %title(['Ratio rC_{50}/fC_{50} vs Rsq, area:' char(areas{i})])
    %xlabel('R^2')
    %ylabel('rC50/fC50')
    
    figure(14);if i==1;clf;end
    subplot(2,2,i)
    histogram(baseline(cut)./(baseline(cut)+Rmax(cut)),-0.05:0.1:1.05)
    title(['BL/(BL+R_{max}), area:' char(areas{i})])
    xlabel('BL/(BL+R_{max})')
    ylabel('#cells')
    
end

figure(15);clf
c_areas = categorical({'V1' 'PM' 'all'},{'V1' 'PM' 'all'});
nTot = sum(nCut./rCut);
nCut2 = [nCut nTot];
rCut2 = [rCut sum(nCut)/nTot];
bar(c_areas,rCut2)
ylim([0 1])
xlabel('Area')
ylabel('frac. cutoff')
title('Con-fit cutoffs in each area')

figure(16);clf;
boxplot(yC50,x)
hold on
errorbar(1:4,C50mean,C50SEM,'x')
plot([1 4],[0.82 0.82],'-k', 'LineWidth',2)
plot(2.5,0.85,'*k')
hold off
set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['LM (n=' num2str(nCut(2)) ')'],['AL (n=' num2str(nCut(3)) ')'],['PM (n=' num2str(nCut(4)) ')']})
title('C_{50} by area')
ylabel('C_{50}')
legend('mean+SEM')
ylim([0 1])

figure(17);clf;
boxplot(yBLr,x)
hold on
%errorbar(1:4,BLrmean,BLrSEM,'x')
errorbar(1:2,BLrmean,BLrSEM,'x')
plot([1 4],[0.93 0.93],'-k', 'LineWidth',2)
plot(2.5,0.95,'*k')
plot([2 4],[0.88 0.88],'-k', 'LineWidth',2)
plot(3,0.9,'*k')
plot([3 4],[0.83 0.83],'-k', 'LineWidth',2)
plot([3.4 3.5 3.6],[0.85 0.85 0.85],'*k')
hold off
set(gca,'XTick',1:2,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['PM (n=' num2str(nCut(2)) ')']})
%set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(nCut(1)) ')'],['LM (n=' num2str(nCut(2)) ')'],['AL (n=' num2str(nCut(3)) ')'],['PM (n=' num2str(nCut(4)) ')']})
title('BL/(BL+R_{max}) by area')
ylabel('BL/(BL+R_{max})')
legend('mean+SEM')
ylim([-0.1 1])

%% stats on C50, BLr
[p,~,statC50] = anova1(yC50,x)
[results, means] = multcompare(statC50,'CType','hsd')

[p,~,statBLr] = anova1(yBLr,x)
[results, means] = multcompare(statBLr,'CType','hsd')
%%

for i=1:4
    for j=i+1:4
        comb = [i j];
        df = nCut(i)+nCut(j)-2;
        tC50 = (C50mean(i)-C50mean(j))./sqrt(C50SEM(i).^2+C50SEM(j).^2);
        tBLr = (BLrmean(i)-BLrmean(j))./sqrt(BLrSEM(i).^2+BLrSEM(j).^2);
        
        fprintf('\n%s v %s\n Deg. f = %d\ntC50 = %.4f\ntBLr = %.4f\ntCrit = %.4f (0.05), %.4f (0.01), %.4f (0.001)\n',areas{i},areas{j},df,tC50,tBLr,tinv(0.05,df),tinv(0.01,df),tinv(0.001,df))
    end
end