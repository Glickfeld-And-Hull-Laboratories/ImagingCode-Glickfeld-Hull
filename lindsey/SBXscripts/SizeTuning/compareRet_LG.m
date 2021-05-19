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
ds = 'retAndSzTuning_ExptList';
%cellType = 'SOM';
rc = behavConstsAV;
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
lg_fn = fullfile(fn_base, 'home\lindsey');
fndata = fullfile(lg_fn, 'Analysis\2P');
fnout = fullfile(fndata,'SizeTuning','RetSummary');
eval(ds)
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);
nExp = size(expt,2)-2;
fprintf(['Retinotopy comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading retinotopy data\n')
lbub_fits_all = [];
goodfit_ind_all = [];
lbub_fits_all_ret2 = [];
goodfit_ind_all_ret2 = [];

t = 0;
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expt(i).date;
    mouse = expt(i).mouse;
    RetImgFolder = cell2mat(expt(i).retFolder);
    ret_str = ['runs-' RetImgFolder];
    datemouse = [date '_' mouse];
    datemouseretrun = [datemouse '_' ret_str];
    filename = fullfile(fndata, datemouse, datemouseretrun, [datemouseretrun '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[datemouseretrun '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')


    nCellsExp(i) = size(lbub_fits,1);
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_all = [goodfit_ind_all tempinds];
    
    if ~isempty(expt(i).ret2Folder)
        RetImgFolder = cell2mat(expt(i).ret2Folder);
        ret_str = ['runs-' RetImgFolder];
        datemouse = [date '_' mouse];
        datemouseretrun = [datemouse '_' ret_str];
        filename = fullfile(fndata, datemouse, datemouseretrun, [datemouseretrun '_lbub_fits.mat']);
        if ~exist(filename, 'file')
            fprintf([[datemouseretrun '_lbub_fits.mat'] ' not found! Please remove from list\n'])
        end
        load(filename, 'lbub_fits', 'goodfit_ind')
        lbub_fits_all_ret2 = cat(1,lbub_fits_all_ret2,lbub_fits);
        tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
        goodfit_ind_all_ret2 = [goodfit_ind_all_ret2 tempinds];
    else
        lbub_fits_all_ret2 = cat(1,lbub_fits_all_ret2,nan(size(lbub_fits)));
    end
    fprintf('done\n')
end

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit)\n',nCellsTot,length(goodfit_ind_all))
expInd = [];
expArea = [];
expGenotype = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
    expArea = [expArea; repmat(expt(i).img_loc, [nCellsExp(i),1])];
    expGenotype = [expGenotype; repmat(expt(i).driver, [nCellsExp(i),1])];
end


%% compare genotypes and areas
areas = unique(expArea);
nArea = length(areas);

genotypes = unique(expGenotype);
nGenotype = length(genotypes);

nExp = zeros(nArea,nGenotype);
nCells = zeros(nArea,nGenotype);

close all
legStrs = strings(nArea,nGenotype); legStrs2=legStrs;
means = nan(nArea,nGenotype);
for ii = 1:nGenotype
    fprintf(['Genotype #' num2str(ii) ' : ' char(genotypes(ii)) '\n'])
    expIndg = find(cellfun(@(x) strcmp(x,genotypes(ii)), [expt.driver], 'UniformOutput', 1));
    indg = intersect(find(ismember(expInd,expIndg)),goodfit_ind_all);
    x = [];
    y_geo = [];
    y_mean = [];
    y_std = [];
    for i = 1:nArea
        fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
        expInda = find(cellfun(@(x) strcmp(x,areas(i)), [expt.img_loc], 'UniformOutput', 1));
        inda = intersect(find(ismember(expInd,expInda)),goodfit_ind_all);

        ind = intersect(inda,indg);

        nExpi = length(intersect(expIndg,expInda));
        nCellsi = length(ind);
        nExp(i,ii) = nExpi;
        nCells(i,ii) = nCellsi;
        if length(ind)>0
            lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
            sigmax = lbub_fits(:,2,4);
            sigmay = lbub_fits(:,3,4);
            RFsize = 2*sqrt(2*log(2))*geomean([sigmax sigmay],2);
            RFrads_geo{i,ii} = RFsize;
            means(i,ii) = mean(RFsize);
            x = [x; i*ones(size(RFsize))];
            y_geo = [y_geo; RFsize];
            y_mean = [y_mean mean(RFrads_geo{i})];
            y_std = [y_std std(RFrads_geo{i})];

            [fi xi] = ksdensity(RFrads_geo{i,ii});
            fnorm(:,i,ii) = fi/max(fi)*0.3;
            xinorm(:,i,ii) = xi;

            legStrs(i,ii)=sprintf('%s (n-cells =%d, n-exp=%d)',cell2mat(areas(i)),nCells(i,ii),nExp(i,ii));

            figure(1);if ii==1 & i==1;clf;end %figure 1 = proportions of model2 by con
            subplot(2, nGenotype, ii)
            scatter(i.*ones(size(RFsize)), RFsize, 'ok')
            hold on
            errorbar(i,mean(RFsize,1), std(RFsize,[],1)./sqrt(nCellsi),'or');
            title({sprintf('Genotype:%s',cell2mat(genotypes(ii)))})
            ylabel('FWHM (deg)')
            ylim([0 60])
            xlim([0 nArea+1])
            set(gca,'Xtick',1:nArea,'XTickLabel', areas)
            subplot(2, nArea, nArea+i)
            scatter(ii.*ones(size(RFsize)), RFsize, 'ok')
            hold on
            errorbar(ii,mean(RFsize,1), std(RFsize,[],1)./sqrt(nCellsi),'or');
            title({sprintf('Area:%s',cell2mat(areas(i)))})
            ylabel('FWHM (deg)')
            ylim([0 60])
            xlim([0 nGenotype+1])
            set(gca,'Xtick',1:nGenotype,'XTickLabel', genotypes)
        end
    end
    expUse = find(nExp(:,ii));
    figure(1)
    subplot(2, nGenotype, ii)
    legend(legStrs(expUse,ii))
    
    if length(expUse)>1
        [p,~,stat_geo] = anova1(y_geo,x,'Display','off');
        [results, av] = multcompare(stat_geo,'CType','hsd');
    end
    figure(4)
    for i = expUse'
        ax = subplot(2,nGenotype, ii);
        colors = get(ax,'ColorOrder');
        hold on
        h5(i)=fill([xinorm(:,i,ii);flipud(xinorm(:,i,ii))],[fnorm(:,i,ii)+(nArea+1-i);flipud((nArea+1-i)-fnorm(:,i,ii))],[1 1 1],'EdgeColor','k');
        p(1)=plot([means(i,ii) means(i,ii)],[interp1(xinorm(:,i,ii),fnorm(:,i,ii)+(nArea+1-i),means(i,ii)), interp1(flipud(xinorm(:,i,ii)),flipud((nArea+1-i)-fnorm(:,i,ii)),means(i,ii)) ],'k','LineWidth',2);
        h5(i).FaceColor = colors(i,:);
    end
    figure(4)
    axis([0 70 0.5 length(areas)+0.5]);
    legend off
    ax.YTick = [1:nArea];
    ax.YTickLabel = flipud(areas);
    set(gca,'box','off','TickDir','out')
    ylabel('Area')
    xlabel('RF size (deg)')
    title({sprintf('Genotype:%s',cell2mat(genotypes(ii)))})
    figure(4)
    for i = expUse'
        ax = subplot(2,nArea,nArea+i);
        colors = get(ax,'ColorOrder');
        hold on
        h5(i)=fill([xinorm(:,i,ii);flipud(xinorm(:,i,ii))],[fnorm(:,i,ii)+(nGenotype+1-ii);flipud((nGenotype+1-ii)-fnorm(:,i,ii))],[1 1 1],'EdgeColor','k');
        p(1)=plot([means(i,ii) means(i,ii)],[interp1(xinorm(:,i,ii),fnorm(:,i,ii)+(nGenotype+1-ii),means(i,ii)), interp1(flipud(xinorm(:,i,ii)),flipud((nGenotype+1-ii)-fnorm(:,i,ii)),means(i,ii)) ],'k','LineWidth',2);
        h5(i).FaceColor = colors(i,:);
        figure(4)
        axis([0 70 0.5 length(genotypes)+0.5]);
        legend off
        ax.YTick = [1:nGenotype];
        ax.YTickLabel = flipud(genotypes);
        set(gca,'box','off','TickDir','out')
        ylabel('Genotype')
        xlabel('RF size (deg)')
        title({sprintf('Area:%s',cell2mat(areas(i)))})
    end
end
figure(1)
print(fullfile(fnout, 'RFsize_scatter.pdf'),'-dpdf','-bestfit')
figure(4)
print(fullfile(fnout, 'RFsize_violin.pdf'),'-dpdf','-bestfit')



%%
RFsize_20_avg = nan(nGenotype,nArea);
RFsize_10_avg = nan(nGenotype,nArea);
figure(6)
start = 1;
for i = 1:nArea
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    expInda = find(cellfun(@(x) strcmp(x,areas(i)), [expt.img_loc], 'UniformOutput', 1));
    for ii = 1:nGenotype
        fprintf(['Genotype #' num2str(ii) ' : ' char(genotypes(ii)) '\n'])
        expIndg = find(cellfun(@(x) strcmp(x,genotypes(ii)), [expt.driver], 'UniformOutput', 1));
        expIndga = intersect(expIndg,expInda);
        goodfit_both = intersect(find(ismember(expInd,expIndga)),intersect(goodfit_ind_all, goodfit_ind_all_ret2));


        lbub_fits = lbub_fits_all(goodfit_both,:,:); %cell,par,val (low up mean true stdev)
        sigmax = lbub_fits(:,2,4);
        sigmay = lbub_fits(:,3,4);
        RFsize_20 = 2*sqrt(2*log(2))*geomean([sigmax sigmay],2);
        RFsize_20_avg(ii,i) = mean(RFsize_20);
        lbub_fits = lbub_fits_all_ret2(goodfit_both,:,:); %cell,par,val (low up mean true stdev)
        sigmax = lbub_fits(:,2,4);
        sigmay = lbub_fits(:,3,4);
        RFsize_10 = 2*sqrt(2*log(2))*geomean([sigmax sigmay],2);
        RFsize_10_avg(ii,i) = mean(RFsize_10);

        subplot(nArea,nGenotype,start) 
        plot(repmat([10 20], [length(goodfit_both) 1])', [RFsize_10 RFsize_20]','k')
        hold on
        errorbar([10 20], [mean(RFsize_10) mean(RFsize_20)], [std(RFsize_10)./sqrt(length(goodfit_both)) std(RFsize_20)./sqrt(length(goodfit_both))],'-or')
        xlim([5 25])
        xlabel('Stim diameter')
        ylabel('RF FWHM')
        ylim([0 40])
        title([char(genotypes(ii)) ' ' char(areas(i)) ' - n = ' num2str(length(goodfit_both))])
        start= start+1;
    end
end
figure(5)
print(fullfile(fnout, 'RFsize_byStimSize.pdf'),'-dpdf','-bestfit')