function [goodfit_ind edge_ind dist_vec] = retFitWithBootstrap(stimStruct,respMat,baseMat,Nshuf,stimPos,fnout);
    %stimStruct- structure with fields for Az and El values for all trials
    %respMat- response window FR [nCells x nTrials] 
    %baseMat- baseline window FR [nCells x nTrials] 
    %Nshuf- number of shuffles for bootstrap
    %fnout- output directory
    %stimPos - [Az El]
    
    Azs = stimStruct.stimAzimuth;
    Els = stimStruct.stimElevation;
    Az = unique(Azs);
    El = unique(Els);
    nAz = length(Az);
    nEl = length(El);
    nStim = nAz*nEl;
    
    Ind_struct = cell(1,nStim);
    start = 1;    
    for iEl = 1:nEl
        indE = find(Els == El(iEl));
        for iAz = 1:nAz
            indA = find(Azs == Az(iAz));
            ind = intersect(indE,indA);
            Ind_struct{start} = ind;
            start = start+1;
        end
    end

    Fit_struct = [];
    [AzAz, ElEl] = meshgrid(Az,El);
    grid2.AzAz = AzAz;
    grid2.ElEl = ElEl;

    dAz = median(diff(Az));
    dEl = median(diff(El));
    Az_vec00 = Az(1):(dAz/10):Az(end);
    El_vec00 = El(1):(dEl/10):El(end);
    [AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
    grid2.AzAz00 = AzAz00;
    grid2.ElEl00 = ElEl00;

    nCells = size(respMat,1);
    p_ttest = zeros(nCells,nStim);
    h_ttest = zeros(nCells,nStim);
    h_all = zeros(1,nCells);

    fprintf('Begin shuffling...\n')
    for count_shuf = 0:Nshuf
        fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
        Im_mat_USE = zeros(nCells, nStim);
        for iCond = 1:nStim
            ind_all = Ind_struct{iCond};
            if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
            else
                ind_all_1 = ind_all;
                [h_ttest(:,iCond), p_ttest(:,iCond)] = ttest(respMat(:,ind_all), baseMat(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
            end
            Im_mat_USE(:,iCond) = mean(respMat(:,ind_all_1)-baseMat(:,ind_all_1),2);
        end
        
        ifig = 1;
        start = 1;
        for iCell = 1:nCells
            if count_shuf == 0
                if sum(h_ttest(iCell,:),2) == 0
                    ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/2));
                    if length(ind_p)<2
                        ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/3));
                        if length(ind_p)<3
                            ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/4));
                            if length(ind_p)<4
                                h_all(1,iCell) = 0;
                            else
                                h_all(1,iCell) = 1;
                            end
                        else
                            h_all(1,iCell) = 1;
                        end
                    else
                        h_all(1,iCell) = 1;
                    end
                else
                    h_all(1,iCell) = 1;
                end
            end
            if count_shuf>0
                if h_all(1,iCell) == 0
                    continue
                end
            end
            a = Im_mat_USE(iCell,:);
            if max(a,[],2) > 0
                b = reshape(a',length(Az),length(El));
                data = b';
                if count_shuf == 0
                    PLOTIT_FIT = 1;
                    SAVEALLDATA = 1;
                    Fit_2Dellipse_ret_lbub
                    eval(['Fit_struct(iCell).True.s_',' = s;']);
                else
                    SAVEALLDATA = 0;
                    PLOTIT_FIT = 0;
                    Fit_2Dellipse_ret_lbub
                    eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
                end
            end
            if count_shuf == 0 && start == 64
                set(gcf, 'Position', [0 0 800 1000]);
                fn_out = fullfile(fnout, ['RFfits' num2str(ifig) '.pdf']);
                print(fn_out,'-dpdf')
            end
        end
        if count_shuf == 0  
            set(gcf, 'Position', [0 0 800 1000]);
            fn_out = fullfile(fnout, ['RFfits' num2str(ifig) '.pdf']);
            print(fn_out,'-dpdf')
        end
    end
    fprintf('Shuffling done, saving fit results\n')
    
    fn_out = fullfile(fnout,'Fit_struct.mat');
    save(fn_out, 'Fit_struct')
    
    resp_ind = find(h_all); % h_all indicates responsive cell (by t-test against baseline)
    
    fprintf('Assessing goodness of fit\n')
    if Nshuf>1
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).True)
                eval('tmp = Fit_struct(iCell).True.s_.x;');
                eval('tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_10];');
                eval('tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_10];');
                % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
                fit_true_vec(iCell,:) = tmp;
            end
        end
        
        fit_shuf_vec = NaN(nCells,10,Nshuf);
        for count_shuf = 1:Nshuf
            for iCell = 1:nCells
                if ~isempty(Fit_struct(iCell).Shuf)
                    eval('tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;');
                    eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50];');
                    eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50];');
                    eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10];');
                    eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];');
                    % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
                    fit_shuf_vec(iCell,:,count_shuf) = tmp;
                end
            end
        end
        
        Npars = size(fit_shuf_vec,2);
        lbub_fits = NaN(nCells,Npars,5);
        alpha_bound = .025;
        ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
        ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
        for iCell = 1:nCells
            for count2 = 1:Npars
                tmp = squeeze(fit_shuf_vec(iCell,count2,:));
                [i,j] = sort(tmp); % sort in order
                lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
                lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
                lbub_fits(iCell,count2,3) = mean(i); %mean
                lbub_fits(iCell,count2,5) = std(i); %stdev
            end
            lbub_fits(iCell,:,4) = fit_true_vec(iCell,:); % true (no shuffle)
        end
    end
    
    lbub_diff = lbub_fits(:,:,2)-lbub_fits(:,:,1);
    
    goodfit_ind = [];
    for iCell = 1:nCells
        if lbub_diff(iCell,4)<unique(diff(Az))*2
            if lbub_diff(iCell,5)<unique(diff(Az))*2
                goodfit_ind = [goodfit_ind iCell];
            end
        end
    end
    fprintf(['#Good cells = ' num2str(length(goodfit_ind)) ' (first pass)...\nNow checking for RFs at ret perimeter\n'])
    
    edge_ind = goodfit_ind;
    for i=1:length(goodfit_ind)
        if sum(round(lbub_fits(goodfit_ind(i),4,4))==[min(Azs) max(Azs)])
            continue
        elseif sum(round(lbub_fits(goodfit_ind(i),5,4))==[min(Els) max(Els)])
            continue
        end
        edge_ind(i) = 0;
    end
    edge_ind(edge_ind==0) = [];
    fprintf(['#Edge cells = ' num2str(length(edge_ind)) ' (final)\nSaving good fits\n'])

    dAz = fit_true_vec(:,5) - stimPos(1);
    dEl = fit_true_vec(:,4) - stimPos(2);
    dist_vec = sqrt(dAz.^2 + dEl.^2);
    
    fn_out = fullfile(fnout, 'lbub_fits.mat');
    save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'edge_ind', 'resp_ind','dist_vec')
end