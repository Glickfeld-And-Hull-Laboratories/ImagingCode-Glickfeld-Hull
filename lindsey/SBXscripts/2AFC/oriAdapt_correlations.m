
naOri = length(aOris);
ntOri = length(tOris);
r_adapt = zeros(nCells,nCells,4,naOri);
r_stim = zeros(nCells,nCells,length(tOris),naOri+1);
geom_adapt = zeros(nCells,nCells,4,naOri);
geom_stim = zeros(nCells,nCells,length(tOris),naOri+1);
adapt_trN = zeros(1,naOri);
stim_trN = zeros(3,ntOri);
mask_label_mat = zeros(nCells);
for iOri = 1:naOri
    ind = intersect(find(aGratingContrast),find(aGratingOri==aOris(iOri)));
    adapt_trN(iOri) = length(ind);
    for iCell = 1:nCells
        for jCell = 1:nCells
            if iOri == 1
                mask_label_mat(iCell,jCell) = mask_label(iCell)+mask_label(jCell);
            end
            for i = 1:4
                r_adapt(iCell,jCell,i,iOri) = triu2vec(corrcoef(adapt_cyc_resp(iCell,ind,i),adapt_cyc_resp(jCell,ind,i)));
                imean = mean(adapt_cyc_resp(iCell,ind,i),2);
                jmean = mean(adapt_cyc_resp(jCell,ind,i),2);
                if imean<0
                    imean = 0;
                end
                if jmean<0
                    jmean = 0;
                end
                geom_adapt(iCell,jCell,i,iOri) = geomean([imean jmean]);
            end
            for tOri = 1:ntOri
                ind_ori = intersect(ind,find(tGratingOri == tOris(tOri)));
                stim_trN(iOri,tOri) = length(ind_ori);
                r_stim(iCell,jCell,tOri,iOri) = triu2vec(corrcoef(stim_resp(iCell,ind_ori),stim_resp(jCell,ind_ori)));
                imean = mean(stim_resp(iCell,ind_ori),2);
                jmean = mean(stim_resp(jCell,ind_ori),2);
                if imean<0
                    imean = 0;
                end
                if jmean<0
                    jmean = 0;
                end
                geom_stim(iCell,jCell,tOri,iOri) = geomean([imean jmean]);
                if iOri == 1
                    ind_ori = intersect(find(aGratingContrast==0),find(tGratingOri == tOris(tOri)));
                    stim_trN(end,tOri) = length(ind_ori);
                    r_stim(iCell,jCell,tOri,end) = triu2vec(corrcoef(stim_resp(iCell,ind_ori),stim_resp(jCell,ind_ori)));
                    imean = mean(stim_resp(iCell,ind_ori),2);
                    jmean = mean(stim_resp(jCell,ind_ori),2);
                    if imean<0
                        imean = 0;
                    end
                    if jmean<0
                        jmean = 0;
                    end
                    geom_stim(iCell,jCell,tOri,end) = geomean([imean jmean]);
                end
            end
        end
    end
end

edges = [0:0.02:0.12];
geom_adapt_bin = zeros(naOri,4,length(edges)-1,2);
r_adapt_bin = zeros(naOri,4,length(edges)-1,2);
start = 1;
figure;
for iOri = 1:naOri
    ind = intersect(adapt_resp_ind{iOri},find(mask_label==0));
    for i = 1:4
        geom_temp = triu2vec(geom_adapt(ind,ind,i,iOri));
        r_temp = triu2vec(r_adapt(ind,ind,i,iOri));
        [n edges bin] = histcounts(geom_temp,edges);
        for ii = 1:length(n)
            ind_bin = find(bin == ii);
            geom_adapt_bin(iOri,i,ii,1) = mean(geom_temp(ind_bin));
            geom_adapt_bin(iOri,i,ii,2) = std(geom_temp(ind_bin))./sqrt(length(ind_bin));
            r_adapt_bin(iOri,i,ii,1) = mean(r_temp(ind_bin));
            r_adapt_bin(iOri,i,ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
        end
        subplot(naOri,4,start)
        errorbar(squeeze(geom_adapt_bin(iOri,i,:,1)),squeeze(r_adapt_bin(iOri,i,:,1)),squeeze(r_adapt_bin(iOri,i,:,2)),squeeze(r_adapt_bin(iOri,i,:,2)),squeeze(geom_adapt_bin(iOri,i,:,2)),squeeze(geom_adapt_bin(iOri,i,:,2)))
        xlim([0 0.12])
        xlabel('Geom. mean dF/F')
        ylim([-0.2 0.5])
        ylabel('Correlation')
        title(['Adapt ' num2str(i)])
        hold on
        start = start+1;
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' All pyr pairs'])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_corrByGeomean_Adapt.pdf']),'-dpdf','-fillpage')

figure;
for i = 1:4
    errorbar(squeeze(geom_adapt_bin(iOri,i,:,1)),squeeze(r_adapt_bin(iOri,i,:,1)),squeeze(r_adapt_bin(iOri,i,:,2)),squeeze(r_adapt_bin(iOri,i,:,2)),squeeze(geom_adapt_bin(iOri,i,:,2)),squeeze(geom_adapt_bin(iOri,i,:,2)))
    xlim([0 0.12])
    xlabel('Geom. mean dF/F')
    ylim([-0.2 0.5])
    ylabel('Correlation')
    hold on
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' All pyr pairs'])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_corrByGeomean_allAdapt.pdf']),'-dpdf','-fillpage')

edges = [0:0.02:0.12];
geom_stim_bin = zeros(naOri,ntOri,length(edges)-1,2);
r_stim_bin = zeros(naOri,ntOri,length(edges)-1,2);
start = 1;
figure;
for iOri = 1:naOri+1
    ind = intersect(stim_resp_ind,find(mask_label==0));
    for i = 1:ntOri
        geom_temp = triu2vec(geom_stim(ind,ind,i,iOri));
        r_temp = triu2vec(r_stim(ind,ind,i,iOri));
        [n edges bin] = histcounts(geom_temp,edges);
        for ii = 1:length(n)
            ind_bin = find(bin == ii);
            geom_stim_bin(iOri,i,ii,1) = mean(geom_temp(ind_bin));
            geom_stim_bin(iOri,i,ii,2) = std(geom_temp(ind_bin))./sqrt(length(ind_bin));
            r_stim_bin(iOri,i,ii,1) = mean(r_temp(ind_bin));
            r_stim_bin(iOri,i,ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
        end
        subplot(naOri+1,ntOri,start)
        errorbar(squeeze(geom_stim_bin(iOri,i,:,1)),squeeze(r_stim_bin(iOri,i,:,1)),squeeze(r_stim_bin(iOri,i,:,2)),squeeze(r_stim_bin(iOri,i,:,2)),squeeze(geom_stim_bin(iOri,i,:,2)),squeeze(geom_stim_bin(iOri,i,:,2)))
        xlim([0 0.12])
        xlabel('Geom. mean dF/F')
        ylim([-0.2 0.5])
        ylabel('Correlation')
        if iOri == naOri+1
            title(['Control: ' num2str(tOris(i))])
        else
            title(['Adapt ' num2str(aOris(iOri)) ': '  num2str(tOris(i))])
        end
        hold on
        start = start+1;
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' All pyr pairs'])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_corrByGeomean_Target.pdf']),'-dpdf','-fillpage')


figure;
cell_combo = [0 0 1; 0 1 1];
subplot(1,2,1)
for ic = 1:length(cell_combo)
    geom_temp = [];
    r_temp = [];
    for iOri = 1:naOri
        ind1 = intersect(adapt_resp_ind{iOri},find(mask_label==cell_combo(1,ic)));
        ind2 = intersect(adapt_resp_ind{iOri},find(mask_label==cell_combo(2,ic)));
        for i = 1:4
            if ic ~=2
                geom_temp = [geom_temp; triu2vec(geom_adapt(ind1,ind2,i,iOri))];
                r_temp = [r_temp; triu2vec(r_adapt(ind1,ind2,i,iOri))];
            else
                geom_temp = [geom_temp; reshape(geom_adapt(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
                r_temp = [r_temp; reshape(r_adapt(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
            end
        end
    end
    edges = [0:0.02:0.12];
    geom_adapt_bin = zeros(length(edges)-1,2);
    r_adapt_bin = zeros(length(edges)-1,2);
    [n edges bin] = histcounts(geom_temp,edges);
    for ii = 1:length(n)
        ind_bin = find(bin == ii);
        geom_adapt_bin(ii,1) = mean(geom_temp(ind_bin));
        geom_adapt_bin(ii,2) = std(geom_temp(ind_bin))./sqrt(length(ind_bin));
        r_adapt_bin(ii,1) = mean(r_temp(ind_bin));
        r_adapt_bin(ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
    end

    errorbar(squeeze(geom_adapt_bin(:,1)),squeeze(r_adapt_bin(:,1)),squeeze(r_adapt_bin(:,2)),squeeze(r_adapt_bin(:,2)),squeeze(geom_adapt_bin(:,2)),squeeze(geom_adapt_bin(:,2)))
    hold on
end
xlim([0 0.12])
ylim([-0.5 0.5])
cell_combo_name = cell(size(cell_combo,2),1);
for i = 1:3
    for ii = 1:2
        if cell_combo(ii,i)
            str = 'PV';
        else
            str = 'Py';
        end
        if ii == 1
            str = [str ' - '];
        end
        cell_combo_name{i} = [cell_combo_name{i} str];
    end
end
legend(cell_combo_name','location','southeast')
title('Adapt')


subplot(1,2,2)
for ic = 1:length(cell_combo)
    geom_temp = [];
    r_temp = [];
    for iOri = 1:naOri+1
        ind1 = intersect(stim_resp_ind,find(mask_label==cell_combo(1,ic)));
        ind2 = intersect(stim_resp_ind,find(mask_label==cell_combo(2,ic)));
        for i = 1:ntOri
            if ic ~=2
                geom_temp = [geom_temp; triu2vec(geom_stim(ind1,ind2,i,iOri))];
                r_temp = [r_temp; triu2vec(r_stim(ind1,ind2,i,iOri))];
            else
                geom_temp = [geom_temp; reshape(geom_stim(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
                r_temp = [r_temp; reshape(r_stim(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
            end
        end
    end
    edges = [0:0.02:0.12];
    geom_stim_bin = zeros(length(edges)-1,2);
    r_stim_bin = zeros(length(edges)-1,2);
    [n edges bin] = histcounts(geom_temp,edges);
    for ii = 1:length(n)
        ind_bin = find(bin == ii);
        geom_stim_bin(ii,1) = mean(geom_temp(ind_bin));
        geom_stim_bin(ii,2) = std(geom_temp(ind_bin))./sqrt(length(ind_bin));
        r_stim_bin(ii,1) = mean(r_temp(ind_bin));
        r_stim_bin(ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
    end

    errorbar(squeeze(geom_stim_bin(:,1)),squeeze(r_stim_bin(:,1)),squeeze(r_stim_bin(:,2)),squeeze(r_stim_bin(:,2)),squeeze(geom_stim_bin(:,2)),squeeze(geom_stim_bin(:,2)))
    hold on
end
    xlim([0 0.12])
    ylim([-0.5 0.5])
    title('Target')
suptitle([expt(iexp).mouse ' ' expt(iexp).date])
print(fullfile([LG_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_corrByGeomean_AllCombo.pdf']),'-dpdf','-fillpage')
    

    
sig_corr = zeros(nCells);
pref_diff = zeros(nCells);
[pref_val pref_ori] = max(squeeze(vonMisesFitAllCellsAllBoots(:,1,:)),[],1);
for iCell = 1:nCells
    for jCell = 1:nCells
        sig_corr(iCell,jCell) = triu2vec(corrcoef(avgResponseEaOri(iCell,:),avgResponseEaOri(jCell,:)));
        pref_diff(iCell,jCell) = abs(pref_ori(iCell)-pref_ori(jCell));
        if pref_diff(iCell,jCell)>90
            pref_diff(iCell,jCell) = 180-pref_diff(iCell,jCell);
        end
    end
end
    
figure;
subplot(1,2,1)
for ic = 1:length(cell_combo)
    sig_temp = [];
    r_temp = [];
    for iOri = 1:naOri
        ind1 = intersect(adapt_resp_ind{iOri},find(mask_label==cell_combo(1,ic)));
        ind2 = intersect(adapt_resp_ind{iOri},find(mask_label==cell_combo(2,ic)));
        for i = 1:4
            if ic ~=2
                sig_temp = [sig_temp; triu2vec(sig_corr(ind1,ind2))];
                r_temp = [r_temp; triu2vec(r_adapt(ind1,ind2,i,iOri))];
            else
                sig_temp = [sig_temp; reshape(sig_corr(ind1,ind2), [length(ind1).*length(ind2) 1])];
                r_temp = [r_temp; reshape(r_adapt(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
            end
        end
    end
          
    edges = [-1:0.2:1];
    sig_bin = zeros(length(edges)-1,2);
    r_adapt_bin = zeros(length(edges)-1,2);
    [n edges bin] = histcounts(sig_temp,edges);
    for ii = 1:length(n)
        ind_bin = find(bin == ii);
        sig_bin(ii,1) = mean(sig_temp(ind_bin));
        sig_bin(ii,2) = std(sig_temp(ind_bin))./sqrt(length(ind_bin));
        r_adapt_bin(ii,1) = mean(r_temp(ind_bin));
        r_adapt_bin(ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
    end
    errorbar(squeeze(sig_bin(:,1)),squeeze(r_adapt_bin(:,1)),squeeze(r_adapt_bin(:,2)),squeeze(r_adapt_bin(:,2)),squeeze(sig_bin(:,2)),squeeze(sig_bin(:,2)))
    hold on
end
    xlim([-1 1])
    xlabel('Signal correlation')
    ylim([-0.5 0.5])
    ylabel('Noise correlation')
    title('Adapt')
    legend(cell_combo_name','location','southeast')
    subplot(1,2,2)
for ic = 1:length(cell_combo)
    sig_temp = [];
    r_temp = [];
    for iOri = 1:naOri+1
        ind1 = intersect(stim_resp_ind,find(mask_label==cell_combo(1,ic)));
        ind2 = intersect(stim_resp_ind,find(mask_label==cell_combo(2,ic)));
        for i = [1:ntOri]
            if ic ~=2
                sig_temp = [sig_temp; triu2vec(sig_corr(ind1,ind2))];
                r_temp = [r_temp; triu2vec(r_stim(ind1,ind2,i,iOri))];
            else
                sig_temp = [sig_temp; reshape(sig_corr(ind1,ind2), [length(ind1).*length(ind2) 1])];
                r_temp = [r_temp; reshape(r_stim(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
            end
        end
    end
          
    edges = [-1:0.2:1];
    sig_bin = zeros(length(edges)-1,2);
    r_stim_bin = zeros(length(edges)-1,2);
    [n edges bin] = histcounts(sig_temp,edges);
    for ii = 1:length(n)
        ind_bin = find(bin == ii);
        sig_bin(ii,1) = mean(sig_temp(ind_bin));
        sig_bin(ii,2) = std(sig_temp(ind_bin))./sqrt(length(ind_bin));
        r_stim_bin(ii,1) = mean(r_temp(ind_bin));
        r_stim_bin(ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
    end
    errorbar(squeeze(sig_bin(:,1)),squeeze(r_stim_bin(:,1)),squeeze(r_stim_bin(:,2)),squeeze(r_stim_bin(:,2)),squeeze(sig_bin(:,2)),squeeze(sig_bin(:,2)))
    hold on
end
    xlim([-1 1])
    xlabel('Signal correlation')
    ylim([-0.5 0.5])
    ylabel('Noise correlation')
    title('Target')
    
  figure;
cell_combo = [0 0 1; 0 1 1];
for ic = 1:length(cell_combo)
    pref_temp = [];
    r_temp = [];
    for iOri = 1:naOri+1
        ind1 = intersect(stim_resp_ind,find(mask_label==cell_combo(1,ic)));
        ind2 = intersect(stim_resp_ind,find(mask_label==cell_combo(2,ic)));
        for i = [1:6]
            if ic ~=2
                pref_temp = [pref_temp; triu2vec(pref_diff(ind1,ind2))];
                r_temp = [r_temp; triu2vec(r_stim(ind1,ind2,i,iOri))];
            else
                pref_temp = [pref_temp; reshape(pref_diff(ind1,ind2), [length(ind1).*length(ind2) 1])];
                r_temp = [r_temp; reshape(r_stim(ind1,ind2,i,iOri), [length(ind1).*length(ind2) 1])];
            end
        end
    end
          
    edges = [0:15:90];
    pref_bin = zeros(length(edges)-1,2);
    r_stim_bin = zeros(length(edges)-1,2);
    [n edges bin] = histcounts(pref_temp,edges);
    for ii = 1:length(n)
        ind_bin = find(bin == ii);
        pref_bin(ii,1) = mean(pref_temp(ind_bin));
        pref_bin(ii,2) = std(pref_temp(ind_bin))./sqrt(length(ind_bin));
        r_stim_bin(ii,1) = mean(r_temp(ind_bin));
        r_stim_bin(ii,2) = std(r_temp(ind_bin))./sqrt(length(ind_bin));
    end
    errorbar(squeeze(pref_bin(:,1)),squeeze(r_stim_bin(:,1)),squeeze(r_stim_bin(:,2)),squeeze(r_stim_bin(:,2)),squeeze(pref_bin(:,2)),squeeze(pref_bin(:,2)))
    hold on
end
    xlim([0 90])
    ylim([-0.5 0.5])
    legend(num2str(cell_combo'))
    title('Target')    