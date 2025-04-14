%% Analyze raw marmoset data sent from Nicholas Priebe, 9/23/2024
clc; clear all; close all;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
svName = 'randDirFourPhase_CrossOri';
bnsz = 10; %spike binsize = 10 ms

doPlot=0;
doGratPlot=0;

expts = strvcat('g01', 'g06', 'g12', 'g17'); %Four acute penetrations in V1, 1 marmoset
nexp = length(expts);

iexp = 4;

    fprintf([expts(iexp,:) '\n'])

    % Load data from experiment. Each experiment file contains one
    % variable, resp, that is nCells by nDir by nPhase by Grating/Plaid by
    % nTrials by Time (in 10 ms bins). For example, 310x12x4x2x13x120.
    % Grating type only has the first phase filled in.
    load(fullfile(base, 'Data\fromNicholas\CrossOri_randDirFourPhase_V1_marmoset', ([expts(iexp,:) '.mat'])))
    
    nCells = size(resp,1);
    nDirs = size(resp,2);
    nPhas = size(resp,3);
    nTrials = size(resp,5);
    nStim = (size(resp,3)+1)*size(resp,2);

    resp_cell = resp(:,:,:,:,:,21:120);
    base_cell = resp(:,:,:,:,:,1:20);
    avg_resp_dir(:,:,:,:,1) = mean(sum(resp_cell,6)/100,5); % Average across time and trial
    avg_resp_dir(:,:,:,:,2) = std(sum(resp_cell,6)/100,0,5)/nTrials; % SEM for each stim condition
    resp_dir_tc = mean(resp,5); % Average across trial, to look at time 
    resp_dir_tr = mean(resp,6); % Average across time, to look at trial 


    % Find cells that are significantly responsive to 1 stimulus
    resp_cell_trials(:,:,:,:,:,1) = sum(resp_cell,6)/100; % convert to one spike rate per trial (in Hz)
    base_cell_trials(:,:,:,:,:,1) = (sum(resp_cell,6)*5)/100; % multiply by 5 to convert to Hz
    
    for id = 1:nDirs
        [h_resp(:,id,1,1), p_resp(:,id,1,1)] = ttest2(squeeze(resp_cell_trials(:,id,1,1,:)),squeeze(base_cell_trials(:,id,1,1,:)),'dim',2,'tail','right','alpha', 0.05./nStim);
        for ip = 1:nPhas
            [h_resp(:,id,ip,2), p_resp(:,id,ip,2)] = ttest2(squeeze(resp_cell_trials(:,id,ip,2,:)),squeeze(base_cell_trials(:,id,ip,2,:)),'dim',2,'tail','right','alpha', 0.05./nStim);
        end
    end
    resp_ind_dir = find(sum(h_resp(:,:,1,1),2)); %sig responsive to gratings
    

    % figure; %trying to figure out if i have stim on + stim off or just stim on
    %     for ic = 1:36
    %         subplot(6,6,ic)
    %         plot(0:119,squeeze(resp_dir_tc(ic,1,1,1,:)))
    %     end


    % Run ANOVA across stim direction
    grat_resp = squeeze(resp_dir_tr(:,:,1,1,:)); % Make matrix nCells x nDir x nTrials
    plaid_resp = squeeze(resp_dir_tr(:,:,:,2,:)); % Make matrix nCells x nDir x nPhas x nTrials
    plaid_resp = permute(plaid_resp,[1,3,2,4]);

    group_mat = zeros(nDirs,nTrials);
    for id = 1:nDirs
        group_mat(id,:) = id;
    end

    all_grat_resp = reshape(grat_resp,[nCells,nDirs*nTrials]);
    all_plaid_resp = reshape(plaid_resp,[nCells,nPhas,nDirs*nTrials]);
    group = reshape(group_mat,[1,nDirs*nTrials]);

    p_anova_dir = zeros(1,nCells);
    for iCell = 1:nCells
        p_anova_dir(iCell) = anova1(all_grat_resp(iCell,:), group, 'off'); %direction selective to gratings
    end
    p_dir = find(p_anova_dir<0.05);

    % Determine direction selectivity
    for iCell = 1:nCells
        [max_val max_ind] = max(avg_resp_dir(iCell,:,1,1,1));
        null_ind = max_ind+(nDirs./2);
        null_ind(find(null_ind>nDirs)) = null_ind(find(null_ind>nDirs))-nDirs;
        min_val = avg_resp_dir(iCell,null_ind,1,1,1);
        if min_val < 0; min_val = 0; end
        DSI(iCell) = (max_val-min_val)./(max_val+min_val);
        DSI_maxInd(iCell) = max_ind; 
    end

    DSI_ind = find(DSI>0.5); %direction selective to gratings
    OSI_ind = find(DSI<0.5);


    % Determine pattern and component direction selectivity
    int = nDirs;
    % component = avg_resp_dir(:,:,1,1,1)+circshift(avg_resp_dir(:,:,1,1,1),-120./int,2);
    % pattern = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);

    component = circshift(avg_resp_dir(:,:,1,1,1),+60./int,2)+circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);
    pattern = avg_resp_dir(:,:,1,1,1);
    
    comp_corr = zeros(nPhas,nCells);
    patt_corr = zeros(nPhas,nCells);
    comp_patt_corr = zeros(nPhas,nCells);
    plaid_corr = zeros(1,nCells);
    plaid_corr1 = zeros(1,nCells);
    plaid_corr2 = zeros(1,nCells);
    plaid_corr3 = zeros(1,nCells);
    plaid_corr4 = zeros(1,nCells);
    plaid_corr5 = zeros(1,nCells);
    plaid_corr6 = zeros(1,nCells);
    
    for iCell = 1:nCells
        for ip = 1:nPhas
            comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
            patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
            comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
        end
        plaid_corr1(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,2,2,1)));
        plaid_corr2(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,3,2,1)));
        plaid_corr3(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr4(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,3,2,1)));
        plaid_corr5(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr6(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,3,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr(1,iCell) = (plaid_corr1(1,iCell)+plaid_corr2(1,iCell)+plaid_corr3(1,iCell)+plaid_corr4(1,iCell)+plaid_corr5(1,iCell)+plaid_corr6(1,iCell))/6;
    end
    
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nDirs-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nDirs-3));
    
    ZcZp_diff = Zc-Zp;
    ind1 = intersect(find(Zp(1,:)>1.28),find(Zp(1,:)-Zc(1,:)>1.28));
    ind2 = intersect(find(Zp(2,:)>1.28),find(Zp(2,:)-Zc(2,:)>1.28));
    ind3 = intersect(find(Zp(3,:)>1.28),find(Zp(3,:)-Zc(3,:)>1.28));
    ind4 = intersect(find(Zp(4,:)>1.28),find(Zp(4,:)-Zc(4,:)>1.28));
    ind = {ind1,ind2,ind3,ind4};
    ind_arr = [ind1,ind2,ind3,ind4];
    ind_pds = unique(ind_arr);


    % Calculate pairwise dist between Zp Zc points
    ZpZc = [];
    for iCell = 1:nCells
        ZpZc(:,1,iCell) = Zp(:,iCell);
        ZpZc(:,2,iCell) = Zc(:,iCell);
    end
    
    ZpZcPWdist = double.empty(6,0);
    for iCell = 1:nCells
        ZpZcPWdist(:,iCell) = pdist(ZpZc(:,:,iCell));
    end
    

    % Calculate PCI fit, get amplitude and baseline
    PCI = (Zp-Zc);
    phase = [0 90 180 270];
    phase_range = 0:1:359;

    for iCell = 1:nCells
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
    end


    % Direction tuning curve fit
    y_fits = zeros(360,nCells);
    stimDirs = [0:30:330];
    dirs = deg2rad(0:1:359);
    for iCell = 1:nCells
        [dir_b_hat_all(iCell,1), k1_hat_all(iCell,1), R1_hat_all(iCell,1), R2_hat_all(iCell,1), u1_hat_all(iCell,1), u2_hat_all(iCell,1), dir_sse_all(iCell,1),dir_R_square_all(iCell,1)] = miaovonmisesfit_dir(deg2rad(stimDirs),avg_resp_dir(iCell,:,1,1,1));
        dir_yfit_all(:,iCell) = b_hat_all(iCell,1)+R1_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1))-1))+R2_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1))-1));
    end


    if ~exist(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]))
        mkdir(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]))
    end
    save(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_fitsSG.mat']), 'resp_ind_dir','nCells','nTrials','nDirs','avg_resp_dir','p_dir','DSI','plaid_corr','Rp','Rc','Zp','Zc','ZpZcPWdist','yfit_all','amp_hat_all','b_hat_all','sse_all','R_square_all','dir_yfit_all','k1_hat_all','dir_sse_all','dir_R_square_all');


%% set inclusion criteria
resp_ind = intersect(resp_ind_dir,find(DSI>0.5));


if doPlot == 1
%%
figure; 
movegui('center')
for i = 1:4
    subplot(4,4,i)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{1}),Zp(i,ind{1}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 0'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+4)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{2}),Zp(i,ind{2}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 90'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+8)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{3}),Zp(i,ind{3}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 180'); end
    plotZcZpBorders
end
for i = 1:4
    subplot(4,4,i+12)
    scatter(Zc(i,resp_ind), Zp(i,resp_ind),'.')
    hold on
    scatter(Zc(i,ind{4}),Zp(i,ind{4}),'.');
    xlabel('Zc')
    ylabel('Zp')
    ylim([-4 8])
    xlim([-4 8])
    hold on
    if i==1; title('pattern cells at 270'); end
    plotZcZpBorders
end
sgtitle('Pattern direction selective cells at four phases')
print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_ZpZc_unshifted.pdf']),'-dpdf', '-fillpage') 


%% Plot grating tuning curves by cell

ind=resp_ind;

if doGratPlot == 1

respgrat = avg_resp_dir(:,:,1,1,1);
semgrat = avg_resp_dir(:,:,1,1,2);
[prefresp, prefdir] = max(respgrat,[],2);
vecshift = zeros(size(respgrat,1),1)+6; 
vecshift = vecshift - prefdir;

avg_resp_grat =[];
std_resp_grat = [];


for iCell = 1:size(respgrat,1)
    avg_resp_grat(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
    sem_resp_grat(iCell,:) = circshift(semgrat(iCell,:),vecshift(iCell),2);
end

figure;
start = 1;
n = 1;
x=[-150:30:180];
for iCell = 1:nCells
    subplot(5,4,start)
        plot(x, avg_resp_grat(iCell,:),'k') 
        shadedErrorBar(x,avg_resp_grat(iCell,:),sem_resp_grat(iCell,:))
        xlabel('grating direction')
        ylabel('sp/s')
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)]); 
        end
    start = start+1;
    if start >20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_DirectionTuningGratings_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_DirectionTuningGratings_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end
close all

else
end

%% Plot pattern index modulation by cell

PCI = (Zp-Zc);

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;
for iCell = 1:length(ind)
    ic=ind(iCell);
    subplot(5,4,start)
        scatter(phase,PCI(:,ic),8,'filled');
        hold on
        [b_hat_all(ic,1), amp_hat_all(ic,1), per_hat_all(ic,1),pha_hat_all(ic,1),sse_all(ic,1),R_square_all(ic,1)] = sinefit_PCI(deg2rad(phase),PCI(:,ic));
        yfit_all(ic,:,1) = b_hat_all(ic,1)+amp_hat_all(ic,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(ic,1) + 2.*pi/pha_hat_all(ic,1)));
        idx = ic==ind; %does ic equal any value in the index?
        if any(idx)  %if yes, plot in red. if not, plot in black
            plot(phase_range, yfit_all(ic,:,1),'k'); 
            subtitle(['cell ' num2str(ic) ', Rsq ' num2str(R_square_all(ic),'%.3f'), ', SSE ' num2str(sse_all(ic),'%.2f')],'fontweight','bold')
        else
            plot(phase_range, yfit_all(ic,:,1),'k:');
            subtitle(['cell ' num2str(ic) ', Rsq ' num2str(R_square_all(ic),'%.2f'), ', SSE ' num2str(sse_all(ic),'%.2f')])
        end
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
    start = start+1;
    if start >20
        sgtitle(['marmosetV1_' expts(iexp,:) ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['marmosetV1_' expts(iexp,:) ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')  
    end        
end
close all


%% Plot polar plots by cell

figure;
start = 1;
n = 1;

avg_resp_plddir = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);

x=[-150:30:180];
x_rad = deg2rad(x);
for iCell =1:length(ind)
    ic = ind(iCell);
    subplot(5,4,start)
        for im = 1:nPhas
            polarplot([x_rad x_rad(1)], [avg_resp_dir(ic,:,im,2,1) avg_resp_dir(ic,1,im,2,1)])
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_plddir(ic,:,1,1,1) avg_resp_plddir(ic,1,1,1,1)],'k', 'LineWidth',2) 
        idx = ic==ind;
        if any(idx)
            subtitle(['cell ' num2str(ic)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(ic)])
        end
    start = start+1;    
    if start>20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Polar plots'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PolarPlots_US_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Polar plots'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PolarPlots_US_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end     
close all

%% Plot Zp Zc classifications by cell

figure;
start = 1;
n = 1;

for iCell = 1:length(ind)
    ic = ind(iCell);
    subplot(5,4,start)
        for im = 1:nPhas
            scatter(Zc(im,ic), Zp(im,ic),8,'filled')
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        if ic ==1; legend('0 deg','90 deg','180 deg', '270 deg'); end;
        idx = ic==ind;
        if any(idx)
            subtitle(['cell ' num2str(ic)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(ic)])
        end
        plotZcZpBorders; set(gca,'TickDir','out'); axis square
    start = start+1;    
    if start>20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_US_' expts(iexp,:) '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_US_' expts(iexp,:) '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
close all


else 
end


