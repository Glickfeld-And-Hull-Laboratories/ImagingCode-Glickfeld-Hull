close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';


load(fullfile(base, (['Data\fromNicholas\lindsey_frg1.mat'])))
%% Make polar plots
x=[-150:30:180];
x_rad = deg2rad(x);

avg_resp_plaid(:,:,1) = pl90;
avg_resp_plaid(:,:,2) = pl270;
avg_resp_grat = gr;

figure(1);
for iCell =1:3
    subplot(5,4,iCell)
        for im = 1:2
            polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,im) avg_resp_plaid(iCell,1,im)])
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:) avg_resp_grat(iCell,1)],'k', 'LineWidth',2) 
end 

%% Plot Zp Zc

figure(1);
s=0;
for iCell = [5 20 24]
    s=s+1;
    subplot(5,4,s+3)
        for im = 1:4
            scatter(Zcd(im,iCell), Zpd(im,iCell))
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        subtitle(['cell ' num2str(iCell)])
        plotZcZpBorders
        set(gca,'TickDir','out'); axis square
end

%%

stimDirs = [0 30 60 90 120 150 180 210 240 270 300 330];
nStimDir = 12;
nMaskPhas = 2;
nCells = 3;

% avg_resp_dir = avg_resp_grat;
% int = unique(diff(stimDirs));
% component = avg_resp_dir(:,:)+circshift(avg_resp_dir(:,:),-120./int,2);
% pattern = circshift(avg_resp_dir(:,:),-60./int,2);
avg_resp_dir = avg_resp_grat;
int = unique(diff(stimDirs));
component = circshift(avg_resp_dir(:,:),+60./int,2)+circshift(avg_resp_dir(:,:),-60./int,2);
pattern = avg_resp_dir(:,:);

%Calculate pattern and component prediction
comp_corr = zeros(nMaskPhas,nCells);
patt_corr = zeros(nMaskPhas,nCells);
comp_patt_corr = zeros(nMaskPhas,nCells);

for iCell = 1:nCells
    for ip = 1:nMaskPhas
        comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_plaid(iCell,:,ip),component(iCell,:)));
        patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_plaid(iCell,:,ip),pattern(iCell,:)));
        comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
    end
end

Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
Zp = (0.5.*log((1+Rp)./(1-Rp))).*sqrt(1./(nStimDir-3));
Zc = (0.5.*log((1+Rc)./(1-Rc))).*sqrt(1./(nStimDir-3));

%% Fit PCI with sinusoid

PCI = (Zpd-Zcd);
nCells = size(PCI,2);

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;

for iCell = 1:nCells
    [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
    yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
end


s=0;
for iCell = [5 20 24]
    s=s+1;    
    subplot(5,4,s+6)
        scatter(phase,PCI(:,iCell),'LineWidth',1.25);
        hold on
        plot(phase_range, yfit_all(iCell,:,1),'k'); 
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-6 6])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
        subtitle([num2str(iCell)])
end

%% Plot population-- PCI by variance (amplitude)

pattern_ind = (Zpd-Zcd);    
pattind = mean(pattern_ind,1);
pattpeak = max(pattern_ind,[],1);

ind1 = intersect(find(Zpd(1,:)>1.28),find(Zpd(1,:)-Zcd(1,:)>1.28));
ind2 = intersect(find(Zpd(2,:)>1.28),find(Zpd(2,:)-Zcd(2,:)>1.28));
ind3 = intersect(find(Zpd(3,:)>1.28),find(Zpd(3,:)-Zcd(3,:)>1.28));
ind4 = intersect(find(Zpd(4,:)>1.28),find(Zpd(4,:)-Zcd(4,:)>1.28));
pattern1 = ind1';
pattern2 = ind2';
pattern3 = ind3';
pattern4 = ind4';
p = [pattern1; pattern2; pattern3; pattern4];
[C,ia,ic] = unique(p);    
a_counts = accumarray(ic,1);

p1 = C(find(a_counts==1),:);
p2 = C(find(a_counts==2),:);
p3 = C(find(a_counts==3),:);
p4 = C(find(a_counts==4),:);

c1 = [0.9375    0.7813    0.7813];
c2 = [0.9023    0.5742    0.5625];
c3 = [0.8320    0.3672    0.3398];
c4 = [0.7266    0.1094    0.1094];

figure(2);
    subplot(2,2,1)
    scatter(pattind,amp_hat_all,[],'filled')
    hold on
    scatter(pattind(p1),amp_hat_all(p1),[],c1,'filled')
    scatter(pattind(p2),amp_hat_all(p2),[],c2,'filled')
    scatter(pattind(p3),amp_hat_all(p3),[],c3,'filled')
    scatter(pattind(p4),amp_hat_all(p4),[],c4,'filled')
    xlabel('Average pattern index (Zp-Zc)')
    ylabel('Spatial variance (amp)')
    ylim([0 6])
    xlim([-4 6])
    set(gca,'TickDir','out'); axis square

    subplot(2,2,2)
    scatter(pattpeak,amp_hat_all,[],'filled')
    hold on
    scatter(pattpeak(p1),amp_hat_all(p1),[],c1,'filled')
    scatter(pattpeak(p2),amp_hat_all(p2),[],c2,'filled')
    scatter(pattpeak(p3),amp_hat_all(p3),[],c3,'filled')
    scatter(pattpeak(p4),amp_hat_all(p4),[],c4,'filled')
    xlabel('Peak pattern index (Zp-Zc)')
    ylabel('Spatial variance (amp)')
    ylim([0 6])
    xlim([-4 6])
    set(gca,'TickDir','out'); axis square

    
%% figure 3 - population Zp Zc by phase

Zp_all = Zpd;
Zc_all = Zcd;

    ind1 = intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28));
    ind2 = intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28));
    ind4 = intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28));
    pattern1 - ind1;
    pattern2 = ind2;
    pattern4 = ind4;
figure(3);    
    subplot(1,4,1)
        scatter(Zc_all(1,:),Zp_all(1,:),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(1,pattern1),Zp_all(1,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        
   subplot(1,4,2)
        scatter(Zc_all(2,:),Zp_all(2,:),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(2,pattern1),Zp_all(2,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        
   subplot(1,4,3)
        scatter(Zc_all(3,:),Zp_all(3,:),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(3,pattern1),Zp_all(3,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square
        
   subplot(1,4,4)
        scatter(Zc_all(4,:),Zp_all(4,:),15,[0.7 0.7 0.7])
        hold on
        scatter(Zc_all(4,pattern1),Zp_all(4,pattern1),15)
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square


    
    %%
    figure(1); print(fullfile(outDir, [svName '_marmoset_examplecells.pdf']),'-dpdf', '-fillpage') 
    figure(2); print(fullfile(outDir, [svName '_marmoset_population.pdf']),'-dpdf', '-fillpage') 
    figure(3); print(fullfile(outDir, [svName '_marmoset_populationZpZc.pdf']),'-dpdf', '-fillpage') 

%% TEST -- plot my calculated Zp Zc
 
figure;
s=0;
for iCell = [1 2 3]
    s=s+1;
    subplot(3,3,s+3)
        for im = 1:2
            scatter(Zc(im,iCell), Zp(im,iCell))
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        subtitle(['cell ' num2str(iCell)])
        plotZcZpBorders
        set(gca,'TickDir','out'); axis square
end



%% TEST -- test pattern and component predictions

figure;
for iCell =1:3
    subplot(1,3,iCell)
        for im = 1:2
            polarplot([x_rad x_rad(1)], [avg_resp_plaid(iCell,:,im) avg_resp_plaid(iCell,1,im)],'k')
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat(iCell,:) avg_resp_grat(iCell,1)],'y','LineWidth',2)
        polarplot([x_rad x_rad(1)], [component(iCell,:) component(iCell,1)],'b')
        polarplot([x_rad x_rad(1)], [pattern(iCell,:) pattern(iCell,1)],'r') 
end 