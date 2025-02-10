
%% experimenting with time course of PDS origination
%cells of interest for g17 -- 130 (pat270) and 170 (pat180&270)

clc; clear all; close all;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase']);
svName = 'randDirFourPhase_CrossOri';
bnsz = 10; %spike binsize = 10 ms

doPlot=0;
doGratPlot=0;

expts = strvcat('g01', 'g06', 'g12', 'g17'); %Four acute penetrations in V1, 1 marmoset
nexp = length(expts);

iexp = 2';
    load(fullfile(base, 'Data\fromNicholas\CrossOri_randDirFourPhase_V1_marmoset', ([expts(iexp,:) '.mat'])))
    
    nCells = size(resp,1);
    nDirs = size(resp,2);
    nPhas = size(resp,3);
    nTrials = size(resp,5);
    nStim = (size(resp,3)+1)*size(resp,2);

    resp_cell = resp(:,:,:,:,:,21:120);
    base_cell = resp(:,:,:,:,:,1:20);
    avg_resp_dir(:,:,:,:,1) = mean(sum(resp_cell,6)/100,5); % Average across time and trial

start=1;
for ibin = 1:100
    avg_resp_dir(:,:,:,:,ibin) = mean(sum(resp_cell(:,:,:,:,:,start:ibin),6)/100,5);
    start=start+1;

    int = nDirs;
        
        component = avg_resp_dir(:,:,1,1,ibin)+circshift(avg_resp_dir(:,:,1,1,ibin),-120./int,2);
        pattern = circshift(avg_resp_dir(:,:,1,1,ibin),-60./int,2);
        
        comp_corr = zeros(nPhas,nCells,100);
        patt_corr = zeros(nPhas,nCells,100);
        comp_patt_corr = zeros(nPhas,nCells,100);
        plaid_corr = zeros(1,nCells,100);
        plaid_corr1 = zeros(1,nCells,100);
        plaid_corr2 = zeros(1,nCells,100);
        plaid_corr3 = zeros(1,nCells,100);
        plaid_corr4 = zeros(1,nCells,100);
        plaid_corr5 = zeros(1,nCells,100);
        plaid_corr6 = zeros(1,nCells,100);
        
        for iCell = 1:nCells
            for ip = 1:nPhas
                comp_corr(ip,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),component(iCell,:)));
                patt_corr(ip,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,ip,2,1),pattern(iCell,:)));
                comp_patt_corr(ip,iCell,ibin) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
            end
            plaid_corr1(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,2,2,ibin)));
            plaid_corr2(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,3,2,ibin)));
            plaid_corr3(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,4,2,ibin)));
            plaid_corr4(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,3,2,ibin)));
            plaid_corr5(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,4,2,ibin)));
            plaid_corr6(1,iCell,ibin) = triu2vec(corrcoef(avg_resp_dir(iCell,:,3,2,1),avg_resp_dir(iCell,:,4,2,ibin)));
            plaid_corr(1,iCell,ibin) = (plaid_corr1(1,iCell,ibin)+plaid_corr2(1,iCell,ibin)+plaid_corr3(1,iCell,ibin)+plaid_corr4(1,iCell,ibin)+plaid_corr5(1,iCell,ibin)+plaid_corr6(1,iCell,ibin))/6;
        end
        
        Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
        Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
        Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nDirs-3));
        Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nDirs-3));
    
        % Rp_all(i,:,:) = Rp;
        % Rc_all(i,:,:) = Rc;
        % Zp_all(i,:,:) = Zp;
        % Zc_all(i,:,:) = Zc;
end

%%

% Rp(Rp==0)= NaN;
Rp_time = squeeze(mean(Rp,1,"omitnan"));

figure;
for ic = 93
    subplot(4,5,1)
        plot(1:100,Rp_time(ic,:))
end
























%% this was for when I was looking at individual cells
% cell 130 at plaid 3
% cell 170 at 3 and 4
% cell 14 at 3
% cell 41 at 4
% cell 178 at 3

figure;
subplot(3,4,5)
    plot(1:3,Rp_all(:,4,130))   
    ylabel('Rp'); ylim([-0.3 0.7])
subplot(3,4,6)
    plot(1:3,Rp_all(:,3,14))
    ylabel('Rp'); ylim([-0.3 0.7])
subplot(3,4,7)
    plot(1:3,Rp_all(:,4,41))
    ylabel('Rp'); ylim([-0.3 0.7])
subplot(3,4,8)
    plot(1:3,Rp_all(:,3,178))
    ylabel('Rp'); ylim([-0.3 0.7])



pattInd130 = (Rp_all(:,4,130)).^2-(Rc_all(:,4,130)).^2;
pattInd14 = (Rp_all(:,3,14)).^2-(Rc_all(:,3,14)).^2;
pattInd41 = (Rp_all(:,4,41)).^2-(Rc_all(:,4,41)).^2;
pattInd178 = (Rp_all(:,3,178)).^2-(Rc_all(:,3,178)).^2;

subplot(3,4,1)
    plot(1:3,pattInd130)
    ylabel('Rp^2 - Rc^2')
    ylim([-0.5 0.5])
    subtitle('cell 130')
subplot(3,4,2)
    plot(1:3,pattInd14)
    ylim([-0.5 0.5])
    ylabel('Rp^2 - Rc^2')
    subtitle('cell 14')
subplot(3,4,3)
    plot(1:3,pattInd41)
    ylim([-0.5 0.5])
    ylabel('Rp^2 - Rc^2')
    subtitle('cell 41')
subplot(3,4,4)
    plot(1:3,pattInd178)
    ylim([-0.5 0.5])
    ylabel('Rp^2 - Rc^2')
    subtitle('cell 178')


PCI130 = Zp_all(:,4,130)-Zc_all(:,4,130);
PCI14 = Zp_all(:,3,14)-Zc_all(:,3,14);
PCI41 = Zp_all(:,4,41)-Zc_all(:,4,41);
PCI178 = Zp_all(:,3,178)-Zc_all(:,3,178);


subplot(3,4,9)
    plot(1:3,PCI130)
    xlabel('0-100ms, 100-200ms, 200-300ms')
    ylabel('Zp-Zc'); ylim([-3 3])
subplot(3,4,10)
    plot(1:3,PCI14)
    xlabel('0-100ms, 100-200ms, 200-300ms')
    ylabel('Zp-Zc'); ylim([-3 3])
subplot(3,4,11)
    plot(1:3,PCI41)
    xlabel('0-100ms, 100-200ms, 200-300ms')
    ylabel('Zp-Zc'); ylim([-3 3])
subplot(3,4,12)
    plot(1:3,PCI178)
    xlabel('0-100ms, 100-200ms, 200-300ms')
    ylabel('Zp-Zc'); ylim([-3 3])
sgtitle('g17')


print(fullfile(outDir, [svName '_timecoursePDSresponse.pdf']),'-dpdf', '-fillpage') 


