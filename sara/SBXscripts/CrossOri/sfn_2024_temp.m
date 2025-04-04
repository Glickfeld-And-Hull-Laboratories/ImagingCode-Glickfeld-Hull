
    outDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase\WholeCellfromLiam']);
    svName = 'randDirFourPhase_CrossOri';


    figure(1);
    print(fullfile(outDir, [svName '_cell1_traces.pdf']),'-dpdf', '-fillpage') 
    figure(4);
    set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_cell1_linearity.pdf']),'-dpdf', '-fillpage') 


    figure(2);
    print(fullfile(outDir, [svName '_cell3_traces.pdf']),'-dpdf', '-fillpage') 
    figure(7);
    set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_cell3_linearity.pdf']),'-dpdf', '-fillpage') 


    figure(3);
    print(fullfile(outDir, [svName '_cell4_traces.pdf']),'-dpdf', '-fillpage') 
    figure(8);
    set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_cell4_linearity.pdf']),'-dpdf', '-fillpage') 


    figure(5);
    axis square; set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_zMeanDSI.pdf']),'-dpdf', '-fillpage') 
    figure(6);
    axis square; set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_zModDSI.pdf']),'-dpdf', '-fillpage') 

    figure(9);
    axis square; set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_zMeanOSI.pdf']),'-dpdf', '-fillpage') 
    figure(10);
    axis square; set(gca,'TickDir','out');
    print(fullfile(outDir, [svName '_zModOSI.pdf']),'-dpdf', '-fillpage') 

%%



    outDir = ([base '\Figures\CrossOri_meetings\fromNicholas']);
    svName = 'randDirFourPhase_CrossOri';

    figure(1);
    print(fullfile(outDir, [svName '_FeedforwardModel_exmod3polar.pdf']),'-dpdf', '-fillpage')

    figure(2);
    print(fullfile(outDir, [svName '_FeedforwardModel_exstable3polar.pdf']),'-dpdf', '-fillpage')





    %%


    figure;
% 
% avg_resp_plddir = circshift(avg_resp_dir(:,:,1,1,1),-60./int,2);

x=[-150:30:180];
x_rad = deg2rad(x);

for ic = 1:5
    subplot(5,4,ic)
            polarplot([x_rad x_rad(1)], [resps3stable(ic,:) resps3stable(ic,1)])
            hold on
    subplot(5,4,ic+8)
            polarplot([x_rad x_rad(1)], [resps3mod(ic,:) resps3stable(ic,1)])
            hold on
 
end

   





    % Determine pattern and component direction selectivity
    int = 12;
    nDirs=12;
    
    avg_resp_dir = resps3mod;

    component = avg_resp_dir(5,:)+circshift(avg_resp_dir(5,:),-120./int,2);
    pattern = circshift(avg_resp_dir(5,:),-60./int,2);

    nCells=1;
    nPhas=4;
    
    comp_corr = zeros(nPhas,1);
    patt_corr = zeros(nPhas,1);
    comp_patt_corr = zeros(nPhas,nCells);
    plaid_corr = zeros(1,nCells);
    plaid_corr1 = zeros(1,nCells);
    plaid_corr2 = zeros(1,nCells);
    plaid_corr3 = zeros(1,nCells);
    plaid_corr4 = zeros(1,nCells);
    plaid_corr5 = zeros(1,nCells);
    plaid_corr6 = zeros(1,nCells);
    

        for ip = 1:nPhas
            comp_corr(ip) = triu2vec(corrcoef(avg_resp_dir(ip,:),component));
            patt_corr(ip) = triu2vec(corrcoef(avg_resp_dir(ip,:),pattern));
            comp_patt_corr(ip) = triu2vec(corrcoef(component,pattern));
        end

    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nDirs-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nDirs-3));


    



PCI = (Zp-Zc);

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;

    subplot(5,4,start)
        scatter(phase,PCI,8,'filled');
        hold on
        [b_hat_all, amp_hat_all, per_hat_all,pha_hat_all,sse_all,R_square_all] = sinefit_PCI(deg2rad(phase),PCI);
        yfit_all(:,1) = b_hat_all+amp_hat_all.*(sin(2*pi*deg2rad(phase_range)./per_hat_all+ 2.*pi/pha_hat_all));
        plot(phase_range, yfit_all(:,1),'k'); 
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
  
        base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'modelV1_PCImodulation_mod.pdf']),'-dpdf', '-fillpage')  



%% Liam whole cell examples --> fits

cell1_z = [-2.1916 -2.8041 -2.0681 -2.8555];
cell3_z = [-1.6978 -1.5192 -1.3404 -1.2407];
cell4_z = [-2.8892 -2.9306 -4.1904 -3.1166];


PCI =cell4_z;

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;

    subplot(5,4,start)
        scatter(phase,PCI,8,'filled');
        hold on
        [b_hat_all, amp_hat_all, per_hat_all,pha_hat_all,sse_all,R_square_all] = sinefit_PCI(deg2rad(phase),PCI);
        yfit_all(:,1) = b_hat_all+amp_hat_all.*(sin(2*pi*deg2rad(phase_range)./per_hat_all+ 2.*pi/pha_hat_all));
        plot(phase_range, yfit_all(:,1),'k'); 
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
  
        base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'wholecellV1_PCImodulation_cell4.pdf']),'-dpdf', '-fillpage')  

%% Liam population plots

r_val = allrs(1:10);
mZ = mean(allzpzcs,2);
maZ = max(allzpzcs,[],2);

figure(1);
subplot(4,4,1)
    histogram(r_val,10)
    hold on
    xlim([-0.3 1])
     set(gca,'TickDir','out')
    subtitle('r values')
subplot(4,4,2)
    scatter(mZ(1:10),r_val,10,'filled')
    hold on
    xlabel('mean zp-zc')
    ylabel('r')
     set(gca,'TickDir','out')

subplot(4,4,5)
    scatter(maZ(1:10),r_val,10,'filled')
    hold on
    xlabel('max zp-zc')
    ylabel('r')
     set(gca,'TickDir','out')

%%


PCI = allzpzcs(1:10,:)';

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;

yfit_all = zeros(10,360);

for ic = 1:10
    % subplot(5,4,iCell)
        % scatter(phase,PCI(:,ic),8,'filled');
        % hold on
        [b_hat_all(ic,1), amp_hat_all(ic,1), per_hat_all(ic,1),pha_hat_all(ic,1),sse_all(ic,1),R_square_all(ic,1)] = sinefit_PCI(deg2rad(phase),PCI(:,ic));
        yfit_all(ic,:) = b_hat_all(ic,1)+amp_hat_all(ic,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(ic,1) + 2.*pi/pha_hat_all(ic,1)));
            % plot(phase_range, yfit_all(ic,:,1),'k:');
            % subtitle(['cell ' num2str(ic) ', Rsq ' num2str(R_square_all(ic),'%.2f'), ', SSE ' num2str(sse_all(ic),'%.2f')])
        % ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        % xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square

end



pattpeak = max(yfit_all,[],2);
figure(1);
subplot(3,4,3)
    scatter(-1*amp_hat_all,pattpeak,7,'filled')
    hold on
        xlim([-2.5 0])
        ylim([-5 5])
        set(gca,'TickDir','out'); axis square
subplot(4,4,6)
    scatter(amp_hat_all(1:10),r_val,10,'filled')
    hold on
    xlabel('amp')
    ylabel('r')
     set(gca,'TickDir','out')
    print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'wholecellV1_population_likeNicholas.pdf']),'-dpdf', '-fillpage')  


    %%

    figure(3);
    set(gca,'TickDir','out');
    print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'wholecellV1_cell1_linearity.pdf']),'-dpdf', '-fillpage')  

    figure(4);
    set(gca,'TickDir','out');
    print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'wholecellV1_cell3_linearity.pdf']),'-dpdf', '-fillpage')  

    figure(5);
    set(gca,'TickDir','out');
    print(fullfile(base, 'Figures\CrossOri_meetings\fromNicholas', [ 'wholecellV1_cell4_linearity.pdf']),'-dpdf', '-fillpage')  
