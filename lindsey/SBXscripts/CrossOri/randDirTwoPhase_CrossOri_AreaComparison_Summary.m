clear all
close all
clc
doRedChannel = 0;
ds = 'CrossOriRandDirTwoPhase_ExptList';
area_list = strvcat('V1','LM');%,'AL','RL','PM');

str = {'hiSF','lowSF'};
a =2;

narea = length(area_list);
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
leg_str = cell(1,narea);
Zc_all_all = [];
Zp_all_all = [];
ZcZp_all_all = [];
resp_ind_all_all = [];
totCells = 0;
area_ind = [];
pca_area = [];
maskPhas = [0 90];
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir,['randDirTwoPhase_Summary_' str{a} '_' area_list(iA,:) '.mat']))
    areaSummary(iA).name = area_list(iA,:);
    areaSummary(iA).mice = unique(mouse_list,'rows');
    areaSummary(iA).nmice = size(unique(mouse_list,'rows'),1);
    areaSummary(iA).nexp = size(mouse_list,1);
    areaSummary(iA).totCells = length(Zc_all);
    areaSummary(iA).respCells = length(resp_ind_all);
    figure(1)
    ind = resp_ind_all;
    leg_str{1,iA} = [area_list(iA,:) '- ' num2str(length(ind)) ' cells'];
    for ip = 1:2
        subplot(2,3,1+((ip-1).*3))
        cdfplot(Zc_all(ip,ind))
        hold on
        subplot(2,3,2+((ip-1).*3))
        cdfplot(Zp_all(ip,ind))
        hold on
        subplot(2,3,3+((ip-1).*3))
        errorbar(mean(Zc_all(ip,ind),2),mean(Zp_all(ip,ind),2),std(Zc_all(ip,ind),[],2)./sqrt(length(ind)),std(Zc_all(ip,ind),[],2)./sqrt(length(ind)),std(Zp_all(ip,ind),[],2)./sqrt(length(ind)),std(Zp_all(ip,ind),[],2)./sqrt(length(ind)),'o')
        hold on
    end
    
    
    figure(2)
    for ip = 1:2
    subplot(2,narea,iA+((ip-1).*narea))
    scatter(Zc_all(ip,resp_ind_all),Zp_all(ip,resp_ind_all),'ok')
    hold on
    Zc_use = intersect(resp_ind_all,intersect(find(Zc_all(ip,:)>1.28),find(Zc_all(ip,:)-Zp_all(ip,:)>1.28)));
    scatter(Zc_all(ip,Zc_use),Zp_all(ip,Zc_use),'ob')
    Zp_use = intersect(resp_ind_all,intersect(find(Zp_all(ip,:)>1.28),find(Zp_all(ip,:)-Zc_all(ip,:)>1.28)));
    scatter(Zc_all(ip,Zp_use),Zp_all(ip,Zp_use),'or')
    xlim([-4 8])
    ylim([-4 8])
    xlabel('Zc')
    ylabel('Zp')
    title([area_list(iA,:) '- ' num2str(maskPhas(ip)) ' deg'])
    axis square
    plotZcZpBorders
    end
    
end
figure(1)
movegui('center')
subplot(2,3,1)
legend(leg_str{1,:},'location','southeast')
xlabel('Zc')
xlim([-2 6])
title('Phase = 0')
subplot(2,3,2)
xlabel('Zp')
xlim([-2 6])
title('')
subplot(2,3,3)
xlabel('Zc')
ylabel('Zp')
xlim([-0.5 2])
ylim([-0.5 2])
subplot(2,3,4)
xlabel('Zc')
xlim([-2 6])
title('Phase = 90')
subplot(2,3,5)
xlabel('Zp')
xlim([-2 6])
title('')
subplot(2,3,6)
xlabel('Zc')
ylabel('Zp')
xlim([-0.5 2])
ylim([-0.5 2])
print(fullfile(summaryDir, ['randDirTwoPhase_' str{a} '_allArea_summary.pdf']),'-dpdf', '-fillpage') 

figure(2)
movegui('center')
print(fullfile(summaryDir, ['randDirTwoPhase_' str{a} 'allArea_summary_ZcZp_scatters.pdf']),'-dpdf', '-fillpage')
