close all;clear all; clc;
dataset = 'exp_list_ket_1wk_tjw'; %experiment list to pick files from
eval(dataset); %load dataset
mouse_mat = unique({expt.mouse});
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\ket\';
nmouse = length(mouse_mat);
kDiff_ket = [];
kDiff_sal = [];
maxDiff_ket = [];
maxDiff_sal = [];
oriDiff_ket = [];
oriDiff_sal = [];
oriVar_ket = [];
oriVar_sal = [];
n_sal = 0;
n_ket = 0;
x = defaultPlotColors;
figure;
for i = 1:nmouse
    load(fullfile(fnout,[cell2mat(mouse_mat(i)) '_matchedData.mat']))
    ind_ket = find(strcmp(drug_list,'ket'),1,'first');
    ind_sal = find(strcmp(drug_list,'sal'),1,'first');
    ind_rec = find(strcmp(drug_list,'rec'),1,'first');
    ind_base = find(strcmp(drug_list,'base'),1,'first');
    if ~isempty(ind_ket)
        kDiff_temp = (alld_k{ind_ket}-alld_k{ind_base})./(alld_k{ind_ket}+alld_k{ind_base});
        kDiff_ket = [kDiff_ket kDiff_temp];
        maxDiff_temp = (alld_max{ind_ket}-alld_max{ind_base})./(alld_max{ind_ket}+alld_max{ind_base});
        maxDiff_ket = [maxDiff_ket maxDiff_temp];
        oriDiff_temp = d1_comp_prefori{ind_ket-1};
        oriDiff_ket = [oriDiff_ket oriDiff_temp];
        oriVar_temp = alld_score_prefori_within{ind_ket};
        oriVar_ket = [oriVar_ket oriVar_temp];
        n_ket = n_ket+1;
        subplot(2,2,1)
        h = cdfplot(kDiff_temp);
        set(h, 'Color',x(1,:))
        xlabel('Change in Tuning width (k)')
        ylabel('Fraction of cells')
        title([])
        hold on
        subplot(2,2,2)
        h = cdfplot(maxDiff_temp);
        set(h, 'Color',x(1,:))
        xlabel('Change in max dF/F')
        ylabel('Fraction of cells')
        title([])
        hold on
        subplot(2,2,3)
        h = cdfplot(oriDiff_temp);
        set(h, 'Color',x(1,:))
        xlabel('Change in pref ori (vs Base)')
        ylabel('Fraction of cells')
        xlim([0 90])
        title([])
        hold on
        subplot(2,2,4)
        h = cdfplot(oriVar_temp);
        set(h, 'Color',x(1,:))
        title([])
        xlabel('Change in pref ori (within sess)')
        ylabel('Fraction of cells')
        xlim([0 90])
        hold on
        fprintf([cell2mat(mouse_mat(i)) '- ketamine \n'])
    end
    if ~isempty(ind_sal)
        kDiff_temp = (alld_k{ind_sal}-alld_k{ind_base})./(alld_k{ind_sal}+alld_k{ind_base});
        kDiff_sal = [kDiff_sal kDiff_temp];
        maxDiff_temp = (alld_max{ind_sal}-alld_max{ind_base})./(alld_max{ind_sal}+alld_max{ind_base});
        maxDiff_sal = [maxDiff_sal maxDiff_temp];
        oriDiff_temp = d1_comp_prefori{ind_sal-1};
        oriDiff_sal = [oriDiff_sal oriDiff_temp];
        oriVar_temp = alld_score_prefori_within{ind_sal};
        oriVar_sal = [oriVar_sal oriVar_temp];
        n_sal = n_sal+1;
        subplot(2,2,1)
        h = cdfplot(kDiff_temp);
        set(h, 'Color','k')
        xlabel('Change in Tuning width (k)')
        ylabel('Fraction of cells')
        title([])
        hold on
        subplot(2,2,2)
        h = cdfplot(maxDiff_temp);
        set(h, 'Color','k')
        hold on
        xlabel('Change in max dF/F')
        ylabel('Fraction of cells')
        title([])
        subplot(2,2,3)
        h = cdfplot(oriDiff_temp);
        set(h, 'Color','k')
        xlabel('Change in pref ori (vs Base)')
        ylabel('Fraction of cells')
        xlim([0 90])
        title([])
        hold on
        subplot(2,2,4)
        h = cdfplot(oriVar_temp);
        set(h, 'Color','k')
        title([])
        xlabel('Change in pref ori (within sess)')
        ylabel('Fraction of cells')
        xlim([0 90])
        hold on
        fprintf([cell2mat(mouse_mat(i)) '- saline \n'])
    end
end
sgtitle(['Ketamine- n=' num2str(n_ket) '; Saline- n=' num2str(n_sal)])
print(fullfile(fnout,'KetamineSummary_eaMouse.pdf'),'-dpdf','-bestfit')

figure;
subplot(2,2,1)
cdfplot(kDiff_ket)
hold on
h = cdfplot(kDiff_sal);
set(h, 'Color','k')
legend('Ket','Sal','location','southeast')
xlabel('Change in Tuning width (k)')
ylabel('Fraction of cells')
title([])
subplot(2,2,2)
cdfplot(maxDiff_ket)
hold on
h = cdfplot(maxDiff_sal);
set(h, 'Color','k')
xlabel('Change in max dF/F')
ylabel('Fraction of cells')
title([])
subplot(2,2,3)
cdfplot(oriDiff_ket)
hold on
h = cdfplot(oriDiff_sal);
set(h, 'Color','k')
xlabel('Change in pref ori (vs Base)')
ylabel('Fraction of cells')
xlim([0 90])
title([])
subplot(2,2,4)
cdfplot(oriVar_ket)
hold on
h = cdfplot(oriVar_sal);
set(h, 'Color','k')
title([])
xlabel('Change in pref ori (within sess)')
ylabel('Fraction of cells')
xlim([0 90])
sgtitle(['Ketamine- n=' num2str(n_ket) '; Saline- n=' num2str(n_sal)])
print(fullfile(fnout,'KetamineSummary.pdf'),'-dpdf','-bestfit')
