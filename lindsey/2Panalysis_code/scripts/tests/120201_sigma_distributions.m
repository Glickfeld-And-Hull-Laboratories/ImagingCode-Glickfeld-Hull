start= 1;
figure;
for iArea = 1:3
subplot(3,2,start)
scatter(Goodfits(iArea).sigma_SF, Goodfits(iArea).dF)
xlabel('sigma SF')
ylabel('dF/F')
subplot(3,2,start+1)
scatter(Goodfits(iArea).sigma_TF, Goodfits(iArea).dF)
xlabel('sigma TF')
ylabel('dF/F')
start = start+2;
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_sigma_vs_dF_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);



figure;
edges_sigma = [0:.5:3.5];
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).sigma_SF, edges_sigma);
    for ibin = 1:length(edges_sigma)
        ind_ibin = find(bin==ibin);
        df_avg(ibin,1) = mean(Goodfits(iArea).dF(ind_ibin),1);
        df_avg(ibin,2) = std(Goodfits(iArea).dF(ind_ibin),1)./sqrt(length(ind_ibin));
    end
    [n bin] = histc(Goodfits(iArea).sigma_TF, edges_sigma);
    for ibin = 1:length(edges_sigma)
        ind_ibin = find(bin==ibin);
        df_avg(ibin,3) = mean(Goodfits(iArea).dF(ind_ibin),1);
        df_avg(ibin,4) = std(Goodfits(iArea).dF(ind_ibin),1)./sqrt(length(ind_ibin));
    end
    subplot(2,2, 1);
    errorbar(edges_sigma,df_avg(:,1),df_avg(:,2), col(iArea,:));
    xlabel('sigma SF')
    ylabel('dF/F')
    ylim([0 .8])
    hold on
    subplot(2,2, 2);
    errorbar(edges_sigma,df_avg(:,3),df_avg(:,4), col(iArea,:));
    xlabel('sigma TF')
    ylabel('dF/F')
    ylim([0 .8])
    hold on
end

edges_sigma = [0:.5:3.5];
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).sigma_SF, edges_sigma);
    for ibin = 1:length(edges_sigma)
        ind_ibin = find(bin==ibin);
        df_avg(ibin,1) = mean(Goodfits(iArea).speed(ind_ibin),1);
        df_avg(ibin,2) = std(Goodfits(iArea).speed(ind_ibin),1)./sqrt(length(ind_ibin));
    end
    [n bin] = histc(Goodfits(iArea).sigma_TF, edges_sigma);
    for ibin = 1:length(edges_sigma)
        ind_ibin = find(bin==ibin);
        df_avg(ibin,3) = mean(Goodfits(iArea).speed(ind_ibin),1);
        df_avg(ibin,4) = std(Goodfits(iArea).speed(ind_ibin),1)./sqrt(length(ind_ibin));
    end
    subplot(2,2, 3);
    errorbar(edges_sigma,df_avg(:,1),df_avg(:,2), col(iArea,:));
    xlabel('sigma SF')
    ylabel('speed')
    ylim([0 500])
    hold on
    subplot(2,2, 4);
    errorbar(edges_sigma,df_avg(:,3),df_avg(:,4), col(iArea,:));
    xlabel('sigma TF')
    ylabel('speed')
    ylim([0 500])
    hold on
end
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_all_areas_sigma_vs_dFandspeed.pdf']);
        print(gcf, '-dpdf', fn_out);        