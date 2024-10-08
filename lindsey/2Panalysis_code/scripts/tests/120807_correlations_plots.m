fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);
col = strvcat('c', 'k', 'r');
vec0 = zeros(2,7);
vec0(1,:) = [.5 1 2 4 8 16 32];
vec0(2,:) = [.01 .02 .04 .08 .16 .32 .64];    
vec1 = interp2(vec0');
TF_vec = zeros(8,1);
SF_vec = zeros(8,1);
TF_vec(1,1)= vec1(2,1);
TF_vec(2,1)= 1.0001;
TF_vec(3:6,1) = vec1(4:2:10,1);
TF_vec(7,1) = 14.999;
TF_vec(8,1) = vec1(12,1);
SF_vec(1,1)= vec1(2,3);
SF_vec(2,1)= .020001;
SF_vec(3:6,1) = vec1(4:2:10,3);
SF_vec(7,1) = .3199;
SF_vec(8,1) = vec1(12,3);
edges = log2(TF_vec./flipud(SF_vec));
plot_edges = [edges(2) mean(edges(2:3)) mean(edges(3:4)) mean(edges(4:5)) mean(edges(5:6)) mean(edges(6:7)) edges(7)];
real_avg = zeros(7,3);
real_sem = zeros(7,3);
shuf_avg = zeros(7,3);
shuf_sem = zeros(7,3);
sub_avg = zeros(7,3);
sub_sem = zeros(7,3);
blank_avg = zeros(7,3);
blank_sem = zeros(7,3);
figure;
for iArea = 1:3
    [n_speed bin_speed] = histc(log2(Goodfits(iArea).speed), edges);
    for ibin = 1:7
        ind_speed= find(bin_speed==ibin);
        real_avg(ibin,iArea) = nanmean(Corr(iArea).real(ind_speed));
        real_sem(ibin,iArea) = nanstd(Corr(iArea).real(ind_speed),[],1)./sqrt(length(ind_speed));
        shuf_avg(ibin,iArea) = nanmean(Corr(iArea).shuf(ind_speed));
        shuf_sem(ibin,iArea) = nanstd(Corr(iArea).shuf(ind_speed),[],1)./sqrt(length(ind_speed));
        sub_avg(ibin,iArea) = nanmean(Corr(iArea).sub(ind_speed));
        sub_sem(ibin,iArea) = nanstd(Corr(iArea).sub(ind_speed),[],1)./sqrt(length(ind_speed));
        blank_avg(ibin,iArea) = nanmean(Corr(iArea).blank(ind_speed));
        blank_sem(ibin,iArea) = nanstd(Corr(iArea).blank(ind_speed),[],1)./sqrt(length(ind_speed));
    end

    subplot(2,2,1)
    errorbar(plot_edges, real_avg(:,iArea), real_sem(:,iArea), col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    ylabel('corr coeff');
    xlim([log2(2) log2(1000)])
    ylim([0 .3])
    title('Real Corr')
    subplot(2,2,2)
    errorbar(plot_edges, shuf_avg(:,iArea), shuf_sem(:,iArea), col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    ylabel('corr coeff');
    xlim([log2(2) log2(1000)])
    ylim([0 .3])
    title('Shuffled Corr')
    subplot(2,2,3)
    errorbar(plot_edges, sub_avg(:,iArea), sub_sem(:,iArea), col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    ylabel('corr coeff');
    xlim([log2(2) log2(1000)])
    ylim([0 .15])
    title('Real-Shuffled Corr')
    subplot(2,2,4)
    errorbar(plot_edges, blank_avg(:,iArea), blank_sem(:,iArea), col(iArea,:));
    hold on
    xlabel('log2(speed)');
    axis square
    ylabel('corr coeff');
    xlim([log2(2) log2(1000)])
    ylim([0 .3])
    title('Blank Corr')
end

base = '\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120807'
fn_out = fullfile(base, [matrix '_' num2str(P) 'P_' inj '_speed_vs_corr_histograms.ps']);
        print(gcf, '-depsc2', fn_out);
fn_out = fullfile(base,  [matrix '_' num2str(P) 'P_' inj '_speed_vs_corr_histograms.pdf']);
        print(gcf, '-dpdf', fn_out);
