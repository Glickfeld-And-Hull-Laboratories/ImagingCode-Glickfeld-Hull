
%make a subset of normalized difference for the SST cells only, then make
% find how many are facilitated or suppressed by more than 1 std from
% baseline
%reselect norm_diff)red to be only the "red all" subset
norm_diff_red = norm_diff(:,:,red_all);

facil_red=norm_diff_red(:,:,:)>=1;
supp_red=norm_diff_red(:,:,:)<=-1;

N=length(red_all);
facil_table_stat = sum(facil_red(1,:,:),3)/N;
supp_table_stat = sum(supp_red(1,:,:),3)/N;
facil_table_loc = sum(facil_red(2,:,:),3)/N;
supp_table_loc = sum(supp_red(2,:,:),3)/N;

% Number of bootstrap samples
num_bootstraps = 1000;
% bootsrapping 95% CI for suppression, stat
observed_data =squeeze(supp_red(1,:,:)); %only taking the stationary values
% Bootstrap resampling
bootstrap_samples = zeros(num_bootstraps, size(observed_data,1),size(observed_data,2));
for i = 1:num_bootstraps
    for iCon = 1:nCon
        bootstrap_samples(i, iCon,:) = datasample(observed_data(iCon,:), size(observed_data,2));
    end
end
% Calculate fraction suppressed for each bootstrap sample
bootstrap_statistics = sum(bootstrap_samples, 3)/N;
% Calculate the 95% confidence interval
confidence_interval_supp_stat = prctile(bootstrap_statistics, [2.5, 97.5])

clear observed_data bootstrap_statistics bootstrap_samples

% bootsrapping 9% CI for facilitation, stat
observed_data =squeeze(facil_red(1,:,:)); %only taking the stationary values
% Bootstrap resampling
bootstrap_samples = zeros(num_bootstraps, size(observed_data,1),size(observed_data,2));
for i = 1:num_bootstraps
    for iCon = 1:nCon
        bootstrap_samples(i, iCon,:) = datasample(observed_data(iCon,:), size(observed_data,2));
    end
end
% Calculate fraction suppressed for each bootstrap sample
bootstrap_statistics = sum(bootstrap_samples, 3)/N;
% Calculate the 95% confidence interval
confidence_interval_facil_stat = prctile(bootstrap_statistics, [2.5, 97.5])



% bootsrapping 9% CI for suppression, running
observed_data =squeeze(supp_red(2,:,:)); %only taking the stationary values
% Bootstrap resampling
bootstrap_samples = zeros(num_bootstraps, size(observed_data,1),size(observed_data,2));
for i = 1:num_bootstraps
    for iCon = 1:nCon
        bootstrap_samples(i, iCon,:) = datasample(observed_data(iCon,:), size(observed_data,2));
    end
end
% Calculate fraction suppressed for each bootstrap sample
bootstrap_statistics = sum(bootstrap_samples, 3)/N;
% Calculate the 95% confidence interval
confidence_interval_supp_loc = prctile(bootstrap_statistics, [2.5, 97.5])

clear observed_data bootstrap_statistics bootstrap_samples

% bootsrapping 9% CI for facilitation, running
observed_data =squeeze(facil_red(2,:,:)); %only taking the stationary values
% Bootstrap resampling
bootstrap_samples = zeros(num_bootstraps, size(observed_data,1),size(observed_data,2));
for i = 1:num_bootstraps
    for iCon = 1:nCon
        bootstrap_samples(i, iCon,:) = datasample(observed_data(iCon,:), size(observed_data,2));
    end
end
% Calculate fraction suppressed for each bootstrap sample
bootstrap_statistics = sum(bootstrap_samples, 3)/N;
% Calculate the 95% confidence interval
confidence_interval_facil_loc = prctile(bootstrap_statistics, [2.5, 97.5])



figure;
subplot(1,2,1)
b=bar([1,2,3],[supp_table_stat; supp_table_loc],'grouped','FaceColor',"#00ffff",'EdgeColor', [1 1 1])
b(1).FaceColor="#70D0F6"
b(2).FaceColor="#0C8ABB"
hold on
error_lengths = [supp_table_stat - confidence_interval_supp_stat(1,:); confidence_interval_supp_stat(2,:) - supp_table_stat];
%errorbar([.85,1.85,2.85], supp_table_stat, error_lengths(1, :), error_lengths(2, :), 'k.', 'LineWidth', 1, 'CapSize', 3);
error_lengths = [supp_table_loc - confidence_interval_supp_loc(1,:); confidence_interval_supp_loc(2,:) - supp_table_loc];
%errorbar([1.15,2.15,3.15], supp_table_loc, error_lengths(1, :), error_lengths(2, :), 'k.', 'LineWidth', 1, 'CapSize', 3);
xticklabels({'25','50','100'})
ylim([0 .6])
title('Suppressed')
ylabel(["Fraction SST cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

subplot(1,2,2)
b=bar([1,2,3],[facil_table_stat; facil_table_loc],'FaceColor',"#a329cc",'EdgeColor', [1 1 1])
b(1).FaceColor="#C983B1"
b(2).FaceColor="#883367"
hold on
error_lengths = [facil_table_stat - confidence_interval_facil_stat(1,:); confidence_interval_facil_stat(2,:) - facil_table_stat];
%errorbar([.85,1.85,2.85], facil_table_stat, error_lengths(1, :), error_lengths(2, :), 'k.', 'LineWidth', 1, 'CapSize', 3);
error_lengths = [facil_table_loc - confidence_interval_facil_loc(1,:); confidence_interval_facil_loc(2,:) - facil_table_loc];
%errorbar([1.15,2.15,3.15], facil_table_loc, error_lengths(1, :), error_lengths(2, :), 'k.', 'LineWidth', 1, 'CapSize', 3);
xticklabels({'25','50','100'})
ylim([0 .6])
title('Facilitated')
ylabel(["Fraction HTP+ cells"]) 
xlabel(["Contrast(%)"])
set(gca,'TickDir','out')
box off

x0=5;
y0=5;
width=3;
height=1.5;
set(gcf,'units','inches','position',[x0,y0,width,height])