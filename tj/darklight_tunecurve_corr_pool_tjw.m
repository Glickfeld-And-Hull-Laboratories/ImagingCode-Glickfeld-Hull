%%
%clear everything
clear all
clear all global
clc
close all

%%
%find folders to load and experiment info

fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P'; %folder to load files from
newfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\multi_day_drift'; %folder to save files to
realfnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\tj\Analysis\Analysis\2P\darklight\pooled\drift'; %folder to save files to
dataset = 'exp_list_darklight_actual_tjw'; %experiment list to pick files from
eval(dataset); %load dataset

%%

abcd_corrs_all = [];
corr_7d_all = [];
corr_14d_all = [];
corr_21d_all = [];

list = [1 8 15 22 26 30 34];
figure;

for iexp = list
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    img_area = expt(iexp).img_loc{1};
    img_layer = expt(iexp).img_loc{2};
    corrs = load(fullfile(newfnout, [mouse '_' img_area '_' img_layer '_tunecurve_corr']));
    abcd_corrs = corrs.abcd_corrs;    
    scatter(1,[abcd_corrs(1,2), abcd_corrs(2,3), abcd_corrs(3,4)])
    hold on
    scatter(2, [abcd_corrs(1,3), abcd_corrs(2,4)])
    scatter(3, abcd_corrs(1,4))
    corr_7d = [abcd_corrs(1,2); abcd_corrs(2,3); abcd_corrs(3,4)];
    corr_7d_all = [corr_7d_all; corr_7d]; 
    corr_14d = [abcd_corrs(1,3); abcd_corrs(2,4)];
    corr_14d_all = [corr_14d_all; corr_14d];
    corr_21d = [abcd_corrs(1,4)];
    corr_21d_all = [corr_21d_all; corr_21d];
end

xlim([0 4])
ylim([0 1])
xlabel('Days Apart')
ylabel('Correlation')
xticklabels({'','','7','','14','','21','',''})

plot(my_fit)

%%

mean_7d = mean(corr_7d_all);
stderr_7d = std(corr_7d_all)/sqrt(length(corr_7d_all));
mean_14d = mean(corr_14d_all);
stderr_14d = std(corr_14d_all)/sqrt(length(corr_14d_all));
mean_21d = mean(corr_21d_all);
stderr_21d = std(corr_21d_all)/sqrt(length(corr_21d_all));


%%

%%

%%

corr_labels_7d = cell(length(corr_7d_all),1);
corr_labels_7d(:) = {7};
corr_labels_14d = cell(length(corr_14d_all),1);
corr_labels_14d(:) = {14};
corr_labels_21d = cell(length(corr_21d_all),1);
corr_labels_21d(:) = {21};
corr_labels_all = cell2mat([corr_labels_7d; corr_labels_14d; corr_labels_21d]);

corr_table = table(corr_labels_all,[corr_7d_all; corr_14d_all; corr_21d_all], 'VariableNames', {'days_between', 'population_corr'});

%%

my_equation = 'a*exp(1)^(b*x)+c'

my_fit = fit(corr_table.days_between, corr_table.population_corr, my_equation)

figure; plot(my_fit)


%%
fake_x = [1;1;1;1;2;2;2;2;2;2;3;3;3;3;4;4;5;5;6;7;7;10;10;11;12;13;14;14;14;14;15;15;16;17;18;21;21;21;21;22;22;23];

fake_equation = 'a+b*exp(1)^(c*x)'

my_fit_fake = fit(corr_table.days_between, corr_table.population_corr, fake_equation)

fake_fit = fit(fake_x, corr_table.population_corr, my_equation)

figure; plot(fake_fit)

%%
figure;
subplot(1,2,1)
plot(my_fit)
subplot(1,2,2)
plot(my_fit_fake)