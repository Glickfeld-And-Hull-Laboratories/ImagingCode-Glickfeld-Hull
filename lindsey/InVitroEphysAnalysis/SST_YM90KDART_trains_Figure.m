close all
clear all

fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\InVitroRecordings\';
fn_out_summary = fullfile(fn_out,'SST_YM90KDART_Trains');

load(fullfile(fn_out_summary,'SST_YM90KDART_Trains_Summary.mat'))

sz = size(data_norm2base_all);
figure;
for i = 2
    for ii = 1:2
        errorbar(1:10,mean(data_norm2base_all(match(i,ii),:,:),3),std(data_norm2base_all(match(i,ii),:,:),[],3)./sqrt(sz(3)-1))
        hold on
    end
    anova_mat = reshape(permute(data_norm2base_all([match(i,1) match(i,2)],:,:), [2 3 1]),[sz(2) sz(3)*2])';
    p_n2B = anova2(anova_mat,sz(3),'off');
    title([num2str(freq_list(match(i,1))) ' Hz; pStim = ' num2str(chop(p_n2B(1),2)) ' pDart = ' num2str(chop(p_n2B(2),2)) ' pInt = ' num2str(chop(p_n2B(3),2))])
    ylim([0 6])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
legend({'Control','YM90K-DART'})
ylabel('Norm to baseline EPSC amplitude')
print(fullfile(fn_out_summary,'Figure','normToBaseline_50Hz.pdf'),'-dpdf','-fillpage')


figure;
for i = 2
    errorbar(1:10,mean(data_norm2con_all(i,:,:),3),std(data_norm2con_all(i,:,:),[],3)./sqrt(sz(3)-1))
    hold on
    ylim([0 1])
    xlim([0 11])
    ylabel('Norm EPSC')
    xlabel('Stimulus #')
    p_n2C = anova1(squeeze(data_norm2con_all(i,:,:))',[],'off');
    title([num2str(freq_list(match(i,1))) ' Hz; p = ' num2str(chop(p_n2C,2))])
    hold on
    plot(1:10,squeeze(data_norm2con_all(i,:,:)),'k')
end
ylabel('Norm to control EPSC amplitude')
xlabel('Stimulus #')
print(fullfile(fn_out_summary,'Figure','normToControl_50Hz.pdf'),'-dpdf','-fillpage')

exp = 2;
ds = 'SST_YM90KDART_trains_ExptList';
eval(ds);
date = expt(exp).date;
abfdate = expt(exp).abfdate;
firstFile = expt(exp).firstFile;
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\InVitroRecordings\';
if firstFile<10
    firstFile_str = ['_000' num2str(firstFile)];
else
    firstFile_str = ['_00' num2str(firstFile)];
end
fn_out_use = fullfile(fn_out,date,[abfdate firstFile_str]);
load(fullfile(fn_out_use,[abfdate firstFile_str '_dataSub.mat']))

figure; 
subplot(2,1,1)
plot(data(:,2)-mean(data(2500:2700,2),1))
hold on
plot(data(:,4)-mean(data(2500:2700,4),1))
xlim([2600 5000])
ylim([-900 100])
subplot(2,1,2)
plot(-(data(:,2)-mean(data(2500:2700,2),1))./data_resp(2,1))
hold on
plot(-(data(:,4)-mean(data(2500:2700,4),1))./data_resp(4,1))
xlim([2600 5000])
ylim([-7 1])
sgtitle([abfdate firstFile_str])

print(fullfile(fn_out_summary,'Figure','exampleTCs_50Hz_exp2.pdf'),'-dpdf','-fillpage')