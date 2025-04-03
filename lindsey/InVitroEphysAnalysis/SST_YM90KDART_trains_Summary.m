close all
clear all

fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\InVitroRecordings\';
ds = 'SST_YM90KDART_trains_ExptList';
fn_out_summary = fullfile(fn_out,'SST_YM90KDART_Trains');
eval(ds);

nexpt = size(expt,2);
data_norm2base_all = [];
data_norm2con_all = [];

for iexp = 1:nexpt
    date = expt(iexp).date;
    abfdate = expt(iexp).abfdate;
    firstFile = expt(iexp).firstFile;
    
    if firstFile<10
        firstFile_str = ['_000' num2str(firstFile)];
    else
        firstFile_str = ['_00' num2str(firstFile)];
    end
    fn_out_use = fullfile(fn_out,date,[abfdate firstFile_str]);

    load(fullfile(fn_out_use,[abfdate firstFile_str '_dataSub.mat']))

    data_norm2base = data_sub./data_sub(:,1);
    data_norm2con = zeros(size(data_sub,1)/2,size(data_sub,2));
    for i = 1:size(match,1)
        data_norm2con(i,:) = data_sub(match(i,2),:)./data_sub(match(i,1),:);
    end

    data_norm2base_all = cat(3,data_norm2base_all,data_norm2base);
    data_norm2con_all = cat(3,data_norm2con_all,data_norm2con);
end

p_n2B = zeros(3,3);
figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 1:2
        errorbar(1:10,mean(data_norm2base_all(match(i,ii),:,:),3),std(data_norm2base_all(match(i,ii),:,:),[],3)./sqrt(nexpt-1))
        hold on
    end
    sz = size(data_norm2base_all);
    anova_mat = reshape(permute(data_norm2base_all([match(i,1) match(i,2)],:,:), [2 3 1]),[sz(2) sz(3)*2])';
    p_n2B(i,:) = anova2(anova_mat,sz(3),'off');
    title([num2str(freq_list(match(i,1))) ' Hz; pStim = ' num2str(chop(p_n2B(i,1),2)) ' pDart = ' num2str(chop(p_n2B(i,2),2)) ' pInt = ' num2str(chop(p_n2B(i,3),2))])
    ylim([0 10])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
subplot(3,1,1)
legend({'Control','YM90K-DART'})
sgtitle('Norm to baseline EPSC amplitude')
print(fullfile(fn_out_summary,'normToBaseline.pdf'),'-dpdf','-fillpage')

colmat = ['k','c'];
figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 1:2
        plot(1:10,squeeze(data_norm2base_all(match(i,ii),:,:)),colmat(ii))
        hold on
    end
    title([num2str(freq_list(match(i,1))) ' Hz'])
    ylim([0 10])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
subplot(3,1,1)
legend({'Control','YM90K-DART'})
sgtitle('Norm to baseline EPSC amplitude')
print(fullfile(fn_out_summary,'normToBaseline_allCells.pdf'),'-dpdf','-fillpage')

p_n2C = zeros(1,3);
figure;
for i = 1:3
    subplot(3,1,i)
    errorbar(1:10,mean(data_norm2con_all(i,:,:),3),std(data_norm2con_all(i,:,:),[],3)./sqrt(nexpt-1))
    hold on
    ylim([0 1])
    xlim([0 11])
    ylabel('Norm EPSC')
    xlabel('Stimulus #')
    p_n2C(i) = anova1(squeeze(data_norm2con_all(i,:,:))',[],'off');
    title([num2str(freq_list(match(i,1))) ' Hz; p = ' num2str(chop(p_n2C(i),2))])
end
sgtitle('Norm to control EPSC amplitude')
print(fullfile(fn_out_summary,'normToControl.pdf'),'-dpdf','-fillpage')

figure;
for i = 1:3
    subplot(3,1,i)
    plot(1:10,squeeze(data_norm2con_all(i,:,:)),'k')
    hold on
    ylim([0 1])
    xlim([0 11])
    ylabel('Norm EPSC')
    xlabel('Stimulus #')
    title([num2str(freq_list(match(i,1))) ' Hz'])
end
sgtitle('Norm to control EPSC amplitude')
print(fullfile(fn_out_summary,'normToControl_allCells.pdf'),'-dpdf','-fillpage')
