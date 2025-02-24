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

figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 1:2
        errorbar(1:10,mean(data_norm2base_all(match(i,ii),:,:),3),std(data_norm2base_all(match(i,ii),:,:),[],3)./sqrt(nexpt-1))
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
print(fullfile(fn_out_summary,'normToBaseline.pdf'),'-dpdf','-fillpage')

figure;
errorbar(1:10,mean(data_norm2con_all,3),std(data_norm2con_all,[],3)./sqrt(nexpt-1))
hold on
ylim([0 1])
xlim([0 11])
ylabel('Norm EPSC')
xlabel('Stimulus #')
legend(num2str(freq_list(1:3)'))
sgtitle('Norm to control EPSC amplitude')
print(fullfile(fn_out_summary,'normToControl.pdf'),'-dpdf','-fillpage')
